import os, json, datetime, sys, re
import rethinkdb as r
from Bio import SeqIO
import numpy as np
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io

def get_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-db', '--database', default='vdb', help="database to download from")
    parser.add_argument('--rethink_host', default=None, help="rethink host url")
    parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
    parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")
    parser.add_argument('-v', '--virus', help="virus name")
    parser.add_argument('--ftype', default='fasta', help="output file format, default \"fasta\", other options are \"json\" and \"tsv\"")
    parser.add_argument('--fstem', default=None, help="default output file name is \"VirusName_Year_Month_Date\"")
    parser.add_argument('--path', default='data', help="path to dump output files to")
    parser.add_argument('--fasta_fields', default=['strain', 'virus', 'accession', 'date', 'region', 'country', 'division', 'location', 'source', 'locus', 'authors'], help="fasta fields for output fasta")

    parser.add_argument('--public_only', default=False, action="store_true", help="include to subset public sequences")
    parser.add_argument('--select', nargs='+', type=str, default=[], help="Select specific fields ie \'--select field1:value1 field2:value1,value2\'")
    parser.add_argument('--present', nargs='+', type=str, default=[], help="Select specific fields to be non-null ie \'--present field1 field2\'")
    parser.add_argument('--interval', nargs='+', type=str, default=[], help="Select interval of values for fields \'--interval field1:value1,value2 field2:value1,value2\'")
    parser.add_argument('--relaxed_interval', default=False, action="store_true", help="Relaxed comparison to date interval, 2016-XX-XX in 2016-01-01 - 2016-03-01")

    parser.add_argument('--pick_longest', default=False, action="store_true",  help ="For duplicate strains, only includes the longest sequence for each locus")
    return parser

class download(object):
    def __init__(self, database, virus, **kwargs):
        '''
        parser for virus, fasta fields, output file names, output file format path, interval
        '''
        self.virus = virus.lower()
        self.viruses_table = virus + "_viruses"
        self.sequences_table = virus + "_sequences"
        self.database = database.lower()
        if self.database not in ['vdb', 'test_vdb', 'test']:
            raise Exception("Can't download from this database: " + self.database)
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        self.rethink_io.check_table_exists(self.database, self.viruses_table)
        self.rethink_io.check_table_exists(self.database, self.sequences_table)

    def count_documents(self, table):
        '''
        return integer count of number of documents in table
        '''
        return r.db(self.database).table(table).count().run()

    def download(self, output=True, **kwargs):
        '''
        download documents from table
        '''
        import time
        start_time = time.time()
        select, present, public_only = self.parse_subset_arguments(**kwargs)
        sequence_count = r.table(self.sequences_table).count().run()
        print(sequence_count, "sequences in table:", self.sequences_table)
        virus_count = r.table(self.viruses_table).count().run()
        print(virus_count, "viruses in table:", self.viruses_table)
        print("Downloading documents from the sequence table: " + self.sequences_table, " and virus table: ", self.viruses_table)
        sequences = self.rethinkdb_download(self.sequences_table, self.viruses_table, present, select, public_only, **kwargs)
        sequences = self.subset(sequences, **kwargs)
        sequences = self.resolve_duplicates(sequences, **kwargs)
        if output:
            self.output(sequences, **kwargs)
        print("--- %s minutes to download ---" % ((time.time() - start_time)/60))

    def parse_subset_arguments(self, select=[], present=[], public_only=False, **kwargs):
        '''
        Parse arguments needed for subsetting on the server side
        Return a tuple of the arguments
        '''

        selections = self.parse_select_argument(select)
        return selections, present, public_only

    def parse_select_argument(self, groupings=[]):
        '''
        :arg groupings like country:brazil,argentina
        parse the 'select' parameter to determine which field name to filter and for what values
        :return: [(grouping name, [group values])] ex. [(country, [brazil, argentina)]
        '''
        selections = []
        if groupings is not None and len(groupings)>0:
            for group in groupings:
                result = group.split(':')
                selections.append((result[0].lower(), result[1].lower().split(',')))
        return selections

    def rethinkdb_download(self, sequence_table, virus_table, presents=[], selections=[], public=False, index='strain', **kwargs):
        '''
        Default command merges documents from the sequence table and virus table
        Chain rethinkdb filter and has_fields commands to the default command
        Return documents from the database that are left after filtering
        '''
        # take each sequence and merge with corresponding virus document
        command = r.table(sequence_table).merge(lambda sequence: r.table(virus_table).get(sequence[index]))
        # Check if documents have fields in `present`
        if len(presents)>0:
            print("Only downloading documents with fields: ", presents)
            command = command.has_fields(r.args(presents))
        # Check if documents have the correct value for certain fields
        if len(selections)>0:
            for sel in selections:
                print("Only downloading documents with field ", sel[0], "equal to", sel[1])
                command = command.filter(lambda doc: r.expr(sel[1]).contains(doc[sel[0]]))
        # Check if documents are public
        if public:
            print("Only downloading public sequences")
            command.filter({'public': True})
        sequences = list(command.run())
        return list(sequences)

    def subset(self, sequences, interval=[], **kwargs):
        '''
        Filter out sequences that are not in the date interval specified
        '''
        intervals = self.parse_select_argument(interval)
        if len(intervals) > 0:
            for sel in intervals:
                if sel[0] in ['collection_date', 'date', 'submission_date'] or sel[0] == 'submission_date':
                    older_date, newer_date = self.check_date_format(sel[1][0], sel[1][1])
                    sequences = filter(lambda doc: self.in_date_interval(doc[sel[0]], older_date, newer_date, **kwargs), sequences)
                    print('Removed documents that were not in the interval specified (' + ' - '.join([older_date, newer_date]) + ') for field \'' + sel[0] + '\', remaining documents: ' + str(len(sequences)))
        return sequences

    def check_date_format(self, older_date, newer_date):
        one_sided_symbols = ['', 'XXXX-XX-XX']
        if newer_date in one_sided_symbols:
            newer_date = str(datetime.datetime.strftime(datetime.datetime.utcnow(),'%Y-%m-%d'))
        if older_date in one_sided_symbols:
            older_date = '0000-00-00'
        if older_date > newer_date:
            raise Exception("Date interval must list the earlier date first")
        if not re.match(r'\d\d\d\d-(\d\d)-(\d\d)$', older_date) or not re.match(r'\d\d\d\d-(\d\d)-(\d\d)$', newer_date):
            raise Exception("Date interval must be in YYYY-MM-DD format with all values defined", older_date, newer_date)
        return(older_date.upper(), newer_date.upper())

    def in_date_interval(self, virus_date, older_date, newer_date, **kwargs):
        '''
        :return: true if the date is in the interval older_date:newer_date, otherwise False
        '''
        virus_gt_old_date = self.date_greater(virus_date.split('-'), older_date.split('-'), **kwargs)
        new_date_gt_virus = self.date_greater(newer_date.split('-'), virus_date.split('-'), **kwargs)
        return virus_gt_old_date and new_date_gt_virus

    def date_greater(self, greater_date, comparison_date, relaxed_interval=False, **kwargs):
        '''
        Dates in YYYY-MM-DD format
        :return: true if greater_date >= comparison_date
        '''
        # compare year
        if greater_date[0] < comparison_date[0]:
            return False
        elif greater_date[0] == comparison_date[0]:
            # compare month
            if greater_date[1] == 'XX' or comparison_date[1] == 'XX':
                return relaxed_interval
            elif greater_date[1] < comparison_date[1]:
                return False
            elif greater_date[1] == comparison_date[1]:
                # compare day
                if greater_date[2] == 'XX' or comparison_date[2] == 'XX':
                    return relaxed_interval
                elif greater_date[2] < comparison_date[2]:
                    return False
        return True

    def resolve_duplicates(self, sequences, pick_longest=True, **kwargs):
        strain_locus_to_doc = {doc['strain']+doc['locus']: doc for doc in sequences}
        if pick_longest:
            print("Resolving duplicate strains and locus by picking the longest sequence")
            for doc in sequences:
                if doc['strain']+doc['locus'] in strain_locus_to_doc:
                    if self.longer_sequence(doc['sequence'], strain_locus_to_doc[doc['strain']+doc['locus']]):
                        strain_locus_to_doc[doc['strain']+doc['locus']] = doc
                else:
                    strain_locus_to_doc[doc['strain']] = doc
        return list(strain_locus_to_doc.values())

    def longer_sequence(self, long_seq, short_seq):
        '''
        :return: true if long_seq is longer than short_seq
        '''
        return len(long_seq) > len(short_seq)

    def write_json(self, data, fname, indent=1):
        '''
        writes as list of viruses (dictionaries)
        '''
        try:
            handle = open(fname, 'w')
        except:
            print("Couldn't open output file")
            print(fname)
            raise FileNotFoundError
        else:
            json.dump(data, handle, indent=indent)
            handle.close()

    def write_fasta(self, viruses, fname, sep='|', fasta_fields=['strain', 'virus', 'accession'], **kwargs):
        try:
            handle = open(fname, 'w')
        except IOError:
            pass
        else:
            for virus in viruses:
            	if 'sequence' in virus:
            	    if virus['sequence']:
                        fields = [str(virus[field]) if (field in virus and virus[field] is not None) else '?'
                                  for field in fasta_fields]
                        handle.write(">"+sep.join(fields)+'\n')
                        handle.write(virus['sequence'] + "\n")
            handle.close()

    def write_tsv(self, viruses, fname, sep='\t', fasta_fields=['strain', 'virus', 'accession'], **kwargs):
        try:
            handle = open(fname, 'w')
        except IOError:
            pass
        else:
            handle.write(sep.join(fasta_fields)+'\n')
            for virus in viruses:
                fields = [str(virus[field]) if (field in virus and virus[field] is not None) else '?'
                          for field in fasta_fields]
                handle.write(sep.join(fields)+'\n')
            handle.close()

    def output(self, documents, path, fstem, ftype, **kwargs):
        fname = path + '/' + fstem + '.' + ftype
        print("Outputing", len(documents), "documents to ", fname)
        if ftype == 'json':
            self.write_json(documents,fname)
        elif ftype == 'fasta':
            self.write_fasta(documents, fname, **kwargs)
        elif ftype == 'tsv':
            self.write_tsv(documents, fname, **kwargs)
        else:
            raise Exception("Can't output to that file type, only json, fasta or tsv allowed")
        print("Wrote to " + fname)

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y_%m_%d'))
    if args.fstem is None:
        args.fstem = args.virus + '_' + current_date
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connVDB = download(**args.__dict__)
    connVDB.download(**args.__dict__)