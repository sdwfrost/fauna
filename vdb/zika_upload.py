import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import parser

class zika_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.grouping_optional_fields = ['lineage']

    def fix_name(self, name):
        tmp_name = name
        tmp_name = tmp_name.replace('Human', '').replace('human', '').replace('H.sapiens_tc', '').replace('Hsapiens_tc', '').replace('H.sapiens-tc', '').replace('Homo_sapiens', '').replace('Homo sapiens', '').replace('Hsapiens', '').replace('H.sapiens', '')
        tmp_name = tmp_name.replace('_Asian', '').replace('_Asia', '').replace('_asian', '').replace('_asia', '')
        tmp_name = tmp_name.replace('Zika_virus', '').replace('Zika virus', '').replace('Zika', '').replace('ZIKV', '')
        tmp_name = tmp_name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('//', '/').replace('__', '_').replace('.', '').replace(',', '')
        tmp_name = re.sub('^/', '', tmp_name)
        try:
            tmp_name = 'V' + str(int(tmp_name))
        except:
            pass
        return tmp_name

    def format_schema(self):
        pass

    def upload_documents(self, **kwargs):
        '''
        Insert viruses into collection
        '''
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        db_relaxed_strains = self.relaxed_strains()
        # Faster way to upload documents, downloads all database documents locally and looks for precense of strain in database
        db_viruses = list(r.db(self.database).table(self.table).run())
        db_strain_to_viruses = {db_v['strain']: db_v for db_v in db_viruses}
        update_viruses = {}
        upload_viruses = {}
        for virus in self.viruses:
            # determine the corresponding database strain name based on relaxed db and virus strain
            db_strain = virus['strain']
            if self.relax_name(virus['strain']) in db_relaxed_strains:
                db_strain = db_relaxed_strains[self.relax_name(virus['strain'])]
            if db_strain in db_strain_to_viruses.keys():  # virus already in database
                update_viruses[db_strain] = virus
            elif db_strain in upload_viruses.keys():  # virus already to be uploaded, need to check for updates to sequence information
                upload_v = upload_viruses[db_strain]
                self.update_document_sequence(upload_v, virus, **kwargs)  # add new sequeunce information to upload_v
            else:  # new virus that needs to be uploaded
                upload_viruses[virus['strain']] = virus
        print("Inserting ", len(upload_viruses), "viruses into database", self.table)
        try:
            r.table(self.table).insert(upload_viruses.values()).run()
        except:
            raise Exception("Couldn't insert new viruses into database")
        print("Checking for updates to ", len(update_viruses), "viruses in database", self.table)
        updated = []
        for db_strain, v in update_viruses.items():  # determine if virus has new information
            document = db_strain_to_viruses[db_strain]
            updated_base = self.update_base(document, v, v.keys(), v['strain'], **kwargs)
            if updated_base:
                document['timestamp'] = v['timestamp']
                updated.append(document)
        try:
            r.table(self.table).insert(updated, conflict="replace").run()
        except:
            raise Exception("Couldn't update viruses already in database")

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = {0:'accession', 2:'strain', 4:'date', 6:'country'}
    # 0        1          2      3  4          5     6
    #>KU501216|Zika_virus|103344|NA|2015_12_01|Human|Guatemala
    setattr(args, 'fasta_fields', fasta_fields)
    connVDB = zika_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
