import logging
import pandas as pd
import numpy as np

from .dprint import *
from .varstash import Var




pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)



class AmpliconSet:

    def __init__(self, upa):
        self.upa = upa

        self._get_obj()



    def _get_obj(self):
        self.obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })['data'][0]['data']

        dprint(self.obj)

        self.amp_mat_upa = self.obj['amplicon_matrix_ref']


    def get_amplicon_matrix_upa(self):
        return self.amp_mat_upa




class AmpliconMatrix:

    def __init__(self, upa):
        self.upa = upa
        self._get_obj()
        self._parse_data()


    def _get_obj(self):
        self.obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })['data'][0]['data']

        dprint(self.obj)

        


    def _parse_data(self):

        logging.info(f"Parsing AmpliconMatrix data from object")

        data = np.array(self.obj['data']['values'], dtype=float)
        rows = self.obj['data']['row_ids']
        cols = self.obj['data']['col_ids']

        dprint('data', run=locals())

        data = pd.DataFrame(data, index=FakeData.tax_path_l, columns=cols) #TODO real tax paths
        data.index.name = "taxonomy"
        data['OTU_Id'] = rows
        data = data[['OTU_Id'] + cols]

        dprint('data', run=locals())

        self.taxon_table_flpth = os.path.join(Var.subdir, 'taxon_table.tsv')

        data.to_csv(self.taxon_table_flpth, sep=',')

        






class FakeData:

    tax_path_l = [
"k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Thermomonas;s__fusca",
"D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Bacillales;D_4__Bacillaceae;D_5__Bacillus;D_6__Bacillus sp. YZ5",
"Bacteria; Chlorobi; Chlorobia; Chlorobiales; Chlorobiaceae; Chlorobaculum; Chlorobaculum thiosulfatiphilum DSM249T",   
"Bacteria; Chlorobi; Chlorobia; Chlorobiales; Chlorobiaceae; Chlorobaculum",
"Bacteria; Chlorobi; Chlorobia; Chlorobiales; Chlorobiaceae; Chlorobaculum; uncultured bacterium",
"Bacteria;Proteobacteria;Alphaproteobacteria;Rhodobacterales;Rhodobacteraceae;Paracoccus;Fervidicola"            
        ]

































