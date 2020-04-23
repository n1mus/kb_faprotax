import os
import logging
import pandas as pd
import numpy as np

from .dprint import dprint
from .varstash import Var




pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)




class AttributeMapping:

    def __init__(self, upa):
        self.upa = upa
        self._get_obj()


    def _get_obj(self):
        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']



    def parse_faprotax_traits(self, groups2records_table_dense_flpth) -> dict:
        '''
        Input: filepath for groups2records_dense.tsv
        Output: dict map from taxonomy to predicted functions
        '''
        g2r_df = pd.read_csv(groups2records_table_dense_flpth, sep='\t', comment='#')
        g2r_df = g2r_df.fillna('').drop_duplicates().set_index('record')

        r2g_d = g2r_df.to_dict(orient='index')
        r2g_d = {record: r2g_d[record]['group'].replace(',', ':') for record in r2g_d}

        return r2g_d


    def update_attribute(self, attr_to_attr_d, attribute_lookup, attribute_add):
        '''
        Add attribute for `attribute_add`, looking up by `attribute_lookup`
        '''
        for i, attr_d in enumerate(self.obj['attributes']):
            if attr_d['attribute'] == attribute_lookup:
                lkp_ind = i
            elif attr_d['attribute'] == attribute_add:
                add_ind = i

        for _, attr_l in self.obj['instances'].items():
            attr_l[add_ind] = attr_to_attr_d.get(attr_l[lkp_ind], '')

        self.obj['attributes'][add_ind]['source'] = 'kb_faprotax'


    def add_attribute_slot(self, attribute):
        
        # check if already exists
        for attr_d in self.obj['attributes']:
            if attr_d['attribute'] == attribute:
                msg = 'Adding attribute slot `%s` to AttributeMapping with name `%s`, ' % (attribute, self.name) + \
                      'but that attribute already exists in object'
                logging.warning(msg)
                Var.warnings.append(msg)
                return

        # append slot to `attributes`
        self.obj['attributes'].append({
            'attribute': attribute,
            })

        # append slots to `instances` 
        for _, attr_l in self.obj['instances'].items():
            attr_l.append('')


    def save(self):
        
        info = Var.dfu.save_objects({
            'id': Var.params['workspace_id'],
            "objects": [{
                "type": "KBaseExperiments.AttributeMapping",
                "data": self.obj,
                "name": self.name,
                "extra_provenance_input_refs": [self.upa]
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        dprint('upa_new', run=locals())

        return upa_new







class AmpliconSet:

    def __init__(self, upa):
        self.upa = upa
        self._get_obj()


    def _get_obj(self):
        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.name = obj['data'][0]['info'][1]
        self.amp_mat_upa = obj['data'][0]['data']['amplicon_matrix_ref']
        self.obj = obj['data'][0]['data']


    def get_amplicon_matrix_upa(self):
        return self.amp_mat_upa


    def update_amplicon_matrix_ref(self, amp_mat_upa_new):
        self.obj['amplicon_matrix_ref'] = amp_mat_upa_new


    def save(self, name=None):

        info = Var.dfu.save_objects({
            'id': Var.params['workspace_id'],
            "objects": [{
                "type": "KBaseExperiments.AmpliconSet",
                "data": self.obj,
                "name": name if name else self.name,
                "extra_provenance_input_refs": [self.upa]
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new




class AmpliconMatrix:

    def __init__(self, upa, amp_set: AmpliconSet, test=False):
        self.upa = upa
        self.amp_set = amp_set
        self.test = test

        self._get_obj()
        self._to_OTU_table()


    def _get_obj(self):
        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']
        self.row_attrmap_upa = self.obj['row_attributemapping_ref'] 


    def _to_OTU_table(self):

        logging.info(f"Parsing AmpliconMatrix data from object")

        data = np.array(self.obj['data']['values'], dtype=float)
        row_ids = self.obj['data']['row_ids']
        col_ids = self.obj['data']['col_ids']

        data = pd.DataFrame(
            data, 
            index=self._get_taxonomy_l(row_ids), 
            columns=col_ids
            )
        data.index.name = "taxonomy"
        data['OTU_Id'] = row_ids # TODO get rid of this
        data = data[['OTU_Id'] + col_ids]

        if self.test:
            data = data.iloc[:50]

        self.taxon_table_flpth = os.path.join(Var.sub_dir, 'taxon_table.tsv')

        data.to_csv(self.taxon_table_flpth, sep='\t')

        
    def _get_taxonomy_l(self, amplicon_id_l):
        NUM_TAX_LVL = 11

        amplicon_d = self.amp_set.obj['amplicons']
        taxonomy_l = []

        for amplicon_id in amplicon_id_l:
            taxonomy = amplicon_d[amplicon_id]['taxonomy']['lineage'][:NUM_TAX_LVL]
            taxonomy = '; '.join(taxonomy)
            taxonomy_l.append(taxonomy)

        return taxonomy_l


    def update_row_attributemapping_ref(self, row_attrmap_upa_new):
        self.obj['row_attributemapping_ref'] = row_attrmap_upa_new


    def save(self):

        info = Var.dfu.save_objects({
            'id': Var.params['workspace_id'],
            "objects": [{
                "type": "KBaseMatrices.AmpliconMatrix",
                "data": self.obj,
                "name": self.name,
                "extra_provenance_input_refs": [self.upa]
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new






