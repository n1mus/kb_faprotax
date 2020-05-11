import os
import logging
import pandas as pd
import numpy as np

from .dprint import dprint
from .varstash import Var
from .error import *
from .message import *




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


    def update_attribute(self, ind: int, taxStr_2_traits_d: dict, amp_set):
        '''
        Update attribute at index `ind`
        
        taxStr -> amplicon ids -> instances to insert traits
        '''

        # taxStr -> amplicon ids

        taxStr_2_ids_d = {}

        for id_, amplicon_d in amp_set.obj['amplicons'].items():
            taxStr = '; '.join(amplicon_d['taxonomy']['lineage'])
            if taxStr in taxStr_2_ids_d:
                taxStr_2_ids_d[taxStr].append(id_)
            else:
                taxStr_2_ids_d[taxStr] = [id_]

        dprint('taxStr_2_ids_d', run=locals())
        
        # taxStr -> amplicon ids -> instances to insert traits
        
        for taxStr, traits in taxStr_2_traits_d.items():
            for id_ in taxStr_2_ids_d.get(taxStr, []):
                self.obj['instances'][id_][ind] = traits


    def add_attribute_slot(self, attribute, source) -> int:
        '''
        If attribute not already entered, add slot for it
        Return its index in the attributes/instances
        '''
        
        # check if already exists
        for ind, attr_d in enumerate(self.obj['attributes']):
            if attr_d['attribute'] == attribute: # TODO check if same source for identity
                msg = msg_overwriteAttribute % (attribute, self.name)
                logging.warning(msg)
                Var.warnings.append(msg)
                return ind

        # append slot to `attributes`
        self.obj['attributes'].append({
            'attribute': attribute,
            'source': source,
            })

        # append slots to `instances` 
        for attr_l in self.obj['instances'].values():
            attr_l.append('')

        return len(attr_l) - 1


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
        self.row_attrmap_upa = self.obj.get('row_attributemapping_ref') 


    def _to_OTU_table(self):

        logging.info(f"Parsing AmpliconMatrix data from object")

        data = np.array(self.obj['data']['values'], dtype=float)
        row_ids = self.obj['data']['row_ids']
        col_ids = self.obj['data']['col_ids']

        data = pd.DataFrame(
            data, 
            index=self._get_taxStr_l(row_ids), 
            columns=col_ids
            )
        data.index.name = "taxonomy"
        data['OTU_Id'] = row_ids # TODO get rid of this? doesn't help
        data = data[['OTU_Id'] + col_ids]

        if self.test:
            data = data.iloc[:50]

        self.taxon_table_flpth = os.path.join(Var.run_dir, 'taxon_table.tsv')

        data.to_csv(self.taxon_table_flpth, sep='\t')

        
    def _get_taxStr_l(self, id_l):
        amplicon_d_d = self.amp_set.obj['amplicons']
        taxStr_l = []

        for id_ in id_l:
            taxonomy = amplicon_d_d[id_]['taxonomy']
            try:
                taxonomy = '; '.join(taxonomy['lineage']) # TODO allow missing taxonomies?
            except KeyError:
                raise NoTaxonomyException(msg_missingTaxonomy % id_)
                
            taxStr_l.append(taxonomy)

        return taxStr_l


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






