import os
import logging
import pandas as pd
import numpy as np
from urllib.parse import urlparse
import json

from .dprint import dprint
from .varstash import Var
from .error import *
from .message import *

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)


def write_json(obj, flnm, AmpliconSet=False):
    '''
    For debugging/testing
    '''
    if 'run_dir' not in Var:
        import uuid
        Var.run_dir = os.path.join('/kb/module/work/tmp', str(uuid.uuid4()))
        os.mkdir(Var.run_dir)

    flpth = os.path.join(Var.run_dir, flnm)
    with open(flpth, 'w') as f:
        json.dump(obj, f)

    if AmpliconSet == True:
        dprint('touch %s' % os.path.join(Var.run_dir, '#' + obj['data'][0]['info'][1]), run='cli') # annotate run_dir with name



####################################################################################################
####################################################################################################
class GenomeSet:

    def __init__(self, upa):
        '''
        Instance variables created during constructor:
        * df - starts out with columns `taxonomy`, `url`
        '''

        # get GenomeSet data
        obj = Var.dfu.get_objects({
            'object_refs': [upa]
        })

        self.name = obj['data'][0]['info'][1]
        obj = obj['data'][0]['data']

        self.description = obj['description']
        
        genome_element_l = []

        # gather fields for each Genome
        for _, genome_d in obj['elements'].items():
            upa = genome_d.get('ref')
            data = genome_d.get('data')
            metadata = genome_d.get('metadata')

            genome_element_l.append((upa, data, metadata))

        num_genomes = len(genome_element_l)

        # check for duplicate upas
        upa_l = [genome_element[0] for genome_element in genome_element_l if genome_element[0] != None]
        if len(set(upa_l)) < len(upa_l):
            Var.warnings.append(msg_dupGenomes)

        # sort out which Genomes came as upas, which came as data
        ind_upa = [i for i in range(len(genome_element_l)) if genome_element_l[i][0] != None]
        ind_data = [i for i in range(len(genome_element_l)) if genome_element_l[i][2] != None]

        assert sorted(ind_upa + ind_data) == list(range(num_genomes))

        # make params of upas for goi3
        objects = [{'ref': genome_element_l[i][0]} for i in ind_upa]

        # goi3
        logging.info('Fetching Genome info and metadata')
        oi_l_temp = Var.ws.get_object_info3({'objects': objects, 'includeMetadata': 1})['infos']

        # insert placeholder `None`s into `oi_l` for Genomes that came as data, not upa
        oi_l = []
        oi_gen = (e for e in oi_l_temp)
        for i in range(num_genomes):
            if i in ind_upa:
                oi_l.append(next(oi_gen))
            elif i in ind_data:
                oi_l.append(None)
            else:
                assert False
        
        # sanity check
        try:
            next(oi_gen)
            assert False
        except StopIteration:
            pass
        except:
            assert False

        # make Genome
        genome_l = []
        for genome_element, oi in zip(genome_element_l, oi_l):
            upa, data, metadata = genome_element
            genome_l.append(Genome(upa, data, metadata, oi))

        self.genome_l = genome_l

        # make df
        self._to_table()

    def _to_table(self):
        rows = []
        column_names = ['name', 'taxonomy', 'url']
        for genome in self.genome_l:
            rows.append((genome.name, genome.taxonomy, genome.landing_url))

        df = pd.DataFrame((row for row in rows), columns=column_names)

        self.df = df


    def to_OTU_table(self, flpth, dummy_val=1.0):
        '''
        Write OTU table with dummy values
        `dummy_value` can be nonnegative number, or `'random'`
        '''
        df = self.df.copy()[['taxonomy']] 
        if dummy_val == 'random':
            df['dummy_sample'] = np.random.random(len(df)) + 0.01 # unif[0.01, 1.01)
        elif dummy_val > 0:
            df['dummy_sample'] = [dummy_val] * len(df)
        else:
            raise Exception('`dummy_val` must be gt 0')

        df.to_csv(flpth, sep='\t', float_format='%.1f', index=False)
        
        '''
        with open(flpth, 'w') as f: # TODO use `to_csv`?
            f.write('OTU\tdummy_sample\n')
            for index, row in self.df.iterrows():
                f.write('%s\t1.0\n' % row['taxonomy'])
        '''        
            



####################################################################################################
####################################################################################################
class Genome:

    def __init__(self, upa, data, metadata, oi):
        '''
        Instance variables created at init time:
        * upa
        * data
        * name
        * taxonomy
        * sci_name
        * domain
        '''
        dprint('upa', 'data', 'metadata', run=locals())

        self.upa = upa
        self.data = data

        # Genome by upa
        if upa != None and oi != None:
            self.name = oi[1]
            ver = oi[2].split('-')[-1] # depends on environment
            self.taxonomy = oi[-1].get('Taxonomy')
            self.sci_name = oi[-1]['Name']
            self.domain = oi[-1]['Domain']

        # Genome by embedded `data`
        # TODO improve
        elif data != None and oi == None:
            Var.warnings.append(\
"Genome is defined in GenomeSet rather than its own KBaseGenomes.Genome object.")
            
            self.name = 'Unnamed'
            ver = 'TODO' # depends on GenomeSet and environment
            self.taxonomy = data['taxonomy']
            self.sci_name = data['scientific_name']
            self.domain = data['domain']

        else:
            raise Exception('Either `upa` and `oi` xor `data` should be valid')

    @property
    def landing_url(self):
        if self.upa == None:
            return None
        else:
            parsed = urlparse(Var.kbase_endpoint) 
            return '%s://%s/#dataview/%s' % (parsed.scheme, parsed.netloc, self.upa) 
        

####################################################################################################
####################################################################################################
class AmpliconSet:

    def __init__(self, upa):
        '''
        Instance variables created at init time:
        * upa
        * name
        * amp_mat_upa
        * obj
        '''
        self.upa = upa
        self._get_obj()


    def _get_obj(self):
        logging.info('Loading object info for AmpliconSet `%s`' % self.upa)

        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
        })

        if Var.debug: write_json(obj, 'get_objects_AmpliconSet.json', AmpliconSet=True)

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']
        self.amp_mat_upa = self.obj['amplicon_matrix_ref']


    @staticmethod
    def _concat_tax(tax: list):
        return '; '.join(tax)

        
    def get_taxStr_l(self, id_l: list): # TODO allow missing taxonomies? 
        '''
        Given `id_l`, return corresponding `taxStr_l`
        '''
        amplicons_d = self.obj['amplicons']
        taxStr_l = []

        for id in id_l:
            taxonomy_d = amplicons_d[id]['taxonomy']
            
            if 'lineage' not in taxonomy_d: # TODO disallow `'lineage': {}`?
                raise NoTaxonomyException(msg_missingTaxonomy % id)

            taxStr = self._concat_tax(taxonomy_d['lineage'])   
            taxStr_l.append(taxStr)

        return taxStr_l


    def _get_tax2ids_d(self):
        '''
        Taxonomy is a string, ids is a list
        '''
        tax2ids_d = {}

        for id, amplicon_d in self.obj['amplicons'].items():
            tax = self._concat_tax(amplicon_d['taxonomy']['lineage'])
            if tax in tax2ids_d:
                tax2ids_d[tax].append(id)
            else:
                tax2ids_d[tax] = [id]

        return tax2ids_d


    def get_id2functions_d(self, tax2functions_d: dict): 
        '''
        AmpliconSet has id-tax mapping

        Taxonomies and functions are both strings
        `tax2functions_d` is computed from FAPROTAX output

        Goal:
        [tax2functions]-------[using-tax2ids]-------->[id2functions]
        '''
        id2functions_d = {}
        tax2ids_d = self._get_tax2ids_d()

        for tax, functions in tax2functions_d.items():
            for id in tax2ids_d[tax]:
                id2functions_d[id] = functions

        return id2functions_d



    def update_amplicon_matrix_ref(self, amp_mat_upa_new):
        self.obj['amplicon_matrix_ref'] = amp_mat_upa_new


    def save(self, name=None):

        info = Var.dfu.save_objects({
            'id': Var.params['workspace_id'],
            "objects": [{
                "type": "KBaseExperiments.AmpliconSet",
                "data": self.obj,
                "name": name if name != None else self.name,
                "extra_provenance_input_refs": [self.upa]
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new




####################################################################################################
####################################################################################################
class AmpliconMatrix:

    def __init__(self, upa, amp_set: AmpliconSet):
        '''
        Instance variables created during init:
        * upa
        * amp_set - referring AmpliconSet object
        * name
        * obj
        '''
        self.upa = upa
        self.amp_set = amp_set

        self._get_obj()


    def _get_obj(self):
        logging.info('Loading object info for AmpliconMatrix `%s`' % self.upa)

        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
        })

        if Var.debug: write_json(obj, 'get_objects_AmpliconMatrix.json')

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']

        #dprint('self.upa', 'self.obj', run=locals(), max_lines=None)


    def to_OTU_table(self, flpth=None):
        '''
        `taxonomy` is index
        `OTU_Id` is first column

        This interface is used for testing
        Return df for testing
        '''

        logging.info(f"Parsing AmpliconMatrix data from object")

        data = np.array(self.obj['data']['values'], dtype=float)
        row_ids = self.obj['data']['row_ids']
        col_ids = self.obj['data']['col_ids']

        df = pd.DataFrame(
            data, 
            index=self.amp_set.get_taxStr_l(row_ids), 
            columns=col_ids
            )
        df.index.name = "taxonomy"
        df['OTU_Id'] = row_ids # add OTU_Id for identification purposes (?)
        df = df[['OTU_Id'] + col_ids] # reorder

        if flpth != None:
            df.to_csv(flpth, sep='\t')

        return df


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




####################################################################################################
####################################################################################################
class AttributeMapping:

    def __init__(self, upa):
        '''
        Instance variables created at init time:
        * upa
        * name
        * obj
        '''
        self.upa = upa
        self._get_obj()


    def _get_obj(self):
        logging.info('Loading object info for AttributeMapping `%s`' % self.upa)

        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
        })

        if Var.debug: write_json(obj, 'get_objects_AttributeMapping.json')

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']



    def update_attribute(self, ind: int, id2attr_d: dict):
        '''
        Update attribute at index `ind` using mapping `id2attr_d`
        '''
        for id, attr in id2attr_d.items():
            self.obj['instances'][id][ind] = attr


    def add_attribute_slot(self, attribute, source) -> int:
        '''
        If attribute not already entered, add slot for it
        Return its index in the attributes/instances
        '''
        
        # check if already exists
        for ind, attr_d in enumerate(self.obj['attributes']):
            if attr_d['attribute'] == attribute and attr_d['source'] == source:
                msg = msg_overwriteAttribute % (attribute, source, self.name)
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




