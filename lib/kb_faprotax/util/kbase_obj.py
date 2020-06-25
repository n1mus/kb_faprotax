import os
import logging
import pandas as pd
import numpy as np
from urllib.parse import urlparse

from .dprint import dprint
from .varstash import Var
from .error import *
from .message import *




pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)





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


    def to_OTU_table(self, flpth):
        '''
        Write OTU table with dummy values
        '''
        with open(flpth, 'w') as f: # TODO use `to_csv`?
            f.write('OTU\tdummy_sample\n')
            for index, row in self.df.iterrows():
                f.write('%s\t1.0\n' % row['taxonomy'])
                
            



####################################################################################################
####################################################################################################
class Genome:

    def __init__(self, upa, data, metadata, oi):
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
    def best_taxonomy(self) -> tuple:
        '''
        Experimental. 
        This is beyond the scope of this program.
        It requires deciding which of various possible taxonomies is best
        '''
        # has something for taxonomy
        if self.taxonomy != None:
            return self.taxonomy, 'taxonomy'

        # get taxon

            #return , 'relation engine'

        # return what is required
        else:
            return self.domain + ';' + self.sci_name, 'domain and scientific name'

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




####################################################################################################
####################################################################################################
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







####################################################################################################
####################################################################################################
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




