import re
import os
import logging
import pandas as pd
import numpy as np
from urllib.parse import urlparse
import json

from .dprint import dprint
from .varstash import Var
from .error import *

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)




####################################################################################################
####################################################################################################
class AmpliconMatrix:

    def __init__(self, upa):
        '''
        Instance variables created during init:
        * upa
        * name
        * obj
        '''
        self.upa = upa

        self._get_obj()


    def _get_obj(self):
        logging.info('Loading object info for AmpliconMatrix `%s`' % self.upa)

        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
        })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']

        # comment run_dir with AmpliconMatrix name
        dprint('touch %s' % os.path.join(Var.run_dir, '#' + self.name), run='cli')


    def to_OTU_table(self, tax_l, flpth=None):
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
            index=tax_l,#self.row_attr_map.get_tax_l(row_ids), 
            columns=col_ids
            )
        df.index.name = "taxonomy"
        df['OTU_Id'] = row_ids # add OTU_Id for identification purposes (?)
        df = df[['OTU_Id'] + col_ids] # reorder

        if flpth is not None:
            df.to_csv(flpth, sep='\t')

        return df

    def _map_id2attr_ids(self, id2attr, axis='row'):
        '''
        Parameters
        ----------
        id2attr - AmpliconMatrix row_ids to attribute you want to give AttributeMapping

        Behavior
        --------
        Swap out ids in id2attr so they end up mapping AttributeMapping ids to attributes
        '''
        if f'{axis}_attributemapping_ref' not in self.obj:
            raise Exception(
                'Trying to map AmpliconMatrix %s_ids to %s AttributeMapping ids '
                "when AmpliconMatrix doesn't have %s AttributeMapping"
                % (axis, axis, axis)
            )
        elif f'{axis}_mapping' not in self.obj:
            msg = (
                'Dude this object has a %s_attributemapping_ref '
                'and needs a %s_mapping. Letting it slide for now.'
                % (axis, axis)
            )
            logging.warning(msg)
            var.warnings.append(msg)
            return id2attr

        id2attr = {
            self.obj[f'{axis}_mapping'][id]: attr
            for id, attr in id2attr.items()
        }

        return id2attr

    def save(self, name=None):
        info = Var.dfu.save_objects({
            'id': Var.params['workspace_id'],
            "objects": [{
                "type": "KBaseMatrices.AmpliconMatrix",
                "data": self.obj,
                "name": name if name else self.name, # null case name="" or name=None
                "extra_provenance_input_refs": [self.upa]
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new




####################################################################################################
####################################################################################################
class AttributeMapping:

    def __init__(self, upa, amp_mat):
        '''
        Instance variables created at init time:
        * upa
        * amp_mat
        * name
        * obj
        '''
        self.upa = upa
        self.amp_mat = amp_mat
        self._get_obj()


    def _get_obj(self):
        logging.info('Loading object info for AttributeMapping `%s`' % self.upa)

        obj = Var.dfu.get_objects({
            'object_refs': ['%s;%s' % (amp_mat.upa, self.upa)]
        })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']



    def get_tax_ind_attribute(self, tax_field):
        '''
        Get index and name of taxonomy attribute
        Whether it's from the user-entered or default regexes
        '''
        if len(tax_field) == 1:
            tax_field = tax_field[0]


        else:
            attribute_l = [d['attribute'] for d in self.obj['attributes']]

            for regex in regex_l:
                for i, attribute in enumerate(attribute_l):
                    if re.search(regex, attribute.lower()) is not None:
                        return i, attribute
            return None, None


    def get_tax_l(self, tax_ind, id_l):
        tax_l = []

        for id in id_l:
            tax_l.append(
                self.obj['instances'][id][tax_ind]
            )

        return tax_l



    def map_update_attribute(self, ind: int, id2attr: dict, map_ids_first=True):
        '''
        First map ids of `id2attr` to AttributeMapping ids using amp_mat's row_mapping
        Update attribute at index `ind` using mapping `id2attr`
        '''
        if map_ids_first is True:
            id2attr = self.amp_mat._map_id2attr_ids(id2attr, axis='row')

        for id, attr in id2attr.items():
            self.obj['instances'][id][ind] = attr


    def add_attribute_slot_warn(self, attribute, source) -> int:
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

