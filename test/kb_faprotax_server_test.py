# -*- coding: utf-8 -*-
from configparser import ConfigParser
import os
import logging
import time
import unittest
import uuid
import pandas as pd
import numpy as np
from shutil import rmtree, copytree # TODO compress test files
import json
from unittest.mock import patch
import pandas as pd
from urllib.parse import urlparse
import re

from kb_faprotax.kb_faprotaxImpl import kb_faprotax
from kb_faprotax.kb_faprotaxServer import MethodContext
from kb_faprotax.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace

from kb_faprotax.util.error import * # exception library
from kb_faprotax.util.message import * # warning/exception messages library
from kb_faprotax.util.dprint import dprint, where_am_i
from kb_faprotax.util.varstash import Var # dict-like dot-access app globals
from kb_faprotax.util.kbase_obj import GenomeSet, Genome, AmpliconSet, AmpliconMatrix, AttributeMapping
from util.mock import * # mock business
from util.upa import * # upa library



# TODO allow toggling patching as much as possible, i.e., for unit testing 
# TODO dummy test datasets (to file) = less file io (slow?) ... leave big datasets for ?

######################################
######################################
######### TOGGLE PATCH ###############
######################################
###################################### 
do_patch = True

if do_patch:
    patch_ = patch
    patch_dict_ = patch.dict

else:
    patch_ = lambda *args, **kwargs: lambda f: f
    patch_dict_ = lambda *args, **kwargs: lambda f: f
######################################
######################################
######################################
######################################



class kb_faprotaxTest(unittest.TestCase):

####################################################################################################
####################################################################################################
####################################################################################################

    # doesn't need mocking - fast
    def test_parse_faprotax_functions(self):
        '''Test function used to parse FAPROTAX's predicted functions from its output'''
        from kb_faprotax.util.workflow import parse_faprotax_functions

        flnm = 'groups2records_dense.tsv'
        flpth = os.path.join(testData_dir, 'by_dataset_input/refseq/return/faprotax_output', flnm)
        
        r2g_d = parse_faprotax_functions(flpth)

        taxonomy = '\
cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Moraxellaceae; \
Acinetobacter; Acinetobacter calcoaceticus/baumannii complex; Acinetobacter pittii'
        functions = '\
aerobic_chemoheterotrophy,human_pathogens_septicemia,human_pathogens_pneumonia,human_pathogens_nosocomia,\
human_pathogens_all,animal_parasites_or_symbionts,aromatic_compound_degradation,chemoheterotrophy'

        self.assertTrue(r2g_d[taxonomy] == functions)


    # doesn't need mocking - fast
    def test_run_check(self):
        '''Test function that runs FAPROTAX'''
        from kb_faprotax.util.workflow import run_check

        with self.assertRaises(NonZeroReturnException) as cm:
            run_check('set -o pipefail && ;s |& tee tmp')
            self.assertTrue('`2`') in str(cm.exception) # return code 2

        with self.assertRaises(NonZeroReturnException) as cm:
            run_check('set -o pipefail && tmp |& tee tmp')
            self.assertTrue('`127`') in str(cm.exception) # return code 127

        with self.assertRaises(NonZeroReturnException) as cm:
            run_check('set -o pipefail && echo hi |& tmp')
            self.assertTrue('`127`') in str(cm.exception) # return code 127

        run_check('set -o pipefail && echo hi |& tee tmp') # run correctly


####################################################################################################
########################### GenomeSet input ########################################################
####################################################################################################

    # doesn't need mocking - fast
    def test_GenomeSet_input(self):
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': refseq,
            }
        )

    # doesn't need mocking - fast
    def test_dummy_abundance(self):
        '''
        Test different dummy abundances in FAPROTAX input OTU table
        (The GenomeSet workflow uses dummy values in the input OTU table 
        just to map Genomes to FAPROTAX functions)
        '''
        from kb_faprotax.util.workflow import run_check, parse_faprotax_functions
        from itertools import combinations

        run_dir = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(run_dir)

        gs = GenomeSet(refseq)

        dummy_vals = [1.0, 20, 'random']
        groups2records_dense_flpths = []

        # run FAPROTAX with different dummy values
        for val in dummy_vals:
            otu_table_flpth = os.path.join(run_dir, 'otu_table_' + str(val) + '.tsv')
            log_flpth = os.path.join(run_dir, 'log_' + str(val) + '.txt')
            
            gs.to_OTU_table(otu_table_flpth)
            groups2records_dense_flpth = os.path.join(run_dir, 'groups2records_dense_' + str(val) + '.tsv') 

            groups2records_dense_flpths.append(groups2records_dense_flpth)

            cmd = ' '.join([
                'set -o pipefail &&',
                Var.cmd_flpth,
                '--input_table', otu_table_flpth,
                '--input_groups_file', Var.db_flpth,
                '--out_groups2records_table_dense', groups2records_dense_flpth,
                '--row_names_are_in_column', 'taxonomy',
                '--verbose',
                '|& tee',
            ])

            with open(os.path.join(run_dir, 'cmd_' + str(val) + '.txt'), 'w') as f:
                f.write(cmd)

            run_check(cmd)

        # parse FAPROTAX output for each dummy value
        r2g_d_l = []
        for groups2records_dense_flpth in groups2records_dense_flpths:
            r2g_d = parse_faprotax_functions(groups2records_dense_flpth)
            r2g_d_l.append(r2g_d)

        # check all output is equal
        for combo in combinations(r2g_d_l, 2):
            self.assertTrue(combo[0] == combo[1])

 
    # TODO superficially test ctor
    def test_GenomeSet_methods(self):
        pass
       

    def test_Genome_methods(self):
        '''
        Genome methods are quite simple, 
        they are tested in other places
        '''
        pass

    # doesn't need mocking - fast
    def test_dup_GenomeSet(self):
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': refseq_dup,
            })

        self.assertTrue(msg_dupGenomes in Var.warnings)



####################################################################################################
########################## AmpliconSet input #######################################################
####################################################################################################
    '''
    Because of the control flow of testing,
    there are two main strategies here for patching
    (1) patch *Util/Workspace class 
        with constructor-like 
        in kb_faprotaxImpl,
        good for integration tests
    (2) patch.dict Var.dfu or Var.ws 
        with function side-effect 
        in the module it is used,
        good for unit tests
    '''

    # TODO mock dfu, run_check
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *args, **kwargs: get_mock_kbr())
    def test_AmpliconSet_no_row_AttributeMapping(self):
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': secret_wRDP,
            })
        
    # TODO mock run_check
    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *args: get_mock_dfu('secret'))
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *args, **kwargs: get_mock_kbr())
    def test_AmpliconSet_no_taxonomy_no_row_AttributeMapping(self):
        with self.assertRaises(NoTaxonomyException) as cm:
            ret = self.serviceImpl.run_FAPROTAX(
                self.ctx, {
                    **self.params_ws,
                    'input_upa': secret,
                })
            

    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *args: get_mock_dfu('17770'))
    @patch_('kb_faprotax.util.workflow.run_check', new=get_mock_run_check('17770'))
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *args, **kwargs: get_mock_kbr())
    def test_AmpliconSet_input_against_reference(self):
        '''
        Check results against answers if full pipeline is run (i.e., no mocks)
        
        Some things that have been hardcoded:
        * `source` is `FAPROTAX Functions`
        * FAPROTAX functions delimiter is `','`
        '''
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': _17770,
                }
            )

        from kb_faprotax.util.workflow import run_check

        # only check results against reference
        # if full pipeline has been run
        if (isinstance(Var.dfu, unittest.mock.NonCallableMagicMock) or
            isinstance(Var.kbr, unittest.mock.NonCallableMagicMock) or
            hasattr(run_check, 'side_effect')):
            return

        logging.info('Comparing traits in AttributeMapping to answers')

        # hardcoded
        # should match in workflow.py
        source = 'FAPROTAX Functions'
        dlm = ','

        # answers to compare against
        def parse_answers_file(answers_flpth='/kb/module/test/data/reference/OTUMetaData_reduced.tsv'):
            answers_df = pd.read_csv(
                answers_flpth, sep='\t', header=0, index_col='#OTU ID', usecols=['#OTU ID', 'FAPROTAX Traits']).fillna('')
            answers_d = answers_df.to_dict(orient='index')
            answers_d = {key: value['FAPROTAX Traits'] for key, value in answers_d.items()}
            return answers_d

        # load KBaseReport
        attrmap_upa_new = Var.dfu.get_objects({
            'object_refs': [ret[0]['report_ref']]
        })['data'][0]['data']['objects_created'][0]['ref']

        # load AttributeMapping
        row_attrmap = AttributeMapping(attrmap_upa_new)
        instances_d = row_attrmap.obj['instances']
        attribute_d_l = row_attrmap.obj['attributes']

        # find index in attribute list
        for i, attribute_d in enumerate(attribute_d_l):
            if attribute_d['attribute'] == source:
                ind = i

        # id to attribute
        res_id2attr_d = {id: attr_l[ind] for id, attr_l in instances_d.items()}
        ans_id2attr_d = self.parse_answers_file()

        # check ids equal set
        res_ids = sorted(list(res_id2attr_d.keys()))
        ans_ids = sorted(list(ans_id2attr_d.keys()))
        self.assertTrue(res_ids == ans_ids)

        #
        html_l = []
        for id in res_ids:
            res = sorted(res_id2attr_d[id].split(dlm))
            ans = sorted(ans_id2attr_d[id].split(':'))

            # no dups
            self.assertTrue(len(set(res)) == len(res))

            if res != ans:
                res = set(res)
                ans = set(ans)
                
                ##
                ## check all extras are left-off parent functions
                self.assertTrue(ans.issubset(res))
                extras = list(res - ans)
                # check that extra is a substring of something in `ans`
                # TODO ascertain that left-off is parent of something in `ans`
                for extra in extras:
                    self.assertTrue(len([None for function in ans if extra in function]) > 0)

                # make a diff of extras in `res`
                html = '<p>' + ','.join([res_ if res_ in ans else '<b>' + res_ + '</b>' for res_ in res]) + '</p>'
                html_l.append(html)

        # dedup for visual ease
        html_l = list(set(html_l))

        with open(f'/kb/module/work/tmp/diff.html', 'w') as fp:
            fp.write('\n'.join(html_l))

 
    @patch.dict('kb_faprotax.util.kbase_obj.Var', values={'dfu': get_mock_dfu('17770')})
    def test_AmpliconSet_methods(self):
        '''
        Superficial testing
        '''

        ##
        ## good input

        amp_set = AmpliconSet(_17770)

        id_l = ['fffa52555f0d542613a26955a558d76d', 'fffd0e743b783e58fb70f6a5c5e41901']
        taxStr_l = ['\
D_0__Bacteria; D_1__Proteobacteria; D_2__Gammaproteobacteria; D_3__Steroidobacterales; D_4__Steroidobacteraceae; \
D_5__Steroidobacter; D_6__Steroidobacter:u; D_7__Steroidobacter:u; D_8__Steroidobacter:u; D_9__Steroidobacter:u; \
D_10__Steroidobacter:u', 'D_0__Bacteria; D_1__Proteobacteria; D_2__Deltaproteobacteria; D_3__Sva0485; D_4__Sva0485:u; \
D_5__Sva0485:u; D_6__Sva0485:u; D_7__Sva0485:u; D_8__Sva0485:u; D_9__Sva0485:u; D_10__Sva0485:u']

        self.assertTrue(amp_set.get_taxStr_l([]) == [])
        self.assertTrue(amp_set.get_taxStr_l(id_l) == taxStr_l) 

        tax2ids_d = amp_set._get_tax2ids_d()
        for id, taxStr in zip(id_l, taxStr_l):
            self.assertTrue(id in tax2ids_d[taxStr])

        ##
        ## no taxonomy field input

        mock_dfu.get_objects.side_effect = lambda params: { # this mock is a global
                'data': [
                    {
                        'data': {
                            'amplicon_matrix_ref': '-1/-2/-3',
                            'amplicons': {
                                'amplicon0': {
                                    'consensus_sequence': 'aaggcctt',
                                    'taxonomy': {
                                        'lineage': [], # TODO disallow?
                                    },
                                },
                                'amplicon1': {
                                    'consensus_sequence': 'aaggcctt',
                                    'taxonomy': {
                                        'lineage': ['this', 'is bogus', 'tax'], # TODO disallow?
                                    },
                                },
                                'amplicon2': {
                                    'consensus_sequence': 'aaggcctt',
                                    'taxonomy': {},
                                },
                            },
                        },
                        'info': [
                            -666,
                            'AmpliconSet_NoTaxonomy_Dummy_Name',
                            'Dummy.AmpliconSet-1.0',
                        ]
                    }
                ]
        }

        amp_set = AmpliconSet('du/mm/y')

        with self.assertRaises(NoTaxonomyException) as cm:
            amp_set.get_taxStr_l(['amplicon0', 'amplicon1', 'amplicon2'])
            self.assertEquals('`%s`' % 'amplicon2' in str(cm.exception)) # TODO should amplicon0 and amplion1 trigger exceptions?


    @patch.dict('kb_faprotax.util.kbase_obj.Var', values={'dfu': get_mock_dfu('17770')})
    def test_AmpliconMatrix_methods(self):
        amp_set = AmpliconSet(_17770)
        amp_mat = AmpliconMatrix(_17770_AmpMat, amp_set)

        # superficially test `to_OTU_table`
        df1 = pd.read_csv(
            os.path.join(testData_dir, 'by_dataset_input/17770/return/otu_table.tsv'), sep='\t', index_col='taxonomy') 
        df2 = amp_mat.to_OTU_table()

        self.assertTrue(df1.columns.tolist() == df2.columns.tolist())
        self.assertTrue(df1.index.tolist() == df2.index.tolist()) # index is `taxonomy`
        self.assertTrue(df1['OTU_Id'].tolist() == df2['OTU_Id'].tolist()) # first column is `OTU_Id`

        data1 = df1.iloc[:,1:].values # rest of columns are data
        data2 = df2.iloc[:,1:].values

        self.assertTrue(np.all(np.abs(data1 - data2) < 10**-15)) # python floating point issue 
                                                                 # e.g., 0.1 + 0.1 + 0.1 => 0.30000000000000004


    @patch.dict('kb_faprotax.util.kbase_obj.Var', values={'dfu': get_mock_dfu(None), 'warnings': []})
    def test_AttributeMapping_methods(self):
        # set `get_objects` to return something simpler and independent
        obj = {
            "data": [
                {
                    "data": {
                        "attributes": [
                            {
                                "attribute": "celestial body",
                                "source": "upload",
                            },
                            {
                                "attribute": "ecosystem",
                                "source": "upload",
                            },
                        ],
                        "instances": {
                            "fffffffffffffffffffffffffffffffd": [
                                "Ganymede",
                                "tundra",
                            ],
                            "fffffffffffffffffffffffffffffffe": [
                                "Mars",
                                "rainforest",
                            ],
                            "ffffffffffffffffffffffffffffffff": [
                                "Saturn",
                                "chaparral",
                            ],
                        }
                    },
                    "info": [
                        -1,
                        "AttributeMapping_dummy_name",
                    ]
                }
            ]
        }
        mock_dfu.get_objects.side_effect = lambda params: obj

        id2attr_d = {'f' * 31 + 'd': 'black hole', 'f' * 31 + 'e': 'quasar', 'f' * 32: 'dark matter'}

        attr_map = AttributeMapping('-111/-1/-1')

        ind = attr_map.add_attribute_slot('celestial body', 'unit testing')
        attr_map.update_attribute(ind, id2attr_d)
        self.assertTrue(len(Var.warnings) == 0)
        
        ind = attr_map.add_attribute_slot('celestial body', 'upload')
        attr_map.update_attribute(ind, id2attr_d)
        self.assertTrue(len(Var.warnings) == 1)


                 



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!! SETUP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_faprotax'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error:w

        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_faprotax',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.wsName = 'kb_faprotax_' + str(uuid.uuid4())                                                 
        cls.wsId = cls.wsClient.create_workspace({'workspace': cls.wsName})[0]                      
        cls.params_ws = {                                                                           
            'workspace_id': cls.wsId,                                                               
            'workspace_name': cls.wsName,                                                           
            }                                                                                       
        cls.serviceImpl = kb_faprotax(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_kb_faprotax_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa


    @classmethod
    def list_tests(cls):
        return [key for key, value in cls.__dict__.items() if type(key) == str and key.startswith('test') and callable(value)]

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        dec = '!!!' * 220
        print(dec, "DON'T FORGET TO SEE DIFF, HTML REPORT(S)", dec)
        print('Tests run:', cls.list_tests())


    def shortDescription(self):
        '''Override unittest using test*() docstrings in lieu of test*() method name in output summary'''
        return None


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!! select what to run !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''
When you just want to run certain tests,
e.g., filter to tests in `run_tests`

Comment out parts like `delattr` to deactivate
'''
unit_tests = ['test_parse_faprotax_functions', 'test_run_check',]
AmpliconSet_tests = unit_tests + ['test_AmpliconSet_input_against_reference', 
    'test_AmpliconSet_no_row_AttributeMapping', 'test_AmpliconSet_no_taxonomy_no_row_AttributeMapping',
    'test_AmpliconSet_methods', 'test_AmpliconMatrix_methods', 'test_AttributeMapping_methods']
GenomeSet_tests = unit_tests + ['test_GenomeSet_input', 'test_GenomeSet_methods', 'test_Genome_methods',
    'test_dummy_abundance', 'test_dup_GenomeSet']
run_tests = ['test_AmpliconMatrix_methods']

for key, value in kb_faprotaxTest.__dict__.copy().items():
    if key.startswith('test') and callable(value):
        if key not in AmpliconSet_tests:
            delattr(kb_faprotaxTest, key)
            pass





