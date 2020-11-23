# -*- coding: utf-8 -*-
from configparser import ConfigParser
import os
import logging
import time
import unittest
import uuid
import pandas as pd
import numpy as np
from shutil import rmtree, copytree
import json
from unittest.mock import patch
import pandas as pd
from urllib.parse import urlparse
import re

from kb_faprotax.kb_faprotaxImpl import kb_faprotax
from kb_faprotax.kb_faprotaxServer import MethodContext
from kb_faprotax.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace

from kb_faprotax.util.workflow import run_check
from kb_faprotax.util.error import * # exception library
from kb_faprotax.util.dprint import dprint, where_am_i
from kb_faprotax.util.varstash import Var # dict-like dot-access app globals
from kb_faprotax.util.kbase_obj import GenomeSet, Genome, AmpliconMatrix, AttributeMapping
from util.mock import * # mock business
from util.upa import * # upa library



######################################
######################################
######### TOGGLE PATCH ###############
######################################
###################################### 
do_patch = False

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
########################## unit tests ##############################################################
####################################################################################################

    ####################
    ####################
    # doesn't need mocking - fast
    def test_parse_faprotax_functions(self):
        '''Test function used to parse FAPROTAX's predicted functions from its output'''
        from kb_faprotax.util.workflow import parse_faprotax_functions

        flnm = 'groups2records_dense.tsv'
        flpth = os.path.join(testData_dir, 'by_dataset_input/refseq/return/FAPROTAX_output', flnm)
        
        r2g_d = parse_faprotax_functions(flpth)

        taxonomy = '\
cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Moraxellaceae; \
Acinetobacter; Acinetobacter calcoaceticus/baumannii complex; Acinetobacter pittii'
        functions = '\
aerobic_chemoheterotrophy,human_pathogens_septicemia,human_pathogens_pneumonia,human_pathogens_nosocomia,\
human_pathogens_all,animal_parasites_or_symbionts,aromatic_compound_degradation,chemoheterotrophy'

        self.assertTrue(r2g_d[taxonomy] == functions)


    ####################
    ####################
    # doesn't need mocking - fast
    def test_run_check(self):
        '''Test function that runs FAPROTAX'''

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


 
    ####################
    ####################
    # TODO superficially test ctor
    def test_GenomeSet_methods(self):
        pass
       

    ####################
    ####################
    def test_Genome_methods(self):
        '''
        Genome methods are quite simple, 
        they are tested in other places
        '''
        pass


    ####################
    ####################
    def test_dummy_OTU_table_abundance(self):
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
 
    # TODO bug with patch.dict interfering with each other

    ####################
    ####################
    def test_AmpliconMatrix_methods(self):
        '''
        dprint('Var.tmp # globals', run=globals())
        from kb_faprotax.util.kbase_obj import Var; dprint('Var.tmp # locals', run=locals())
        '''

        # use as cm since patch.dict decs interfering with e/o for some reason
        with patch.dict('kb_faprotax.util.kbase_obj.Var', values={
                        'tmp': 'test_AmpMat_methods', 'dfu': get_mock_dfu('enigma50by30')}):
            amp_mat = AmpliconMatrix(enigma50by30)
            row_attr_map = AttributeMapping(enigma50by30_rowAttrMap, amp_mat)

        # id_l and tax_l correspond to AmpliconMatrix rows
        # which don't necessarily map 1-to-1 [bijectively] with row AttributeMapping
        ind = row_attr_map.get_attr_ind('taxonomy')
        id_l = amp_mat.obj['data']['row_ids']
        tax_l = row_attr_map.get_ordered_tax_l(ind, id_l)

        # superficially test `to_OTU_table`
        df1 = pd.read_csv(
            os.path.join(testData_dir, 'by_dataset_input/enigma50by30/return/otu_table.tsv'), 
            sep='\t', index_col='taxonomy') 
        df2 = amp_mat.to_OTU_table(tax_l)

        self.assertTrue(df1.columns.tolist() == df2.columns.tolist())
        self.assertTrue(df1.index.tolist() == df2.index.tolist()) # index is `taxonomy`
        self.assertTrue(df1['OTU_Id'].tolist() == df2['OTU_Id'].tolist()) # first column is `OTU_Id`

        data1 = df1.iloc[:,1:].values # rest of columns are data
        data2 = df2.iloc[:,1:].values

        self.assertTrue(np.allclose(data1, data2) and np.allclose(data2, data1))

    ####################
    ####################
    def test_AttributeMapping_methods(self):

        # update with
        id2attr = {'amplicon_id_2': 'black hole', 'amplicon_id_5': 'quasar', 'amplicon_id_9': 'dark matter'}

        # use as cm since patch.dict decs interfering with e/o for some reason
        with patch.dict('kb_faprotax.util.kbase_obj.Var', values={
                        'tmp': 'AttrMap_methods', 'dfu': get_mock_dfu('dummy10by8'), 'warnings': [], }):
            amp_mat = AmpliconMatrix(dummy10by8)
            row_attr_map = AttributeMapping(dummy10by8_rowAttrMap, amp_mat)

        # update in new attribute/upload slot
        with self.subTest():
            ind0 = row_attr_map.add_attribute_slot_warn('celestial body', 'unit testing')
            row_attr_map.map_update_attribute(ind0, id2attr)

            self.assertTrue(ind0 == row_attr_map.attributes_length-1, '%d vs %d' % (ind0, row_attr_map.attributes_length-1))
            self.assertTrue(len(Var.warnings) == 0)
            row_attr_map._check_attr_consistent(ind0, 'celestial body', 'unit testing', id2attr)

        # update in existing attribute/upload slot
        with self.subTest():
            ind1 = row_attr_map.add_attribute_slot_warn('celestial body', 'upload')
            row_attr_map.map_update_attribute(ind1, id2attr)

            self.assertTrue(ind1 < row_attr_map.attributes_length-1)
            self.assertTrue(len(Var.warnings) == 1)
            row_attr_map._check_attr_consistent(ind1, 'celestial body', 'upload', id2attr)

        #
        self.assertTrue(row_attr_map.get_attr_ind('gene') is None)

        # gets first attribute match
        self.assertTrue(row_attr_map.get_attr_ind('celestial body') == ind1)

    ####################
    ####################
    def test_AmpliconMatrix_validation(self):
        '''
        Test validation of amplicon table in AmpliconMatrix
        Should be (1) count data,  missing (None) allowed
        Can assume that obj holds data in list of lists of numeric/None
        '''
        
        logging.info('Testing with test_AmpliconMatrix_validation')

        # use as cm since patch.dict decs interfering with e/o for some reason
        with patch.dict('kb_faprotax.util.kbase_obj.Var', values={
                        'tmp': 'AmpMat_validation', 'dfu': get_mock_dfu('dummy10by8')}):
            amp_mat = AmpliconMatrix(dummy10by8)

        
        with self.assertRaisesRegex(ValidationException, '[Ii]nteger'): amp_mat.validate_amplicon_abundance_data()


        amp_mat.obj['data']['values'] = [0.0, 0.0, 1319.0, 1.0] # float
        amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [0, 0, 1319, 1] # int
        amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, 0., 0., 1319., 1.] # float, with missing
        amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, 0, 0, 1319, 1] # int, with missing
        amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, 0, -0., 1319.0, 1] # int/float, with missing
        amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, None, 0, -0, 0.0, 0, 0.] # 0s, with missing
        amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, 0.999999999] # close enough
        amp_mat.validate_amplicon_abundance_data() 

        amp_mat.obj['data']['values'] = [None, -0.0000000001] # close enough
        amp_mat.validate_amplicon_abundance_data() 

        
        amp_mat.obj['data']['values'] = [0.9]
        with self.assertRaisesRegex(ValidationException, '[Ii]nteger'): amp_mat.validate_amplicon_abundance_data() 

        amp_mat.obj['data']['values'] = [-1]
        with self.assertRaisesRegex(ValidationException, '[Nn]egative'): amp_mat.validate_amplicon_abundance_data() 

        amp_mat.obj['data']['values'] = [None, None, None]
        with self.assertRaisesRegex(ValidationException, 'all missing'): amp_mat.validate_amplicon_abundance_data() 

        amp_mat.obj['data']['values'] = [None]
        with self.assertRaisesRegex(ValidationException, 'all missing'): amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, -1]
        with self.assertRaisesRegex(ValidationException, '[Nn]egative'): amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, None, 1.00001]
        with self.assertRaisesRegex(ValidationException, '[Ii]nteger'): amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [-1.0, 0, 1319]
        with self.assertRaisesRegex(ValidationException, '[Nn]egative'): amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, 0, 1, 2, 3, 4.5]
        with self.assertRaisesRegex(ValidationException, '[Ii]nteger'): amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, 0.0, 1.0, 2.0, 3.0, 4.00001] # 4.00001 would pass with np.allclose default rtol
        with self.assertRaisesRegex(ValidationException, '[Ii]nteger'): amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [0.9, -0.9]
        with self.assertRaises(ValidationException): amp_mat.validate_amplicon_abundance_data() 




####################################################################################################
########################### GenomeSet input workflow ###############################################
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
    def test_dup_GenomeSet(self):
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': refseq_dup,
            })

        self.assertTrue(msg_dupGenomes in Var.warnings)



####################################################################################################
########################## AmpliconMatrix input #######################################################
####################################################################################################

# TODO test: missing, rel abund vs raw abund, no tax

    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma50by30_RDPClsf'))
    @patch_('kb_faprotax.kb_faprotaxImpl.GenericsAPI', new=lambda *a, **k: get_mock_gapi())
    @patch_('kb_faprotax.kb_faprotaxImpl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu())
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test_AmpliconMatrix_tax_field(self):
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': enigma50by30_RDPClsf,
                'tax_field': 'RDP Classifier taxonomy, conf=0.777, gene=silva_138_ssu, minWords=default', # TODO don't expose minWords in UI/attribute ?
                'output_amplicon_matrix_name': 'a_name',
            })

        self.assertTrue(len(Var.params_report.objects_created) == 4)


    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma50by30_noAttrMaps_noSampleSet'))
    def test_AmpliconMatrix_noRowAttributeMapping(self):
        with self.assertRaises(NoWsReferenceException) as cm:
            ret = self.serviceImpl.run_FAPROTAX(
                self.ctx, {
                    **self.params_ws,
                    'input_upa': enigma50by30_noAttrMaps_noSampleSet,
                })
            
    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma50by30_noSampleSet'))
    @patch_('kb_faprotax.kb_faprotaxImpl.GenericsAPI', new=lambda *a, **k: get_mock_gapi())
    @patch_('kb_faprotax.kb_faprotaxImpl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu())
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test_AmpliconMatrix_noSampleSet(self):
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': enigma50by30_noSampleSet,
                'tax_field': 'taxonomy',
            })

        self.assertTrue(len(Var.params_report.objects_created) == 3)

    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma17770by511'))
    @patch_('kb_faprotax.util.workflow.run_check', new=get_mock_run_check('enigma17770by511')) #
    @patch_('kb_faprotax.kb_faprotaxImpl.GenericsAPI', new=lambda *a, **k: get_mock_gapi())
    @patch_('kb_faprotax.kb_faprotaxImpl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu())
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test_AmpliconMatrix_against_reference(self):
        '''
        Check results against answers if full pipeline is run (i.e., no mocks)
        
        Some things that have been hardcoded:
        * `source` is `FAPROTAX Functions`
        * FAPROTAX functions delimiter is `','`
        '''
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': enigma17770by511,
                'tax_field': 'taxonomy',
                }
            )

        self.assertTrue(len(Var.params_report.objects_created) == 4)

        ##
        ## Only important to check results if actual run occcurred (as opposed to mocked)

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
        dummy_amp_mat = MagicMock(upa=enigma17770by511)
        row_attrmap = AttributeMapping(attrmap_upa_new, dummy_amp_mat)
        instances_d = row_attrmap.obj['instances']
        attribute_d_l = row_attrmap.obj['attributes']

        # find index in attribute list
        for i, attribute_d in enumerate(attribute_d_l):
            if attribute_d['attribute'] == source:
                ind = i

        # id to attribute
        res_id2attr_d = {id: attr_l[ind] for id, attr_l in instances_d.items()}
        ans_id2attr_d = parse_answers_file()

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
        skipped_tests = list(set(all_tests) - set(cls.list_tests()))
        print('* All tests (%d): %s' % (len(all_tests), all_tests))
        print('* Tests skipped (%d): %s' % (len(skipped_tests), skipped_tests))
        print('* Tests run (%d): %s' % (len(cls.list_tests()), cls.list_tests()))
        

    def shortDescription(self):
        '''Override unittest using test*() docstrings in lieu of test*() method name in output summary'''
        return None


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!! select what to run !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
all_tests = []
for key, value in kb_faprotaxTest.__dict__.items():
    if key.startswith('test') and callable(value):
        all_tests.append(key)
'''
When you just want to run certain tests,
e.g., filter to tests in `run_tests`

Comment out parts like `delattr` to deactivate
'''
unit_tests = [
    'test_parse_faprotax_functions', 'test_run_check',
    'test_AmpliconMatrix_validation', 'test_AmpliconMatrix_methods', 'test_AttributeMapping_methods',
    'test_Genome_methods', 'test_GenomeSet_methods',
]
AmpliconMatrix_integration_tests = [
    'test_AmpliconMatrix_tax_field', 
    'test_AmpliconMatrix_noRowAttributeMapping', 'test_AmpliconMatrix_noSampleSet',
    'test_AmpliconMatrix_against_reference'
]
GenomeSet_integration_tests = [
    'test_GenomeSet_input',
    'test_dup_GenomeSet',
    'test_dummy_OTU_table_abundance', # TODO mock this so it's a unit test
]
run_tests = ['test_AmpliconMatrix_against_reference']

for test in all_tests:
    if test not in run_tests:#(unit_tests + AmpliconMatrix_integration_tests):
            delattr(kb_faprotaxTest, test)
            pass





