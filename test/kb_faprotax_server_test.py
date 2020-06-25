# -*- coding: utf-8 -*-
from configparser import ConfigParser
import os
import logging
import time
import unittest
from unittest.mock import patch, create_autospec
import uuid
import pandas as pd

from kb_faprotax.kb_faprotaxImpl import kb_faprotax
from kb_faprotax.kb_faprotaxServer import MethodContext
from kb_faprotax.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport

from kb_faprotax.util.error import *
from kb_faprotax.util.message import *
from kb_faprotax.util.dprint import dprint
from kb_faprotax.util.varstash import Var
from kb_faprotax.util.kbase_obj import AttributeMapping




##################################
##################################
_17770 = '48666/2/9' # AmpliconSet containing 17770 entries
first50 = "48402/9/2" # AmpliconSet containing first 50 of 17770 entries. row AttributeMapping has all 1770 entries (?)
secret = '49926/6/1' # AmpliconSet from collaborator. No taxonomy or row AttributeMapping. Do not share
secret_wRDP = '49926/9/3' # AmpliconSet from collaborator, with taxonomy and row AttributeMapping. Do not share

refseq = '43623/61/2' # GenomeSet with RefSeq prokaryote reference genomes
##################################
##################################





mock_dfu = create_autospec(DataFileUtil, spec_set=True, instance=True)
mock_kbr = create_autospec(KBaseReport, spec_set=True, instance=True)
#mock_run_check


def set_17770(mock_dfu):
    mock_dfu.reset_mock(return_value=True, side_effect=True)
    return mock_dfu




class kb_faprotaxTest(unittest.TestCase):

####################################################################################################
########################### GenomeSet input ########################################################
####################################################################################################

    #@patch('kb_faprotax.util.kbase_obj.Var.dfu', new=set_17770(mock_dfu))
    def test_GenomeSet_input(self):
        ret = self.serviceImpl.faprotax(
            self.ctx, {
                **self.params_ws,
                'input_upa': refseq,
                }
            )
        return
   

####################################################################################################
########################## AmpliconSet input #######################################################
####################################################################################################

    def test_attribute_and_source_exists(self):
        pass
   

    def test_add_new_attribute(self): # TODO tie to source
        ret = self.serviceImpl.faprotax(
            self.ctx, {
                **self.params_ws,
                'input_upa': secret_wRDP,
            })
        

    def test_no_taxonomy_no_AttributeMapping(self):
        with self.assertRaises(NoTaxonomyException) as cm:
            ret = self.serviceImpl.faprotax(
                self.ctx, {
                    **self.params_ws,
                    'input_upa': secret,
                })
            

    def test_against_reference(self):
        ret = self.serviceImpl.faprotax(
            self.ctx, {
                **self.params_ws,
                'input_upa': _17770,
                'return_test_info': True,
                }
            )

        logging.info('Comparing traits in AttributeMapping to answers')

        # load AttributeMapping
        row_attrmap = AttributeMapping(Var.objects_created[0]['ref'])
        instances_d = row_attrmap.obj['instances']
        attribute_d_l = row_attrmap.obj['attributes']

        # find index in attribute list
        for i, attribute_d in enumerate(attribute_d_l):
            if attribute_d['attribute'] == 'FAPROTAX Traits':
                ind = i

        # id to attribute
        results_d = {id: attr_l[ind] for id, attr_l in instances_d.items()}

        # id to traits
        answers_d = self.parse_answers_file()

        html_l = []

        for id in answers_d:
            assert id in results_d

            res = results_d[id]
            ans = answers_d[id]

            if res != ans:
                res_l = res.split(':')
                ans_l = ans.split(':')
                assert set(ans_l).issubset(res_l)
                
                html = '<p>' + ':'.join([res if res in ans_l else '<b>' + res + '</b>' for res in res_l]) + '</p>'
                html_l.append(html)

        html_l = list(set(html_l))


        with open(f'/kb/module/work/tmp/diff.html', 'w') as fp: # TODO automate the left-off parent function detection so you don't have to look at this every time
            fp.write('\n'.join(html_l))

                  


####################################################################################################

    @staticmethod
    def parse_answers_file():
        answers_flpth = '/kb/module/test/data/OTUMetaData_reduced.tsv'
        answers_df = pd.read_csv(
            answers_flpth, sep='\t', header=0, index_col='#OTU ID', usecols=['#OTU ID', 'FAPROTAX Traits']).fillna('')
        answers_d = answers_df.to_dict(orient='index')
        answers_d = {key: value['FAPROTAX Traits'] for key, value in answers_d.items()}
        return answers_d



####################################################################################################


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
        # it'll result in a NoneType error
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
        dprint('cls.wsId', run=locals())  
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
        print(dec, "DON'T FORGET TO SEE DIFF, HTML", dec)
        print('All tests:', cls.list_tests())




############################ select what to run ####################################################
'''
When you just want to run certain tests,
e.g., filter to tests in `run_tests`

Comment out parts like `delattr` to deactivate
'''
AmpliconSet_tests = []
GenomeSet_tests = ['test_GenomeSet_input']
run_tests = ['test_against_reference']

for key, value in kb_faprotaxTest.__dict__.copy().items():
    if key.startswith('test') and callable(value):
        if key not in run_tests:
            delattr(kb_faprotaxTest, key)
            pass





