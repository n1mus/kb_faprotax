# -*- coding: utf-8 -*-
import os
import logging
import time
import unittest
import uuid
import pandas as pd

from configparser import ConfigParser
from kb_faprotax.kb_faprotaxImpl import kb_faprotax
from kb_faprotax.kb_faprotaxServer import MethodContext
from kb_faprotax.authclient import KBaseAuth as _KBaseAuth
from kb_faprotax.util.dprint import dprint
from kb_faprotax.util.varstash import Var
from kb_faprotax.util.kbase_obj import AttributeMapping

from installed_clients.WorkspaceClient import Workspace


test_amplicon_set_upa = "39332/58/1"
test_amplicon_matrix_upa = "39332/57/2"


enigma_amp_set_upa = "48363/2/1"


params_debug = {
    'skip_run': True,
    'skip_retFiles': True,
    }


class kb_faprotaxTest(unittest.TestCase):

    def test(self):
        ret = self.serviceImpl.faprotax(
            self.ctx, {
                'workspace_name': self.wsName,
                'amplicon_set_upa': enigma_amp_set_upa,
                **self.params_ws,
                #**params_debug,
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


        with open(f'/kb/module/work/tmp/{uuid.uuid4()}.html', 'w') as fp:
            fp.write('\n'.join(html_l))

                  


    @staticmethod
    def parse_answers_file():
        answers_flpth = '/kb/module/test/data/OTUMetaData_reduced.tsv'
        answers_df = pd.read_csv(
            answers_flpth, sep='\t', header=0, index_col='#OTU ID', usecols=['#OTU ID', 'FAPROTAX Traits']).fillna('')
        answers_d = answers_df.to_dict(orient='index')
        answers_d = {key: value['FAPROTAX Traits'] for key, value in answers_d.items()}
        return answers_d



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
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        print('!!!!!' * 500 + ' DON\'T FORGET TO INSPECT DIFF ' + '!!!!!' * 500)
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
