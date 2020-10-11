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
import sys

from kb_faprotax.kb_faprotaxImpl import kb_faprotax
from kb_faprotax.kb_faprotaxServer import MethodContext
from kb_faprotax.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace

from kb_faprotax.util.debug import dprint, where_am_i
from kb_faprotax.util.config import Var # dict-like dot-access app globals
from kb_faprotax.impl.error import * # exception library
from kb_faprotax.impl.kbase_obj import GenomeSet, Genome, AmpliconSet, AmpliconMatrix, AttributeMapping
from .mock import * # mock business
from .upa import * # upa library

dec = '###' * 5                                                                                     
skipped_tests = []                                                                                  
                                                                                                    
def tag_kb_env(e):                                                                                  
    @functools.wraps(f) # preserve wrapped's function name                                          
    def decorator(f):                                                                               
        f.kb_env = e                                                                                
        return f                                                                                    
    return decorator                                                                                
                                                                                                    
def skip_cond(select_run=None, regex=None, exclude_regex=None, kb_env=None):                        
    def decorator(f):                                                                               
        @functools.wraps(f) # preserve wrapped's function name                                      
        def f_new(self, *a, **kw):                                                                  
            if kb_env is not None and not hasattr(f, "kb_env"):                                     
                raise Exception("Tag function (e.g., @tag_kb_env('ci')) with kb_env first to skip with this feature")
                                                                                                    
            if kb_env is not None and f.kb_env is not None and kb_env != f.kb_env:                  
                skipped_tests.append(f.__name__)                                                    
                self.skipTest("Test does not operate in this KBase environment")                    
            if select_run is not None and f.__name__ not in select_run:                             
                skipped_tests.append(f.__name__)                                                    
                print(dec, 'Skipping test %s because not in select_run' % f.__name__, dec)          
                self.skipTest("Test is not in list of select_run")                                  
            if regex is not None and re.search(regex, f.__name__) is None:                          
                skipped_tests.append(f.__name__)                                                    
                print(dec, 'Skipping test %s because regex has no hit' % f.__name__, dec)           
                self.skipTest("Test name does not have hit from regex %s" % regex)                  
                                                                                                    
            print(dec, 'Running test %s' % f.__name__, dec)                                         
            f(self, *a, **kw)                                                                       
        return f_new                                                                                
    return decorator                                                                                
                                                                                                    
kb_env = None
select_run = None #['test_transform_pipeline']                                                      
regex = None #'transform'



# TODO allow toggling patching as much as possible, i.e., for unit testing 
# TODO dummy test datasets (to file) = less file io (slow?) ... leave big datasets for ?

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

class KBFaprotaxTest(unittest.TestCase):

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
    def list(cls):
        '''
        List tests
        '''
        return [key for key, value in cls.__dict__.items() if type(key) == str and key.startswith('test') and callable(value)]

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        dec = '!!!' * 220
        print(dec, "DON'T FORGET TO SEE DIFF, HTML REPORT(S)", dec)
        print('All tests in cls %s: %s' % (cls.__name__, cls.list()))


    def shortDescription(self):
        '''Override unittest using test*() docstrings in lieu of test*() method name in output summary'''
        return None

