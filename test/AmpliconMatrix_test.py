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
from pytest import raises

from kb_faprotax.kb_faprotaxImpl import kb_faprotax
from kb_faprotax.kb_faprotaxServer import MethodContext
from kb_faprotax.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace

from kb_faprotax.util.workflow import run_check
from kb_faprotax.util.error import * # exception library
from kb_faprotax.util.dprint import dprint, where_am_i
from kb_faprotax.util.varstash import Var # dict-like dot-access app globals
from kb_faprotax.util.kbase_obj import GenomeSet, Genome, AmpliconMatrix, AttributeMapping
from mock import * # mock business
from upa import * # upa library
import config



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



class TestCase(config.BaseTest):


# TODO test: missing, rel abund vs raw abund, no tax, attributes at end (spot)

    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma50by30_RDPClsf'))
    @patch_('kb_faprotax.util.workflow.run_check', new=get_mock_run_check('enigma50by30_RDPClsf'))
    @patch_('kb_faprotax.kb_faprotaxImpl.GenericsAPI', new=lambda *a, **k: get_mock_gapi())
    @patch_('kb_faprotax.kb_faprotaxImpl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu())
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test_tax_field(self):
        ret = config.get_serviceImpl().run_FAPROTAX(
            config.ctx, 
            {
                **self.ws,
                'input_upa': enigma50by30_RDPClsf,
                'tax_field': 'RDP Classifier taxonomy, conf=0.777, gene=silva_138_ssu, minWords=default',
                'output_amplicon_matrix_name': 'a_name',
            }
        )

        assert len(Var.params_report.objects_created) == 4 


    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma50by30_noAttrMaps_noSampleSet'))
    def test_noRowAttributeMapping(self):
        with raises(NoWsReferenceException):
            ret = config.get_serviceImpl().run_FAPROTAX(
                config.ctx, 
                {
                    **self.ws,
                    'input_upa': enigma50by30_noAttrMaps_noSampleSet,
                }
            )
            
    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma50by30'))
    @patch_('kb_faprotax.util.workflow.run_check', new=get_mock_run_check('enigma50by30'))
    @patch_('kb_faprotax.kb_faprotaxImpl.GenericsAPI', new=lambda *a, **k: get_mock_gapi())
    @patch_('kb_faprotax.kb_faprotaxImpl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu())
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test(self):
        ret = config.get_serviceImpl().run_FAPROTAX(
            config.ctx, 
            {
                **self.ws,
                'input_upa': enigma50by30,
                'tax_field': 'taxonomy',
                'output_amplicon_matrix_name': 'a_name',
            }
        )
        
        assert len(Var.params_report.objects_created) == 4
        

    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma17770by511'))
    @patch_('kb_faprotax.util.workflow.run_check', new=get_mock_run_check('enigma17770by511')) #
    @patch_('kb_faprotax.kb_faprotaxImpl.GenericsAPI', new=lambda *a, **k: get_mock_gapi())
    @patch_('kb_faprotax.kb_faprotaxImpl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu())
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test_against_reference(self): # TODO
        ret = config.get_serviceImpl().run_FAPROTAX(
            config.ctx, 
            {
                **self.ws,
                'input_upa': enigma17770by511,
                'tax_field': 'taxonomy',
                'output_amplicon_matrix_name': 'a_name',
            }
        )

        assert len(Var.params_report.objects_created) == 4 



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!! SETUP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    @classmethod
    def list_tests(cls):
        return [key for key, value in cls.__dict__.items() if type(key) == str and key.startswith('test') and callable(value)]

    @classmethod
    def tearDownClass(cls):
        super(cls, cls).tearDownClass()

        dec = '!!!' * 220
        print(dec, "DON'T FORGET TO SEE DIFF, HTML REPORT(S)", dec)
        skipped_tests = list(set(all_tests) - set(cls.list_tests()))
        print('* All tests (%d): %s' % (len(all_tests), all_tests))
        print('* Tests skipped (%d): %s' % (len(skipped_tests), skipped_tests))
        print('* Tests run (%d): %s' % (len(cls.list_tests()), cls.list_tests()))
        




#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!! select what to run !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
all_tests = []
for key, value in TestCase.__dict__.items():
    if key.startswith('test') and callable(value):
        all_tests.append(key)

AmpliconMatrix_integration_tests = [
        'test',
    'test_tax_field', 
    'test_noRowAttributeMapping',
    'test_against_reference'
]

run_tests = ['test']

for test in all_tests:
    if test not in run_tests:#(unit_tests + AmpliconMatrix_integration_tests):
            delattr(TestCase, test)
            pass





