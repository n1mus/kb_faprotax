from unittest.mock import patch, create_autospec, Mock
import os
from shutil import rmtree, copytree
import logging
import json

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace

from kb_faprotax.util.dprint import dprint
from kb_faprotax.util.varstash import Var
from kb_faprotax.util.workflow import run_check
from .upa import *


##################################
##################################
testData_dir = '/kb/module/test/data'
##################################
##################################



mock_dfu = create_autospec(DataFileUtil, instance=True) # get_objects, save_objects (for AmpliconSet workflow)
mock_kbr = create_autospec(KBaseReport, instance=True) # create_extended_report
mock_ws = create_autospec(Workspace, instance=True) # get_object_info3
mock_run_check = create_autospec(run_check) # avoid lengthy runs

""" Trying to tease these out from get_mock_dfu
#####
#####
#####
def mock_dfu_save_objects(params):
    params_str = str(params)
    if len(params_str) > 100: params_str = params_str[:100] + ' ...'
    logging.info('Mocking `dfu.save_objects` with `params=%s`' % params_str)
    

#####
#####
#####
def mock_dfu_get_objects(params):
    logging.info('Mocking `dfu.get_objects` with `params=%s`' % str(params))

    upa = params['object_refs'][0]
    flnm = {
        enigma50by30_noAttrMaps_noSampleSet = 'AmpliconMatrix.json',
        enigma50by30_noSampleSet = 'AmpliconMatrix.json',
        enigma50by30 = 'AmpliconMatrix.json',
        }[upa]
    flpth = os.path.join(testData_dir, 'by_dataset_input', dataset, 'get_objects', flnm)

    with open(flpth) as f:
        obj = json.load(f)

    return obj
"""

#####
#####
#####
def get_mock_dfu(dataset):
    '''
    Pass in `dataset=None` to get reset mock_dfu
    It's a global so it can be edited in test-method
    '''
    # reset
    mock_dfu.reset_mock(return_value=True, side_effect=True)
    
    # return clean `mock_dfu`
    if dataset == None:
        return mock_dfu

    ##
    ## mock `save_objects`
    def mock_dfu_save_objects(params):
        params_str = str(params)
        if len(params_str) > 100: params_str = params_str[:100] + ' ...'
        logging.info('Mocking `dfu.save_objects` with `params=%s`' % params_str)

        return [['-1111', 1, 2, 3, '-1111', 5, '-1111']] # UPA made from pos 6/0/4
    
    mock_dfu.save_objects.side_effect = mock_dfu_save_objects

    ##
    ## mock `get_objects`
    def mock_dfu_get_objects(params):
        logging.info('Mocking `dfu.get_objects` with `params=%s`' % str(params))

        upa = params['object_refs'][0]
        flnm = {
                # TODO
            }[upa]
        flpth = os.path.join(testData_dir, 'by_dataset_input', dataset, 'get_objects', flnm)

        with open(flpth) as f:
            obj = json.load(f)

        return obj

    mock_dfu.get_objects.side_effect = mock_dfu_get_objects

    return mock_dfu



#####
#####
#####
def get_mock_kbr(dummy_dataset=None): # allow dummy param since the other mocks take `dataset` arg
    # reset
    mock_kbr.reset_mock(return_value=True, side_effect=True)

    # mock `create_extended_report`
    def mock_create_extended_report(params):
        logging.info('Mocking `kbr.create_extended_report`')
        return {
            'name': 'kbr mock name',
            'ref': 'kbr mock ref',
        }

    mock_kbr.create_extended_report.side_effect = mock_create_extended_report
    
    return mock_kbr



#####
#####
#####
def get_mock_run_check(dataset):
    # reset
    mock_run_check.reset_mock()

    # side effect
    def mock_run_check_(cmd):
        logging.info('Mocking running cmd `%s`' % cmd)
       
        rmtree(Var.return_dir)
        copytree(os.path.join(testData_dir, 'by_dataset_input', dataset, 'return'), Var.return_dir)

    mock_run_check.side_effect = mock_run_check_

    return mock_run_check

