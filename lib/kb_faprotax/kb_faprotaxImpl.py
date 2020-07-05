# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
import uuid
import subprocess
import functools
from dotmap import DotMap


from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


from .util.kbase_obj import AmpliconSet, AmpliconMatrix, AttributeMapping
from .util.dprint import dprint
from .util.varstash import Var, reset # `Var` holds globals, `reset` clears everything but config stuff
from .util.workflow import do_AmpliconSet_workflow, do_GenomeSet_workflow


#END_HEADER


class kb_faprotax:
    '''
    Module Name:
    kb_faprotax

    Module Description:
    A KBase module: kb_faprotax
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/n1mus/kb_faprotax"
    GIT_COMMIT_HASH = "24efa7938f23c8343cac4a8accc754e2e625165b"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR

        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s', level=logging.INFO)

        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.workspace_url = config['workspace-url']

        self.Var = { # carry over into globals `Var`, regardless of resetting, for all API-method runs
            'kbase_endpoint': config['kbase-endpoint'], # contains environment, for constructing Genome landing page url
            'template_flpth': '/kb/module/ui/template/edge_data.tt',
        }
        

        #END_CONSTRUCTOR
        pass


    def run_FAPROTAX(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_FAPROTAX
    

        #
        ##
        ### set up globals `Var`
        ####
        #####

        reset(Var) # reset globals for this API method run

        Var.update({ 
            **self.Var,
            'ws': Workspace(self.workspace_url),
            'dfu': DataFileUtil(self.callback_url), # instantiate here so within runtime of @patch
            'kbr': KBaseReport(self.callback_url, service_ver='dev'), # instantiate here so within runtime of @patch 
            'params': params,
            'run_dir': os.path.join(self.shared_folder, str(uuid.uuid4())),
            'warnings': [],
        })

        os.mkdir(Var.run_dir)

        Var.update({
            'return_dir': os.path.join(Var.run_dir, 'return'),
        })

        os.mkdir(Var.return_dir)



        #
        ##
        ### detect input type
        ####
        #####


        oi = Var.ws.get_object_info3({'objects': [{'ref': params['input_upa']}]})['infos'][0]

        if oi[2].startswith('KBaseSearch.GenomeSet'):
            return do_GenomeSet_workflow()

        elif oi[2].startswith('KBaseExperiments.AmpliconSet'):
            return do_AmpliconSet_workflow()

        else:
            raise Exception('Unknown type `%s` for `input_upa`' % oi[2])




        #END run_FAPROTAX

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_FAPROTAX return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
