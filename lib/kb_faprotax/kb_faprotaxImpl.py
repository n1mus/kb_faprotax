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
from .util.varstash import Var, reset # `Var` holds globals, `reset` clears everything but `Var.debug`
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
    GIT_COMMIT_HASH = "9a891976114f6e5cb6d5500a66f9cb333a7f36b3"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

        self.Var = {
            'shared_folder': self.shared_folder,
            'callback_url': self.callback_url,
            'kbase_endpoint': config['kbase-endpoint'], # contains environment, for constructing Genome landing page url
            'dfu': DataFileUtil(self.callback_url),
            'ws': Workspace(config['workspace-url']),
            'kbr': KBaseReport(self.callback_url, service_ver='dev'),
            'db_flpth': '/kb/module/data/FAPROTAX.txt', # curated database file for FAPROTAX
            'cmd_flpth': '/opt/FAPROTAX_1.2.1/collapse_table.py', # FAPROTAX executable
            'template_flpth': '/kb/module/ui/template/edge_data.tt',
        }


        #END_CONSTRUCTOR


    def faprotax(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN faprotax
    

        #
        ##
        ### set up globals `Var`
        ####
        #####

        reset(Var) # reset globals for this API method run

        Var.update({ 
            **self.Var,
            'params': params,
            'run_dir': os.path.join(self.Var['shared_folder'], str(uuid.uuid4())),
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

        dprint('oi', run=locals())
        
        if oi[2].startswith('KBaseSearch.GenomeSet'):
            return do_GenomeSet_workflow()

        elif oi[2].startswith('KBaseExperiments.AmpliconSet'):
            return do_AmpliconSet_workflow()

        else:
            raise Exception('Unknown type `%s` for `input_upa`' % oi[2])




        #END faprotax

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method faprotax return value ' +
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
