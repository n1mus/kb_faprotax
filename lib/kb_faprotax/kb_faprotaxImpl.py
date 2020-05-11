# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
import uuid
import subprocess
import functools


from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil


from .util.kbase_obj import AmpliconSet, AmpliconMatrix, AttributeMapping
from .util.dprint import dprint
from .util.varstash import Var


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

        Var.update({
            'shared_folder': self.shared_folder,
            'callback_url': self.callback_url,
            'dfu': DataFileUtil(self.callback_url),
            'db_flpth': '/kb/module/data/FAPROTAX.txt',
            'cmd_flpth': '/opt/FAPROTAX_1.2.1/collapse_table.py',
            })



        #END_CONSTRUCTOR
        pass


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
       


        Var.update({ # TODO reset this beginning of API-method run
            'params': params,
            'run_dir': os.path.join(Var.shared_folder, str(uuid.uuid4())),
            'warnings': [],
            })

        os.mkdir(Var.run_dir)

        dprint(params)



        #####
        ##### kbase obj
        #####

        logging.info('Loading AmpliconSet and AmpliconMatrix')

        amp_set = AmpliconSet(params['amplicon_set_upa'])
        amp_mat = AmpliconMatrix(amp_set.get_amplicon_matrix_upa(), amp_set) 
        if amp_mat.row_attrmap_upa:
            logging.info('Loading row AttributeMapping')
            row_attrmap = AttributeMapping(amp_mat.row_attrmap_upa)
        else:
            msg = (
"Input AmpliconSet's associated AmpliconMatrix does not have a row AttributeMapping object to assign traits to. "
"To create a row AttributeMapping, try running the the AttributeMapping uploader or kb_RDP_Classifier first")
            logging.warning(msg)
            Var.warnings.append(msg)






        #
        ##
        ### params
        ####
        #####


        log_flpth = os.path.join(Var.run_dir, 'log.txt')

        out_dir = os.path.join(Var.run_dir, 'faprotax_output')
        sub_tables_dir = os.path.join(out_dir, 'sub_tables')

        os.mkdir(out_dir)
        os.mkdir(sub_tables_dir)

        func_table_flpth = os.path.join(out_dir, 'func_table.tsv')
        report_flpth = os.path.join(out_dir, 'report.txt')
        groups2records_table_flpth = os.path.join(out_dir, 'groups2records.tsv')
        groups2records_table_dense_flpth = os.path.join(out_dir, 'groups2records_dense.tsv')
        group_overlaps_flpth = os.path.join(out_dir, 'group_overlaps.tsv')
        group_definitions_used_flpth = os.path.join(out_dir, 'group_definitions_used.txt')


        cmd = ' '.join([
            Var.cmd_flpth,
            '--input_table', amp_mat.taxon_table_flpth,
            '--input_groups_file', Var.db_flpth,
            '--out_collapsed', func_table_flpth,
            '--out_report', report_flpth,
            '--out_sub_tables_dir', sub_tables_dir,
            '--out_groups2records_table', groups2records_table_flpth,
            '--out_groups2records_table_dense', groups2records_table_dense_flpth,
            '--out_group_overlaps', group_overlaps_flpth,
            '--out_group_definitions_used', group_definitions_used_flpth,
            '--row_names_are_in_column', 'taxonomy',
            '--omit_columns', '1',
            '--group_leftovers_as', 'OTUs_no_func_assign',
            '--verbose',
            '|& tee', log_flpth
            ])





        #
        ##
        ### run
        ####
        #####


        if params.get('skip_run'):
            logging.info('Skip run')


        else:

            logging.info(f'Running FAPROTAX via command `{cmd}`')

            completed_proc = subprocess.run(cmd, shell=True, executable='/bin/bash', stdout=sys.stdout, stderr=sys.stdout)

            if completed_proc.returncode != 0:
                msg = (
"FAPROTAX command `%s` returned with non-zero return code `%d`. Please check logs for more details")
                raise NonZeroReturnException(msg)


        #
        ##
        ### update AttributeMap, AmpliconMatrix, AmpliconSet
        ####
        #####


        if params.get('skip_run'):
            groups2records_table_dense_flpth = '/kb/module/test/data/faprotax_output/groups2records_dense.tsv'

        objects_created = []
        if amp_mat.row_attrmap_upa:

            attribute = 'FAPROTAX Traits'
            source = 'kb_faprotax/faprotax'

            taxStr_2_traits_d = row_attrmap.parse_faprotax_traits(groups2records_table_dense_flpth)
            ind = row_attrmap.add_attribute_slot(attribute, source)
            row_attrmap.update_attribute(ind, taxStr_2_traits_d, amp_set) # AmpliconSet is the rosetta stone between taxStr and traits (supplies id)
            row_attrmap_upa_new = row_attrmap.save()

            amp_mat.update_row_attributemapping_ref(row_attrmap_upa_new)
            amp_mat_upa_new = amp_mat.save()

            amp_set.update_amplicon_matrix_ref(amp_mat_upa_new)
            amp_set_upa_new = amp_set.save(name=params.get('output_amplicon_set_name'))
            
            objects_created = [
                {'ref': row_attrmap_upa_new, 'description': 'Added or updated attribute `%s`' % attribute}, 
                {'ref': amp_mat_upa_new, 'description': 'Updated row AttributeMapping reference'},
                {'ref': amp_set_upa_new, 'description': 'Updated AmpliconMatrix reference'},
            ]



        #
        ##
        ### return 
        ####
        #####


        file_links = [{
            'path': Var.run_dir, 
            'name': 'faprotax_results.zip',
            'description': 'Input, output'
            }]


        params_report = {
            'warnings': Var.warnings,
            'objects_created': objects_created,
            'file_links': file_links,
            'report_object_name': 'kb_faprotax_report',
            'workspace_name': params['workspace_name'],
            }


        if Var.debug and params.get('skip_kbReport'):
            output = {}

        else:
            kbr = KBaseReport(self.callback_url)
            report_output = kbr.create_extended_report(params_report)

            output = {
                'report_name': report_output['name'],
                'report_ref': report_output['ref'],
            }


        if Var.debug and params.get('return_test_info'):
            return {
                **params_report,
            }

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
