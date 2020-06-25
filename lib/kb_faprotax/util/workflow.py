import logging
import os
import sys
import uuid
import subprocess
import functools
import pandas as pd
import numpy as np
import json
import shutil

from .dprint import dprint
from .varstash import Var
from .error import *
from .message import *
from .kbase_obj import AmpliconSet, AmpliconMatrix, AttributeMapping, GenomeSet, Genome







####################################################################################################
def run_check(cmd):
    '''
    Wrap tool-running method for patching
    '''
    logging.info(f'Running FAPROTAX via command `{cmd}`')
    completed_proc = subprocess.run(cmd, shell=True, executable='/bin/bash', stdout=sys.stdout, stderr=sys.stdout)

    if completed_proc.returncode != 0:
        msg = (
"FAPROTAX command `%s` returned with non-zero return code `%d`. Please check logs for more details")
        raise NonZeroReturnException(msg)



####################################################################################################
def parse_faprotax_traits(groups2records_table_dense_flpth, dlm=',') -> dict:
    '''
    Input: filepath for groups2records_dense.tsv
    Output: dict map from taxonomy to predicted functions
    '''
    g2r_df = pd.read_csv(groups2records_table_dense_flpth, sep='\t', comment='#')
    g2r_df = g2r_df.fillna('').drop_duplicates().set_index('record')

    r2g_d = g2r_df.to_dict(orient='index')
    r2g_d = {record: r2g_d[record]['group'] for record in r2g_d}
    if dlm != ',': 
        r2g_d = {record: groups.replace(',', dlm) for record, groups in r2g_d.items()}

    return r2g_d




####################################################################################################
####################################################################################################
def do_GenomeSet_workflow():



    #
    ##
    ### kbase obj
    ####
    #####

    logging.info('Loading GenomeSet or Genome data')

    gs = GenomeSet(Var.params['input_upa'])
    


    #
    ##
    ### params
    ####
    #####


    # the excess files are copied from AmpliconSet
    # and not necessary
    # just keep them for debugging purposes

    otu_table_flpth = os.path.join(Var.return_dir, 'OTU_table_1s.tsv')
    gs.to_OTU_table(otu_table_flpth)

    log_flpth = os.path.join(Var.return_dir, 'log.txt')
    cmd_flpth = os.path.join(Var.return_dir, 'cmd.txt')

    out_dir = os.path.join(Var.return_dir, 'faprotax_output') # for output files
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
        '--input_table', otu_table_flpth,
        '--input_groups_file', Var.db_flpth,
        '--out_collapsed', func_table_flpth,
        '--out_report', report_flpth,
        '--out_sub_tables_dir', sub_tables_dir,
        '--out_groups2records_table', groups2records_table_flpth,
        '--out_groups2records_table_dense', groups2records_table_dense_flpth,
        '--out_group_overlaps', group_overlaps_flpth,
        '--out_group_definitions_used', group_definitions_used_flpth,
        '--row_names_are_in_column OTU',
        '--verbose',
        '|& tee', log_flpth
    ])

    with open(cmd_flpth, 'w') as f:
        f.write(cmd)



    #
    ##
    ### cmd
    ####
    #####

    run_check(cmd)





    #
    ##
    ### results
    ####
    #####


    taxStr_2_traits_d = parse_faprotax_traits(groups2records_table_dense_flpth) # parse FAPROTAX results

    gs.df['functions'] = gs.df.apply(lambda row: taxStr_2_traits_d.get(row['taxonomy'], np.nan), axis=1) # stitch FAPROTAX results onto GenomeSet df


    dprint('gs.df', run=locals())
    
    #
    ##
    ### prep df for DataTables
    ####
    #####

    df = gs.df
    df['genome name'] = df.apply(lambda row: '<a href="%s" target="_blank">%s</a>' % (row['url'], row['name']), axis=1) # add column of Genome name linking to landing page
    df = df[['genome name', 'taxonomy', 'functions']] # filter to final columns
    df.columns = ['Genome Name', 'Taxonomy', 'FAPROTAX Functions'] # capitalize properly
    columns = [{'title': col} for col in df.columns.tolist()] # format for DataTables
    lines = df.values.tolist() # format for DataTables
    


    #
    ##
    ### report
    ####
    #####


    # prepare template
    template_flpth = os.path.join(Var.run_dir, os.path.basename(Var.template_flpth))
    shutil.copyfile(Var.template_flpth, template_flpth)
   
    tmpl_data = {
            'page_title': 'FAPROTAX Functions for GenomeSet: %s' % gs.name,
            'data_array': lines,
            'cols': columns,
    }

    html_links = [{
        'name': 'report.html',
        'template': {
            'template_file': template_flpth,
            'template_data_json': json.dumps(tmpl_data),
        },
        'description': 'Table of Genomes with assigned FAPROTAX traits',
    }]


    report_output = Var.kbr.create_extended_report({
        'warnings': Var.warnings,
        'html_links': html_links,
        'direct_html_link_index': 0,
        'report_object_name': 'kb_faprotax_report',
        'workspace_id': Var.params['workspace_id'],
           
    })

    output = {
        'report_name': report_output['name'],
        'report_ref': report_output['ref'],
    }

    return [output]






    



####################################################################################################
####################################################################################################
def do_AmpliconSet_workflow():

    
    #
    ##
    ### kbase obj
    ####
    #####

    logging.info('Loading AmpliconSet and AmpliconMatrix')

    amp_set = AmpliconSet(Var.params['input_upa'])
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


    log_flpth = os.path.join(Var.return_dir, 'log.txt')
    cmd_flpth = os.path.join(Var.return_dir, 'cmd.txt')

    out_dir = os.path.join(Var.return_dir, 'faprotax_output')
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
        'set -o pipefail &&',
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
        '--verbose',
        '|& tee', log_flpth
        ])

    with open(cmd_flpth, 'w') as f:
        f.write(cmd)



    #
    ##
    ### run
    ####
    #####



    run_check(cmd)


    #
    ##
    ### update AttributeMap, AmpliconMatrix, AmpliconSet
    ####
    #####


    ### just don't patch running for now?
    #if params.get('skip_run'):
    #    groups2records_table_dense_flpth = '/kb/module/test/data/faprotax_output/groups2records_dense.tsv'

    objects_created = []
    if amp_mat.row_attrmap_upa:

        attribute = 'FAPROTAX Traits'
        source = 'kb_faprotax/faprotax'

        taxStr_2_traits_d = parse_faprotax_traits(groups2records_table_dense_flpth)
        ind = row_attrmap.add_attribute_slot(attribute, source)
        row_attrmap.update_attribute(ind, taxStr_2_traits_d, amp_set) # AmpliconSet is the rosetta stone between taxStr and traits (supplies id)
        row_attrmap_upa_new = row_attrmap.save()

        amp_mat.update_row_attributemapping_ref(row_attrmap_upa_new)
        amp_mat_upa_new = amp_mat.save()

        amp_set.update_amplicon_matrix_ref(amp_mat_upa_new)
        amp_set_upa_new = amp_set.save(name=Var.params.get('output_amplicon_set_name'))
        
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
        'path': Var.return_dir, 
        'name': 'faprotax_results.zip',
        'description': 'Input, output, logs to FAPROTAX run'
        }]


    params_report = {
        'warnings': Var.warnings,
        'objects_created': objects_created,
        'file_links': file_links,
        'report_object_name': 'kb_faprotax_report',
        'workspace_id': Var.params['workspace_id'],
        }


    report_output = Var.kbr.create_extended_report(params_report)

    output = {
        'report_name': report_output['name'],
        'report_ref': report_output['ref'],
    }

    return [output]








