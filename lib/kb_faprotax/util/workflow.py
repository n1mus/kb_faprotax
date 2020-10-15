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
from .kbase_obj import AmpliconMatrix, AttributeMapping, GenomeSet, Genome







####################################################################################################
def run_check(cmd):
    '''
    Wrap tool-running method for patching
    '''
    logging.info(f'Running FAPROTAX via command `{cmd}`')
    completed_proc = subprocess.run(cmd, shell=True, executable='/bin/bash', stdout=sys.stdout, stderr=sys.stdout)

    if completed_proc.returncode != 0:
        msg = (
"FAPROTAX command `%s` returned with non-zero return code `%d`. Please check logs for more details" % (cmd, completed_proc.returncode))
        raise NonZeroReturnException(msg)



####################################################################################################
def parse_tax2groups(groups2records_table_dense_flpth, dlm=',') -> dict:
    '''
    In FAPROTAX, a taxonomic path is also known as a 'record', and a
    function is also known as a 'group'.

    The input is a filepath for grouops2records_dense.tsv, which is one of
    the outputs from FAPROTAX. Basically, its index is the records, its columns
    are the groups, and its values are nonnegative floats

    Input: filepath for groups2records_dense.tsv
    Output: dict map from taxonomy to predicted functions
    '''
    
    """
    tax2groups_df = pd.read_csv(groups2records_table_dense_flpth, sep='\t', comment='#')
    tax2groups_df = tax2groups_df.fillna('').drop_duplicates().set_index('record')

    tax2groups_d = tax2groups_df.to_dict(orient='index')
    tax2groups_d = {record: tax2groups_d[record]['group'] for record in tax2groups_d}
    if dlm != ',': 
        tax2groups_d = {record: group.replace(',', dlm) for record, group in tax2groups_d.items()}

    return tax2groups_d
    """
    tax2groups_df = pd.read_csv(groups2records_table_dense_flpth, sep='\t', comment='#').fillna('')
    return tax2groups_df['record'].tolist(), tax2groups_df['group'].tolist()


####################################################################################################
def parse_faprotax_functions(groups2records_table_dense_flpth, dlm=',') -> dict: # TODO does same thing as parse_tax2groups
    '''
    In FAPROTAX, a taxonomic path is also known as a 'record', and a
    function is also known as a 'group'.
    The input is a filepath for grouops2records_dense.tsv, which is one of
    the outputs from FAPROTAX. Basically, its index is the records, its columns
    are the groups, and its values are nonnegative floats
    Input: filepath for groups2records_dense.tsv
    Output: dict map from taxonomy to predicted functions
    '''
    g2r_df = pd.read_csv(groups2records_table_dense_flpth, sep='\t', comment='#')
    g2r_df = g2r_df.fillna('').drop_duplicates().set_index('record')

    r2g_d = g2r_df.to_dict(orient='index')
    r2g_d = {record: r2g_d[record]['group'] for record in r2g_d}
    if dlm != ',': 
        r2g_d = {record: group.replace(',', dlm) for record, group in r2g_d.items()}

    return r2g_d

####################################################################################################
def map_groups2records_to_groups2ids(groups2records_table_flpth, groups2ids_table_flpth, id_l):

    tax2groups_df = pd.read_csv(groups2records_table_flpth, sep='\t', comment='#').fillna('')

    df = tax2groups_df.drop('record', axis=1)
    df.index = id_l
    df.index.name = 'id'

    df.to_csv(groups2ids_table_flpth, sep='\t', header=True, index=True)


    



####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
def do_AmpliconMatrix_workflow():

    
    #
    ##
    ### kbase obj
    ####
    #####

    amp_mat = AmpliconMatrix(Var.params['input_upa'])

    row_attr_map_upa = amp_mat.obj.get('row_attributemapping_ref')
    if row_attr_map_upa is None:
        msg = (
"Input AmpliconMatrix does not have a row AttributeMapping object to grab taxonomy from and assign traits to. "
"To upload a row AttributeMapping with taxonomy, try running attribute mapping import app. "
"To create a row AttributeMapping and assign it with taxonomy, try running kb_RDP_Classifier app first"
        )
        raise NoTaxonomyException(msg)

    row_attr_map = AttributeMapping(row_attr_map_upa)
    amp_mat.row_attr_map = row_attr_map



    #
    ##
    ### id2tax business
    ####
    #####


    ind, tax_attribute = row_attr_map.get_tax_attribute(
        regex_l = [
            'rdp classifier taxonomy',
            'parsed_user_taxonomy',
            'taxonomy',
        ]
    ) # detect attribute that holds tax

    if ind is None:
        raise NoTaxonomyException('Sorry no taxonomy was found!')

    # id_l and tax_l correspond to AmpliconMatrix rows
    id_l = amp_mat.obj['data']['row_ids']
    tax_l = row_attr_map.get_tax_l(ind, id_l)



    #
    ##
    ### params
    ####
    #####


    log_flpth = os.path.join(Var.return_dir, 'log.txt')
    cmd_flpth = os.path.join(Var.return_dir, 'cmd.txt')

    taxon_table_flpth = os.path.join(Var.return_dir, 'otu_table.tsv')
    amp_mat.to_OTU_table(tax_l, taxon_table_flpth)

    Var.out_dir = os.path.join(Var.return_dir, 'FAPROTAX_output')
    sub_tables_dir = os.path.join(Var.out_dir, 'sub_tables')

    os.mkdir(Var.out_dir)
    os.mkdir(sub_tables_dir)

    collapsed_func_table_flpth = os.path.join(Var.out_dir, 'collapsed_func_table.tsv')
    report_flpth = os.path.join(Var.out_dir, 'report.txt')
    groups2records_table_flpth = os.path.join(Var.out_dir, 'groups2records.tsv')
    groups2records_table_dense_flpth = os.path.join(Var.out_dir, 'groups2records_dense.tsv')
    group_overlaps_flpth = os.path.join(Var.out_dir, 'group_overlaps.tsv')
    group_definitions_used_flpth = os.path.join(Var.out_dir, 'group_definitions_used.txt')


    cmd = ' '.join([
        'set -o pipefail &&',
        Var.cmd_flpth,
        '--input_table', taxon_table_flpth,
        '--input_groups_file', Var.db_flpth,
        '--out_collapsed', collapsed_func_table_flpth,
        '--out_report', report_flpth,
        '--out_sub_tables_dir', sub_tables_dir,
        '--out_groups2records_table', groups2records_table_flpth,
        '--out_groups2records_table_dense', groups2records_table_dense_flpth,
        '--out_group_overlaps', group_overlaps_flpth,
        '--out_group_definitions_used', group_definitions_used_flpth,
        '--row_names_are_in_column', 'taxonomy',
        '--omit_columns', '1', # amplicon ID column
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
    ### update AttributeMapping, AmpliconMatrix
    ####
    #####



    attribute = 'FAPROTAX Functions'
    source = 'kb_faprotax/run_FAPROTAX'
    
    '''
    tax2groups = parse_tax2groups(groups2records_table_dense_flpth)
    id2groups = {id: tax for id, (tax, groups) in zip(id_l, tax2groups.items())}
    '''
    tax_order, groups_l = parse_tax2groups(groups2records_table_dense_flpth)
    id2groups = {id: group for id, group in zip(id_l, groups_l)}

    if Var.debug: # check FAPROTAX did not reorder the rows of OTU table
        assert len(tax_order) == len(tax_l), 'Lengths: %s vs %s' % (len(tax_order), len(tax_l))
        for i, (tax1, tax2) in enumerate(zip(tax_order, tax_l)):
            assert tax1 == tax2, 'index %d\n`%s`\n`%s`' % (i, tax1, tax2)
        assert tax_order == tax_l, '`%s`\n`%s`' % (tax_order, tax_l)

    ind = row_attr_map.add_attribute_slot(attribute, source)
    row_attr_map.update_attribute(ind, id2groups)
    row_attr_map_upa_new = row_attr_map.save()

    amp_mat.obj['row_attributemapping_ref'] = row_attr_map_upa_new
    amp_mat_upa_new = amp_mat.save(name=Var.params.get('output_amplicon_matrix_name'))

    Var.objects_created = [
        {'ref': row_attr_map_upa_new, 'description': 'Added or updated attribute `%s`' % attribute}, 
        {'ref': amp_mat_upa_new, 'description': 'Updated row AttributeMapping reference'},
    ]




    #
    ##
    ### make FunctionalProfiles
    ####
    #####

                 
    ### Organism FP ###

    # map groups2records back to groups2ids
    groups2ids_table_flpth = os.path.join(Var.run_dir, 'groups2ids.tsv')
    map_groups2records_to_groups2ids(groups2records_table_flpth, groups2ids_table_flpth, id_l)

    func_prof_amplicon_upa = Var.fpu.import_func_profile(dict(
        workspace_id=Var.params['workspace_id'],
        func_profile_obj_name='%s.FAPROTAX_groups2records' % amp_mat.name,
        original_matrix_ref=amp_mat_upa_new,
        profile_file_path=groups2ids_table_flpth,
        profile_type='amplicon',
        profile_category='organism',
        data_epistemology='predicted',
        epistemology_method='FAPROTAX',
        description='Amplicon vs FAPROTAX functions. In binary for association. Double inference through taxonomic assignment and FAPROTAX',
    ))['func_profile_ref']

    Var.objects_created.append(
        dict(ref=func_prof_amplicon_upa, description='Amplicon FunctionalProfile')
    )

    ### Community FP ###
    if amp_mat.obj.get('sample_set_ref') is None:
        msg = (
            'AmpliconMatrix with name %s has no SampleSet reference, '
            'and so making a community FunctionalProfile is prohibited'
        )
        logging.warning(msg)
        Var.warnings.append(msg)

    else:
        func_prof_sample_upa = Var.fpu.import_func_profile(dict(
            workspace_id=Var.params['workspace_id'],
            func_profile_obj_name='%s.FAPROTAX_func_table' % amp_mat.name,
            original_matrix_ref=amp_mat_upa_new, 
            profile_file_path=collapsed_func_table_flpth,
            profile_type='mg',
            profile_category='community',
            data_epistemology='predicted',
            epistemology_method='FAPROTAX',
            description='Metagenomic sample vs FAPROTAX functions'
        ))['func_profile_ref']


        Var.objects_created.append(
            dict(ref=func_prof_sample_upa, description='Sample FunctionalProfile')
        )
    


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
        'objects_created': Var.objects_created,
        'file_links': file_links,
        'report_object_name': 'kb_faprotax_report',
        'workspace_id': Var.params['workspace_id'],
        }


    report_output = Var.kbr.create_extended_report(params_report)

    output = {
        'report_name': report_output['name'],
        'report_ref': report_output['ref'],
    }

    dprint('output', run=locals())

    return output








####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
def do_GenomeSet_workflow():



    #
    ##
    ### kbase obj
    ####
    #####

    gs = GenomeSet(Var.params['input_upa'])
    


    #
    ##
    ### params
    ####
    #####


    # the excess files are copied from Amplicon workflow
    # and not are not necessary

    otu_table_flpth = os.path.join(Var.return_dir, 'otu_table.tsv')
    gs.to_OTU_table(otu_table_flpth)

    log_flpth = os.path.join(Var.return_dir, 'log.txt')
    cmd_flpth = os.path.join(Var.return_dir, 'cmd.txt')

    Var.out_dir = os.path.join(Var.return_dir, 'FAPROTAX_output') # for output files
    sub_tables_dir = os.path.join(Var.out_dir, 'sub_tables')

    os.mkdir(Var.out_dir)
    os.mkdir(sub_tables_dir)

    collapsed_func_table_flpth = os.path.join(Var.out_dir, 'collapsed_func_table.tsv')
    report_flpth = os.path.join(Var.out_dir, 'report.txt')
    groups2records_table_flpth = os.path.join(Var.out_dir, 'groups2records.tsv')
    groups2records_table_dense_flpth = os.path.join(Var.out_dir, 'groups2records_dense.tsv')
    group_overlaps_flpth = os.path.join(Var.out_dir, 'group_overlaps.tsv')
    group_definitions_used_flpth = os.path.join(Var.out_dir, 'group_definitions_used.txt')


    cmd = ' '.join([
        'set -o pipefail &&',
        Var.cmd_flpth,
        '--input_table', otu_table_flpth,
        '--input_groups_file', Var.db_flpth,
        '--out_collapsed', collapsed_func_table_flpth,
        '--out_report', report_flpth,
        '--out_sub_tables_dir', sub_tables_dir,
        '--out_groups2records_table', groups2records_table_flpth,
        '--out_groups2records_table_dense', groups2records_table_dense_flpth,
        '--out_group_overlaps', group_overlaps_flpth,
        '--out_group_definitions_used', group_definitions_used_flpth,
        '--row_names_are_in_column', 'taxonomy',
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


    tax2groups = parse_tax2groups(groups2records_table_dense_flpth, dlm=', ') # parse FAPROTAX results

    gs.df['functions'] = gs.df.apply(lambda row: tax2groups.get(row['taxonomy'], np.nan), axis=1) # stitch FAPROTAX results onto GenomeSet df


    dprint('gs.df', run=locals())
    
    #
    ##
    ### prep df for DataTables
    ####
    #####

    df = gs.df
    df['genome name'] = df.apply(lambda row: '<a href="%s" target="_blank">%s</a>' % (row['url'], row['name']), axis=1) # add column of Genome name linking to landing page
    df = df[['genome name', 'taxonomy', 'functions']] # filter to final columns
    df.columns = ['Genome Workspace Name', 'Taxonomy', 'FAPROTAX Functions'] # capitalize properly
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

    return output

    



