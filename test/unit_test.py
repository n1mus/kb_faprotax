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

from kb_faprotax.impl.kbase_obj import GenomeSet, Genome, AmpliconSet, AmpliconMatrix, AttributeMapping
from kb_faprotax.impl.error import * # exceptions, messages
from kb_faprotax.util.config import Var
from auxil.harness import KBFaprotaxTest, skip_cond, patch_, patch_dict_
from auxil.mock import * # custom mock logic


class KBFaprotaxUnitTest(KBFaprotaxTest):

    # doesn't need mocking - fast
    def test_parse_faprotax_functions(self):
        '''Test function used to parse FAPROTAX's predicted functions from its output'''
        from kb_faprotax.impl.workflow import parse_faprotax_functions

        flnm = 'groups2records_dense.tsv'
        flpth = os.path.join(testData_dir, 'by_dataset_input/refseq/return/FAPROTAX_output', flnm)
        
        r2g_d = parse_faprotax_functions(flpth)

        taxonomy = '\
cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Moraxellaceae; \
Acinetobacter; Acinetobacter calcoaceticus/baumannii complex; Acinetobacter pittii'
        functions = '\
aerobic_chemoheterotrophy,human_pathogens_septicemia,human_pathogens_pneumonia,human_pathogens_nosocomia,\
human_pathogens_all,animal_parasites_or_symbionts,aromatic_compound_degradation,chemoheterotrophy'

        self.assertTrue(r2g_d[taxonomy] == functions)


    # doesn't need mocking - fast
    def test_run_check(self):
        '''Test function that runs FAPROTAX'''
        from kb_faprotax.impl.workflow import run_check

        with self.assertRaises(NonZeroReturnException) as cm:
            run_check('set -o pipefail && ;s |& tee tmp')
            self.assertTrue('`2`') in str(cm.exception) # return code 2

        with self.assertRaises(NonZeroReturnException) as cm:
            run_check('set -o pipefail && tmp |& tee tmp')
            self.assertTrue('`127`') in str(cm.exception) # return code 127

        with self.assertRaises(NonZeroReturnException) as cm:
            run_check('set -o pipefail && echo hi |& tmp')
            self.assertTrue('`127`') in str(cm.exception) # return code 127

        run_check('set -o pipefail && echo hi |& tee tmp') # run correctly


 
    # TODO superficially test ctor
    def test_GenomeSet_methods(self):
        pass
       

    def test_Genome_methods(self):
        '''
        Genome methods are quite simple, 
        they are tested in other places
        '''
        pass


    # doesn't need mocking - fast
    def test_dummy_OTU_table_abundance(self):
        '''
        Test different dummy abundances in FAPROTAX input OTU table
        (The GenomeSet workflow uses dummy values in the input OTU table 
        just to map Genomes to FAPROTAX functions)
        '''
        from kb_faprotax.impl.workflow import run_check, parse_faprotax_functions
        from itertools import combinations

        run_dir = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(run_dir)

        gs = GenomeSet(refseq)

        dummy_vals = [1.0, 20, 'random']
        groups2records_dense_flpths = []

        # run FAPROTAX with different dummy values
        for val in dummy_vals:
            otu_table_flpth = os.path.join(run_dir, 'otu_table_' + str(val) + '.tsv')
            log_flpth = os.path.join(run_dir, 'log_' + str(val) + '.txt')
            
            gs.to_OTU_table(otu_table_flpth)
            groups2records_dense_flpth = os.path.join(run_dir, 'groups2records_dense_' + str(val) + '.tsv') 

            groups2records_dense_flpths.append(groups2records_dense_flpth)

            cmd = ' '.join([
                'set -o pipefail &&',
                Var.cmd_flpth,
                '--input_table', otu_table_flpth,
                '--input_groups_file', Var.db_flpth,
                '--out_groups2records_table_dense', groups2records_dense_flpth,
                '--row_names_are_in_column', 'taxonomy',
                '--verbose',
                '|& tee',
            ])

            with open(os.path.join(run_dir, 'cmd_' + str(val) + '.txt'), 'w') as f:
                f.write(cmd)

            run_check(cmd)

        # parse FAPROTAX output for each dummy value
        r2g_d_l = []
        for groups2records_dense_flpth in groups2records_dense_flpths:
            r2g_d = parse_faprotax_functions(groups2records_dense_flpth)
            r2g_d_l.append(r2g_d)

        # check all output is equal
        for combo in combinations(r2g_d_l, 2):
            self.assertTrue(combo[0] == combo[1])
 
    @patch.dict('kb_faprotax.impl.kbase_obj.Var', values={'dfu': get_mock_dfu('17770')})
    def test_AmpliconSet_methods(self):
        '''
        Superficial testing
        '''

        ##
        ## good input

        amp_set = AmpliconSet(the_17770)

        id_l = ['fffa52555f0d542613a26955a558d76d', 'fffd0e743b783e58fb70f6a5c5e41901']
        taxStr_l = ['\
D_0__Bacteria; D_1__Proteobacteria; D_2__Gammaproteobacteria; D_3__Steroidobacterales; D_4__Steroidobacteraceae; \
D_5__Steroidobacter; D_6__Steroidobacter:u; D_7__Steroidobacter:u; D_8__Steroidobacter:u; D_9__Steroidobacter:u; \
D_10__Steroidobacter:u', 'D_0__Bacteria; D_1__Proteobacteria; D_2__Deltaproteobacteria; D_3__Sva0485; D_4__Sva0485:u; \
D_5__Sva0485:u; D_6__Sva0485:u; D_7__Sva0485:u; D_8__Sva0485:u; D_9__Sva0485:u; D_10__Sva0485:u']

        self.assertTrue(amp_set.get_taxStr_l([]) == [])
        self.assertTrue(amp_set.get_taxStr_l(id_l) == taxStr_l) 

        tax2ids_d = amp_set._get_tax2ids_d()
        for id, taxStr in zip(id_l, taxStr_l):
            self.assertTrue(id in tax2ids_d[taxStr])

        ##
        ## no taxonomy field input

        mock_dfu.get_objects.side_effect = lambda params: { # this mock is a global
                'data': [
                    {
                        'data': {
                            'amplicon_matrix_ref': '-1/-2/-3',
                            'amplicons': {
                                'amplicon0': {
                                    'consensus_sequence': 'aaggcctt',
                                    'taxonomy': {
                                        'lineage': [], # TODO disallow?
                                    },
                                },
                                'amplicon1': {
                                    'consensus_sequence': 'aaggcctt',
                                    'taxonomy': {
                                        'lineage': ['this', 'is bogus', 'tax'], # TODO disallow?
                                    },
                                },
                                'amplicon2': {
                                    'consensus_sequence': 'aaggcctt',
                                    'taxonomy': {},
                                },
                            },
                        },
                        'info': [
                            -666,
                            'AmpliconSet_NoTaxonomy_Dummy_Name',
                            'Dummy.AmpliconSet-1.0',
                        ]
                    }
                ]
        }

        amp_set = AmpliconSet('du/mm/y')

        with self.assertRaises(NoTaxonomyException) as cm:
            amp_set.get_taxStr_l(['amplicon0', 'amplicon1', 'amplicon2'])
            self.assertEquals('`%s`' % 'amplicon2' in str(cm.exception)) # TODO should amplicon0 and amplion1 trigger exceptions?


    @patch.dict('kb_faprotax.impl.kbase_obj.Var', values={'dfu': get_mock_dfu('17770')})
    def test_AmpliconMatrix_methods(self):
        amp_set = AmpliconSet(the_17770)
        amp_mat = AmpliconMatrix(the_17770_AmpMat, amp_set)

        # superficially test `to_OTU_table`
        df1 = pd.read_csv(
            os.path.join(testData_dir, 'by_dataset_input/17770/return/otu_table.tsv'), sep='\t', index_col='taxonomy') 
        df2 = amp_mat.to_OTU_table()

        self.assertTrue(df1.columns.tolist() == df2.columns.tolist())
        self.assertTrue(df1.index.tolist() == df2.index.tolist()) # index is `taxonomy`
        self.assertTrue(df1['OTU_Id'].tolist() == df2['OTU_Id'].tolist()) # first column is `OTU_Id`

        data1 = df1.iloc[:,1:].values # rest of columns are data
        data2 = df2.iloc[:,1:].values

        self.assertTrue(np.all(np.abs(data1 - data2) < 10**-15)) # python floating point issue 
                                                                 # e.g., 0.1 + 0.1 + 0.1 => 0.30000000000000004


    @patch.dict('kb_faprotax.impl.kbase_obj.Var', values={'dfu': get_mock_dfu(None), 'warnings': []})
    def test_AttributeMapping_methods(self):
        # set `get_objects` to return something simpler and independent
        obj = {
            "data": [
                {
                    "data": {
                        "attributes": [
                            {
                                "attribute": "celestial body",
                                "source": "upload",
                            },
                            {
                                "attribute": "biome",
                                "source": "upload",
                            },
                        ],
                        "instances": {
                            "fffffffffffffffffffffffffffffffd": [
                                "Ganymede",
                                "tundra",
                            ],
                            "fffffffffffffffffffffffffffffffe": [
                                "Mars",
                                "rainforest",
                            ],
                            "ffffffffffffffffffffffffffffffff": [
                                "Saturn",
                                "chaparral",
                            ],
                        }
                    },
                    "info": [
                        -1,
                        "AttributeMapping_dummy_name",
                    ]
                }
            ]
        }
        mock_dfu.get_objects.side_effect = lambda params: obj

        id2attr_d = {'f' * 31 + 'd': 'black hole', 'f' * 31 + 'e': 'quasar', 'f' * 32: 'dark matter'}

        attr_map = AttributeMapping('-111/-1/-1')

        ind = attr_map.get_attribute_slot('celestial body', 'unit testing')
        attr_map.update_attribute(ind, id2attr_d)
        self.assertTrue(len(Var.warnings) == 0)
        
        ind = attr_map.get_attribute_slot('celestial body', 'upload')
        attr_map.update_attribute(ind, id2attr_d)
        self.assertTrue(len(Var.warnings) == 1)



