from auxil.harness import KBFaprotaxTest, skip_cond, patch_, patch_dict_
from auxil.mock import * # custom mocking
from unittest.mock import patch




class AmpliconSetWorkflowTest(KBFaprotaxTest):
    '''
    Because of the control flow of testing,
    there are two main strategies here for patching
    (1) patch *Util/Workspace class 
        with constructor-like 
        in kb_faprotaxImpl,
        good for integration tests
    (2) patch.dict Var.dfu or Var.ws 
        with function side-effect 
        in the module it is used,
        good for unit tests
    '''

    # TODO mock dfu, run_check
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *args, **kwargs: get_mock_kbr())
    def test_AmpliconSet_no_row_AttributeMapping(self):
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': secret_wRDP,
            })
        
    # TODO mock run_check
    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *args: get_mock_dfu('secret'))
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *args, **kwargs: get_mock_kbr())
    def test_AmpliconSet_no_taxonomy_no_row_AttributeMapping(self):
        with self.assertRaises(NoTaxonomyException) as cm:
            ret = self.serviceImpl.run_FAPROTAX(
                self.ctx, {
                    **self.params_ws,
                    'input_upa': secret,
                })
            

    @patch_('kb_faprotax.kb_faprotaxImpl.DataFileUtil', new=lambda *args: get_mock_dfu('17770'))
    @patch_('kb_faprotax.util.workflow.run_check', new=get_mock_run_check('17770'))
    @patch_('kb_faprotax.kb_faprotaxImpl.KBaseReport', new=lambda *args, **kwargs: get_mock_kbr())
    def test_AmpliconSet_input_against_reference(self):
        '''
        Check results against answers if full pipeline is run (i.e., no mocks)
        
        Some things that have been hardcoded:
        * `source` is `FAPROTAX Functions`
        * FAPROTAX functions delimiter is `','`
        '''
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': the_17770,
                }
            )

        from kb_faprotax.util.workflow import run_check

        # only check results against reference
        # if full pipeline has been run
        if (isinstance(Var.dfu, unittest.mock.NonCallableMagicMock) or
            isinstance(Var.kbr, unittest.mock.NonCallableMagicMock) or
            hasattr(run_check, 'side_effect')):
            return

        logging.info('Comparing traits in AttributeMapping to answers')

        # hardcoded
        # should match in workflow.py
        source = 'FAPROTAX Functions'
        dlm = ','

        # answers to compare against
        def parse_answers_file(answers_flpth='/kb/module/test/data/reference/OTUMetaData_reduced.tsv'):
            answers_df = pd.read_csv(
                answers_flpth, sep='\t', header=0, index_col='#OTU ID', usecols=['#OTU ID', 'FAPROTAX Traits']).fillna('')
            answers_d = answers_df.to_dict(orient='index')
            answers_d = {key: value['FAPROTAX Traits'] for key, value in answers_d.items()}
            return answers_d

        # load KBaseReport
        attrmap_upa_new = Var.dfu.get_objects({
            'object_refs': [ret[0]['report_ref']]
        })['data'][0]['data']['objects_created'][0]['ref']

        # load AttributeMapping
        row_attrmap = AttributeMapping(attrmap_upa_new)
        instances_d = row_attrmap.obj['instances']
        attribute_d_l = row_attrmap.obj['attributes']

        # find index in attribute list
        for i, attribute_d in enumerate(attribute_d_l):
            if attribute_d['attribute'] == source:
                ind = i

        # id to attribute
        res_id2attr_d = {id: attr_l[ind] for id, attr_l in instances_d.items()}
        ans_id2attr_d = self.parse_answers_file()

        # check ids equal set
        res_ids = sorted(list(res_id2attr_d.keys()))
        ans_ids = sorted(list(ans_id2attr_d.keys()))
        self.assertTrue(res_ids == ans_ids)

        #
        html_l = []
        for id in res_ids:
            res = sorted(res_id2attr_d[id].split(dlm))
            ans = sorted(ans_id2attr_d[id].split(':'))

            # no dups
            self.assertTrue(len(set(res)) == len(res))

            if res != ans:
                res = set(res)
                ans = set(ans)
                
                ##
                ## check all extras are left-off parent functions
                self.assertTrue(ans.issubset(res))
                extras = list(res - ans)
                # check that extra is a substring of something in `ans`
                # TODO ascertain that left-off is parent of something in `ans`
                for extra in extras:
                    self.assertTrue(len([None for function in ans if extra in function]) > 0)

                # make a diff of extras in `res`
                html = '<p>' + ','.join([res_ if res_ in ans else '<b>' + res_ + '</b>' for res_ in res]) + '</p>'
                html_l.append(html)

        # dedup for visual ease
        html_l = list(set(html_l))

        with open(f'/kb/module/work/tmp/diff.html', 'w') as fp:
            fp.write('\n'.join(html_l))

                 


