from auxil.harness import KBFaprotaxTest, skip_cond, patch_, patch_dict_
from auxil.upa import * # defined UPAs
from kb_faprotax.impl import error


class GenomeSetWorkflowTest(KBFaprotaxTest):

    # doesn't need mocking - fast
    def test_GenomeSet_input(self):
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': refseq,
            }
        )

    # doesn't need mocking - fast
    def test_dup_GenomeSet(self):
        ret = self.serviceImpl.run_FAPROTAX(
            self.ctx, {
                **self.params_ws,
                'input_upa': refseq_dup,
            })

        self.assertTrue(len(Var.warnings) == 1) # TODO fix message



