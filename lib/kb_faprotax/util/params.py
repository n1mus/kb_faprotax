import json

from .dprint import dprint
from .varstash import Var





#####
#####
class Params:


    DEFAULTS = {
            'output_amplicon_matrix_name': None,
    }

    def __init__(self, params):

        ## Validation
        self._validate(params)

        
        ## Custom transformations to internal state ##
        ## This is kind of silly ##
        ## But the code is written before I figure out how the narrative passes things ##

        if params.get('output_amplicon_matrix_name') == '':
            params['output_amplicon_matrix_name'] = None # treat empty string as null case since ui only returns strings for string type
        if type(params.get('tax_field')) is list:
            params['tax_field'] = params['tax_field'][0]
       
        ##
        self.params = params


    def _validate(self, params):
        '''
        None of params go directly to shell
        '''
        # TODO

        VALID = [
            'input_upa',
            'tax_field',
            'output_amplicon_matrix_name',
            #---
            'workspace_id',
            'workspace_name',
        ]

        for p in params:
            if p not in VALID:
                raise Exception(p)



    def __getitem__(self, key):
        '''
        For required params (e.g., input UPAs, workspace stuff)
        Should not use this for default-backed params
        as those can be left off params
        so use `getd` for those
        '''
        return self.params[key]


    def getd(self, key):
        '''
        For default-backed params (e.g., tunable numbers)
        Like `get`
        Return the user-supplied value, or the default value if none was supplied
        '''
        if key not in self.DEFAULTS:
            raise Exception('`params.getd(x)` only applicable to params with defaults')

        return self.params.get(key, self.DEFAULTS[key])


    def __repr__(self) -> str:
        return 'Wrapper for params:\n%s' % (json.dumps(self.params, indent=4))


