from pytest import raises

from kb_faprotax.kb_faprotaxImpl import kb_faprotax
from kb_faprotax.kb_faprotaxServer import MethodContext
from kb_faprotax.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace

from kb_faprotax.util.dprint import dprint, where_am_i
from kb_faprotax.util.file import get_numbered_duplicate
from kb_faprotax.util.error import NonZeroReturnException
from kb_faprotax.util.workflow import run_check
from mock import * # mock business
from upa import * # upa library
import config




####################################################################################################
####################################################################################################
def test_run_check():
    with raises(NonZeroReturnException, match='`127`'):
        run_check('set -o pipefail && tmp |& tee tmp')

    with raises(NonZeroReturnException, match='`127`'):
        run_check('set -o pipefail && echo hi |& tmp')

    run_check('set -o pipefail && echo hi |& tee tmp') # run correctly


####################################################################################################
####################################################################################################
def test_get_numbered_duplicate():                                                          
    # test numbering system                                                                     
    q = 'the_attr'                                                                              
    names = ['null']                                                                            
    assert get_numbered_duplicate(names, q) == 'the_attr'                                       
    names = ['the_attr']                                                                        
    assert get_numbered_duplicate(names, q) == 'the_attr (1)'                                   
    names = ['the_attr', 'the_attr (1)', 'the_attr (2)']                                        
    assert get_numbered_duplicate(names, q) == 'the_attr (3)', get_numbered_duplicate(names, q) 
    names = ['the_attr (1)', 'the_attr (2)']                                                    
    assert get_numbered_duplicate(names, q) == 'the_attr'                                       
    names = ['the_attr', 'the_attr (1)', 'the_attr (5)']                                        
    assert get_numbered_duplicate(names, q) == 'the_attr (2)'                                   
    names = ['the_attr (0)', 'the_attr (1)', 'the_attr (2)']                                    
    assert get_numbered_duplicate(names, q) == 'the_attr'                                       
    names = ['the_attr (-1)', 'the_attr (1)', 'the_attr (2)']                                   
    assert get_numbered_duplicate(names, q) == 'the_attr'                                       
                                                                                                
    # test internal regexing                                                                    
    q = 'the[0-9]_attr'                                                                         
    names = [q]                                                                                 
    assert get_numbered_duplicate(names, q) == q + ' (1)'                                       
                                                                                                
   
