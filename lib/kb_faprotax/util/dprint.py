import functools
import pprint, json
import subprocess
import sys
import os
import time
import logging

MAX_LINES = 70
subprocess.run = functools.partial(
        subprocess.run, stdout=sys.stdout, stderr=sys.stderr, shell=True, executable='/bin/bash')

# TODO time, where
def dprint(*args, run=False, subproc_run_kwargs={}, print_kwargs={}):
    print = functools.partial(__builtins__['print'], **print_kwargs)

    def print_format(arg):
        if isinstance(arg, (dict, list)):
            arg_json = json.dumps(arg, indent=3, default=str)
            if arg_json.count('\n') > MAX_LINES:
                arg_json = '\n'.join(arg_json.split('\n')[0:MAX_LINES] + ['...'])
            print(arg_json)
        else:
            print(arg, end=' ')

    print('##############################################################')
    for arg in args:
        if run:
            print('>> ' + arg)
            if run in ['cli', 'shell']:
                completed_proc = subprocess.run(arg, **subproc_run_kwargs)
                retcode = completed_proc.returncode
            elif isinstance(run, dict):
                print_format(eval(arg, run))
            else:
                assert False
        else:
            print_format(arg)
        print()
    print('--------------------------------------------------------------')
    # return last run
    if run and run in ['cli', 'shell']:
        return retcode

# TODO
def where_am_i(f):
    '''
    Decorator
    '''
    dprint("where am i? in module " + globals()['__name__'] + " method " + f.__qualname__)
    return f
