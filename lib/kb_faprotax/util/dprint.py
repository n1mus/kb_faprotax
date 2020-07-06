import functools
import json
import subprocess
import sys
import os
import time
import inspect
import time as time_

from .varstash import Var


subproc_run = functools.partial(
    subprocess.run, stdout=sys.stdout, stderr=sys.stderr, shell=True, executable='/bin/bash')

TAG_WIDTH = 80
MAX_LINES = 70


def dprint(*args, run=False, where=False, time=False, max_lines=MAX_LINES, exit=False, subproc_run_kwargs={}, print_kwargs={}):
    """
    For debug printing
    Also executes shell/python commands, printing the command and the outcome

    Input:
        args - strings, which can be evaluated with Bash or as python
        run - `shell` or `cli` if the `args` are for Bash, or a namespace dictionary if `args` are 
            python code
        where - include information (file, function, line) about the calling stack frame
        time - include how long it took to process each of `args`
        max_lines - limit how many lines to print (this will be json format)
        exit - exit after printing (for debugging)
    """

    if not Var.debug:
        return

    print = functools.partial(__builtins__['print'], **print_kwargs)

    def print_format(arg):
        if isinstance(arg, (dict, list)):
            arg_json = json.dumps(arg, indent=3, default=str)
            if max_lines != None and max_lines != False and arg_json.count('\n') > max_lines:
                arg_json = '\n'.join(arg_json.split('\n')[0:max_lines] + ['...'])
            print(arg_json)
        else:
            print(arg)

    print('#' * TAG_WIDTH)

    if where:
        last_frame = inspect.stack()[1]
        print("(file `%s`)\n(func `%s`)\n(line `%d`)" % (os.path.basename(last_frame[1]), last_frame[3], last_frame[2]))
    
    for arg in args:
        if time:
            t0 = time_.time()
        if run:
            print('>> ' + arg)
            if run in ['cli', 'shell']:
                completed_proc = subproc_run(arg, **subproc_run_kwargs)
                retcode = completed_proc.returncode
            elif isinstance(run, dict):
                print_format(eval(arg, run))
            else:
                assert False
        else:
            print_format(arg)
        if time:
            t = time_.time() - t0
            print('[%fs]' % t)
    
    print('-' * TAG_WIDTH)

    if exit == True:
        sys.exit(0)
    
    # return last retcode
    if run and run in ['cli', 'shell']:
        return retcode


def where_am_i(f):
    '''Decorator'''
    def f_new(*args, **kwargs):
        dprint("where am i?\n(func `%s`)" % (f.__qualname__), where=False)
        f(*args, **kwargs)
    return f_new




