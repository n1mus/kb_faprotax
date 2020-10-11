from dotmap import DotMap

'''
App globals and config
'''

_config = DotMap({
    'debug': True, # toggle for global debugging behavior
    'db_flpth': '/kb/module/data/FAPROTAX.txt', # curated database file for FAPROTAX
    'cmd_flpth': '/opt/FAPROTAX_1.2.1/collapse_table.py', # FAPROTAX executable
})

Var = DotMap(_config) # globals for app

def reset_Var():
    Var.clear()
    Var.update(_config)
