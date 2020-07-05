from dotmap import DotMap


_config = DotMap({
    'debug': True, # toggle for global debugging behavior
    'db_flpth': '/kb/module/data/FAPROTAX.txt', # curated database file for FAPROTAX
    'cmd_flpth': '/opt/FAPROTAX_1.2.1/collapse_table.py', # FAPROTAX executable
})


def reset(dm: DotMap=None):
    return DotMap(_config)

Var = reset()
