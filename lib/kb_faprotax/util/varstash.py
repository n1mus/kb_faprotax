from dotmap import DotMap


DEBUG = True # toggle for global debugging behavior

Var = DotMap(debug=DEBUG)


def reset(dm: DotMap):
    dm.clear()
    dm.update(debug=DEBUG)
