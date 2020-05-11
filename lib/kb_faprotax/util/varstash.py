class Var:

    DEBUG = True # toggle for global debugging behavior
    
    @classmethod
    def setup(cls):
        cls.debug = cls.DEBUG
        cls._originals = list(cls.__dict__.keys())

    @classmethod
    def update(cls, d: dict):
        for attr_name, attr in d.items():
            setattr(cls, attr_name, attr)
    
    @classmethod
    def reset(cls):
        # TODO check values?
        for attr_name, attr in cls.__dict__.copy():
            if attr_name not in cls.originals:
                delattr(cls, attr_name)


Var.setup()
