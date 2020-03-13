




class Var:

    @classmethod
    def update(cls, d: dict):

        for attr_name, attr in list(d.items()):
            setattr(cls, attr_name, attr)




