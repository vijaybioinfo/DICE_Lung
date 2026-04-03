from .plink import PLINK


class LDFactory:
    METHODS = {"PLINK": PLINK}

    @classmethod
    def build(cls, name: str):
        return cls.METHODS[name]()