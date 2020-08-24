__all__ = ['pysat_dineof']

for inst in __all__:
    exec("from pysatModels.instruments import {x}".format(x=inst))
