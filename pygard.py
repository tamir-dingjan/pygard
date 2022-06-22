from collections import defaultdict
import numpy as np

class Micelle():
    """
    A single micelle in the simulation. Tracks the composition of lipids within the micelle,
    and contains methods to add lipids and to split into daughter micelles.
    """
    
    def __init__(self) -> None:
        self.composition = defaultdict(int)
        pass

    def split(self):
        """
        Split the micelle into daughter micelles.
        """


class Soup():
    """
    The entire system, affectionatly titled the 'soup'. This object keeps track of the 
    total number of distinct amphiphile types available, the intrinsic affinity matrix, and
    the enhancement matrix.
    """

    def __init__(self) -> None:
        self.Ng = 100
        self.K = np.ones(shape=(self.Ng, self.Ng))
        self.B = np.ones(shape=(self.Ng, self.Ng))
