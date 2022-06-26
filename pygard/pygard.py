import numpy as np

class Micelle:
    """
    A single micelle in the simulation. Tracks the composition of lipids within the micelle,
    and contains methods to add lipids and to split into daughter micelles.
    """
    
    def __init__(self, N, Ng):
        # Store the N and Ng values for instantiating daughter micelles.
        self.N = N
        self.Ng = Ng

        # The micelle is a vector of length Ng, with N indices non-zero.
        self.composition = np.zeros(shape=(Ng), dtype=np.int64)

        # Increment a random amphiphile count
        self.rng = np.random.default_rng()
        for i in self.rng.integers(low=0, high=Ng, size=N):
            self.composition[i] += 1            


    def add_amphiphile(self,i):
        self.composition[i] += 1


    def split(self):
        """
        Split the micelle into daughter micelles.
        """
        # For each amphiphile, assign a random number of them to the first daughter
        # The remainder are assigned to the second daughter.
        compFirst = np.zeros(shape=self.composition.shape, dtype=np.int64)

        for i_i, i in enumerate(self.composition):
            if i == 0:
                continue
            else:
                transfer = self.rng.integers(low=0, high=i+1)
                compFirst[i_i] = transfer
            
        compSecond = self.composition - compFirst

        daughterFirst = Micelle(self.N, self.Ng)
        daughterFirst.composition = compFirst

        daughterSecond = Micelle(self.N, self.Ng)
        daughterSecond.composition = compSecond

        return (daughterFirst, daughterSecond)



class Soup:
    """
    The entire system, affectionatly titled the 'soup'. This object keeps track of the 
    total number of distinct amphiphile types available, the intrinsic affinity matrix, and
    the enhancement matrix.

    The enhancement matrix, B, is structured as a square matrix of size (Ng, Ng). It contains the 
    enhancement factors for all the amphiphiles in the system toward all other amphiphiles, 
    including self-to-self enhancements. The dimensions are arranged thusly:

    Vertical axis: soup amphiphiles
    Horizontal axis: micelle amphiphiles

     Soup amphiphile 1 ---> [[  x  x  x  ]
     Soup amphiphile 2 --->  [  x  x  x  ]
     Soup amphiphile 3 --->  [  x  x  x  ]]
                                ^  ^  ^
                                |  |  |
             Micelle amphiphile 1  |  |
                Micelle amphiphile 2  |
                   Micelle amphiphile 3

    When determining the influence of micelle composition on the enhancement of soup
    amphiphiles joining the micelle, the enhancement matrix is weighted by the micellar composition
    such that each row of the matrix reflects the enhancement toward a single amphiphile in the soup
    by every amphiphile in the micelle. 
    
    The intrinsic affinity matrix [also square, of size (Ng, Ng)] is then added to the
    weighted enhancement matrix.

    By taking the sum along the rows, the vector of per-amphiphile
    enhancement rates is produced. The vector element/s with the highest rate indicates the 
    amphiphile/s to be added to the growing micelle.
    """

    def __init__(self, N, Ng):
        self.N = N
        self.Ng = Ng
        self.K = np.ones(shape=(self.Ng, self.Ng))

        rng = np.random.default_rng()
        mu = 3
        sigma = 1
        logNormal = rng.lognormal(mu, sigma, self.Ng**2)
        self.B = np.reshape(logNormal, (self.Ng, self.Ng))

        self.contents = [Micelle(self.N, self.Ng)]
    

    def grow_micelles(self):
        """
        Decide which amphiphile is most favoured to join the micelle, and increment it..
        """
        for micelle in self.contents:

            # Weight B by the micelle composition
            weighted = micelle.composition * self.B

            # Add weighted-B to K and find maximum favoured amphiphile/s
            favouredAmphiphiles = np.sum(weighted + self.K, axis=1)
            amphiphilesToAdd = np.where(favouredAmphiphiles == favouredAmphiphiles.max())[0]

            # Add favoured amphiphile/s to micelle
            for i in amphiphilesToAdd:
                micelle.add_amphiphile(i)


    def split_micelles(self):
        """
        Check if any micelles are at or above size 2 * N, and split them.
        """
        for micelle in self.contents:
            if np.sum(micelle.composition) >= 2 * self.N:
                daughterFirst, daughterSecond = micelle.split()
                self.contents.remove(micelle)
                self.contents.append(daughterFirst)
                self.contents.append(daughterSecond)


    def run_cycle(self):
        """
        Run a single cycle of the system. This will:
        1. Check if any micelles are ready to split. If they are, replace them
        with their daughter micelles.
        2. Add the most favoured amphiphiles to all micelles.

        The split step is done last so that at the end of each cycle, the soup
        contains the same number of micelles as will receive amphiphiles in the next cycle.
        """
        self.grow_micelles()
        self.split_micelles()


    def get_composition(self):
        """
        Report the composition of all the micelles in the soup.
        """
        return np.asarray([x.composition for x in self.contents])

N = 3
Ng = 8
soup = Soup(N,Ng)

soup.B = np.zeros(shape=(Ng, Ng))
soup.B[0] = 10

predictToAdd = 0
amphiphilesToAdd = 0

for i in range(0,10):
    # How many amphiphiles do we think we will be adding this cycle?
    predictToAdd = predictToAdd + len(soup.contents)

    # How many amphiphiles will we REALLY be adding this cycle?    
    for micelle in soup.contents:
        weighted = micelle.composition * soup.B
        favouredAmphiphiles = np.sum(weighted + soup.K, axis=1)
        amphiphilesToAdd += len(np.where(favouredAmphiphiles == favouredAmphiphiles.max())[0])


    soup.run_cycle()
    print("Cycle: %s \t Micelles: %s \t Total amphiphiles at end of cycle: %s \t Predicted to add since start: %s \t Actual add since start: %s \n" % (i+1, len(soup.get_composition()), np.sum(soup.get_composition(), axis=1), predictToAdd, amphiphilesToAdd))
