from pygard import pygard
import numpy as np

def test_micelle_composition_initialisation():
    micelle = pygard.Micelle(5,10)

    assert len(micelle.composition) == 10


def test_micelle_split():
    micelle = pygard.Micelle(5,10)

    daughterFirst, daughterSecond = micelle.split()
    
    assert np.sum(np.asarray([daughterFirst.composition, daughterSecond.composition])) == np.sum(micelle.composition)


def test_micelle_add():
    micelle = pygard.Micelle(5,10)
    compositionPre = micelle.composition[0]
    
    micelle.add_amphiphile(0)
    compositionPost = micelle.composition[0]
    
    assert compositionPre == compositionPost - 1