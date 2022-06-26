from pygard import pygard
import numpy as np
import pytest

def test_soup_initialisation():
    soup = pygard.Soup(5, 10)
    micelle = soup.contents[0]
    
    assert len(micelle.composition) == soup.Ng
    assert np.sum(micelle.composition) == soup.N


def test_soup_grow():
    soup = pygard.Soup(5, 10)
    
    # Determine how many maximally favoured amphiphiles will be added
    favs = np.sum((soup.B * soup.contents[0].composition) + soup.K, axis=0)
    numToAdd = len(np.where(favs == favs.max())[0])

    soup.grow_micelles()

    assert np.sum(soup.contents[0].composition) == soup.N + numToAdd


def test_soup_split():
    soup = pygard.Soup(5, 10)

    soup.contents[0].composition[0] += 5

    soup.split_micelles()
    assert len(soup.contents) == 2
    

def test_soup_composition():
    soup = pygard.Soup(5,10)
    x = soup.get_composition()
    assert np.sum(x) == 5

@pytest.mark.parametrize("num_cycles", (x for x in range(1,10)))
def test_soup_cycle(num_cycles):
    N = 5
    Ng = 10
    soup = pygard.Soup(N, Ng)

    # Use a highly biased B matrix to ensure only a single lipid is added in each cycle
    soup.B = np.zeros(shape=(Ng, Ng))
    soup.B[0] = 10
    
    toAdd = 0
    # Since the micelle splitting is random, use the micelle populations
    # to infer how many lipids were added
    for i in range(1,num_cycles+1):
        toAdd = toAdd + len(soup.contents) # This is currently wrong
        soup.run_cycle()

    assert np.sum(soup.get_composition()) == N + toAdd


    