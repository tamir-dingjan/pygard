from pygard import pygard
import numpy as np
import pytest

def test_soup_initialisation():
    soup = pygard.Soup(5, 10)
    micelle = soup.contents[0]
    
    assert len(micelle.composition) == soup.Ng
    assert np.sum(micelle.composition) == soup.N


@pytest.mark.parametrize(("N,Ng"), [(1,5), (5,10), (10,40), (50,100)])
def test_soup_grow(N, Ng):
    soup = pygard.Soup(N, Ng)
    
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

    