from conway import *

def test_neighbors_at_origin():
    result = [(1,1), (-1,-1), (0,1), (1,0), (-1,1), (1,-1), (-1,0), (0,-1)]
    nb = neighbors((0,0))
    assert( set(result) == set(nb) )    

def test_neighbors_at_negative_quadrant():
    result = [(0, -1), (-2, -1), (-1, 0), (-1, -2), (0, 0), (0, -2), (-2, 0), (-2, -2)]
    nb = neighbors((-1,-1))
    assert( set(result) == set(nb) ) 

def test_blinker():
    blinker = [(-1,0), (0,0), (1,0)]
    result = conway(blinker, generations=2) 

    assert( set(result) == set(blinker) ) 

def test_block():
    block = [(0,0), (0,1), (1,0), (1,1)]
    result = conway(block, generations=2) 

    assert( set(result) == set(block) ) 
