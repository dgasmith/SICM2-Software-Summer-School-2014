"""
Conway's game of life vectorized numpy example.

Probably not the easiest example.
"""
import numpy as np
import scipy.spatial.distance as sd
from scipy.signal import convolve2d

def conway(population, generations=100):
    """Runs Conway's game of life on an initial population."""

    for i in range(generations):
        population = evolve(population)

    return population


def evolve(population_grid):
    """Use convolve to update conways game"""

    cv_array = np.ones((3,3))
    nbrs_count = convolve2d(population_grid, cv_array, mode='same', boundary='wrap') - population_grid
    mask = (nbrs_count == 2)
    mask &= population_grid.astype(np.bool)
    mask |= (nbrs_count == 3)    

    return mask


blinker = [np.array([0,0,0]), np.array([-1,0,1])]
blinker[0] += np.random.randint(1,9,1) 
blinker[1] += np.random.randint(1,9,1) 

grid = np.zeros((10,10))
grid[tuple(blinker)] = 1

end = conway(grid, generations=101)

print 'Start: ', blinker
print 'End: ', end
