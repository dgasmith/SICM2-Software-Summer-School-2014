"""
Conway's game of life vectorized numpy example.

Probably not the easiest example.
"""
import numpy as np
import scipy.spatial.distance as sd

def conway(population, generations=100):
    """Runs Conway's game of life on an initial population."""

    for i in range(generations):
        population = evolve(population)

    return population


def evolve(population):
    """Evolves the population by one generation."""

    # Get a unique set of discrete cells that need to be checked
    cbdist = sd.cdist(population, population, metric='cityblock')
    print (cbdist==1).sum(axis=0)

    exit()

    active_cells = population[:]
    #ipdb.set_trace()
    for cell in population:
        for neighbor in neighbors(cell):
            if neighbor not in active_cells:
                active_cells.append(neighbor)

    # For each cell in the set, test if it lives or dies
    new_population = []
    for cell in active_cells:
        count = sum(neighbor in population for neighbor in neighbors(cell))
        if count == 3 or (count == 2 and cell in population):
            if cell not in new_population:
                new_population.append(cell)

    # Return the new surviving population
    return new_population

def neighbors(popultion):
    """Generate 2D array of all possible neighbors"""
    shift = np.array([[1, 0], [-1, 0], [0, 1], [0, -1],
                     [1, 1], [1, -1], [-1, 1], [-1, -1]])

    return np.array(cell)+shift


#glider = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 2)]
#glider = np.array(glider)
#end = conway(glider)

blinker = np.array([(0, 0), (1, 0), (-1, 0)])
end = conway(blinker)

print 'Start: ', blinker
print 'End: ', end
