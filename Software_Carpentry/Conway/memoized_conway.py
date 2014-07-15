"""
Conway's game of life example, part three.

This code has a bug that makes the game return improper results.
"""

class memoized(object):
   def __init__(self, func):
      self.func = func
      self.cache = {}
      self.count = 0
   def __call__(self, *args):
      if args in self.cache:
         self.count +=1
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value

def conway(population, generations=100):
    """Runs Conway's game of life on an initial population."""

    for i in range(generations):
        population = evolve(population)

    return population


def evolve(population):
    """Evolves the population by one generation."""

    # Get a unique set of discrete cells that need to be checked
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

@memoized
def neighbors(cell):
    x, y = cell
    return [(x + 1, y), (x - 1, y), 
            (x, y + 1), (x, y - 1),
            (x + 1, y + 1), (x + 1, y - 1), 
            (x - 1, y + 1), (x - 1, y - 1)]


glider = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 2)]
blinker = [(-1,0), (0,0), (1,0)]

start = glider
#start = blinker
end = conway(start)

print 'Start: ', start
print 'End: ', end
