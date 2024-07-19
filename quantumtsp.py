DISTANCE_MATRIX = [
    [0, 1, 4, 6],
    [1, 0, 2, 5],
    [4, 2, 0, 3],
    [6, 5, 3, 0]
]

n = len(DISTANCE_MATRIX)
print("n =", n)

from dwave.optimization.generators import traveling_salesperson
model = traveling_salesperson(distance_matrix=DISTANCE_MATRIX)

from dwave.system import LeapHybridNLSampler
sampler = LeapHybridNLSampler()

results = sampler.sample(model, label='TSP')

route, = model.iter_decisions()
print("TOUR:")
print(route.state(0))
