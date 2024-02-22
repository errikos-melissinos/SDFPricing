using Test
import ConsumptionFinance as cf

#* test the pickIndex function
#* the first ```paths``` number of paths should be assigned to the first initial value, then to the second and so on 
paths = 1000
@test cf.pickIndex(1, paths) == 1
@test cf.pickIndex(paths, paths) == 1
@test cf.pickIndex(paths + 1, paths) == 2
@test cf.pickIndex(2paths, paths) == 2
@test cf.pickIndex(2paths + 1, paths) == 3




mod = cf.Model()
(price,) = cf.solveModel(mod)

price(1.0, 0.01)