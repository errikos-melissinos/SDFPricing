using Test
import SDFPricing as sdf

#* test the pickIndex function
#* the first ```paths``` number of paths should be assigned to the first initial value, then to the second and so on 
paths = 1000
@test sdf.pickIndex(1, paths) == 1
@test sdf.pickIndex(paths, paths) == 1
@test sdf.pickIndex(paths + 1, paths) == 2
@test sdf.pickIndex(2paths, paths) == 2
@test sdf.pickIndex(2paths + 1, paths) == 3
