using Test
using SafeTestsets

@time begin
    @time @safetestset "Simple Tests" begin
        include("simple_tests.jl")
    end

    @time @safetestset "Test defineMySolve" begin
        include("define_my_solve.jl")
    end

    @time @safetestset "Test zeroCouponSecurity" begin
        include("solve_tests.jl")
    end
end