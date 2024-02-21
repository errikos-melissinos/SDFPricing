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
        include("zero_coupon_security.jl")
    end
end