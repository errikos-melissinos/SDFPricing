import Pkg
Pkg.test("ConsumptionFinance"; coverage=true)

using Base.Filesystem
for file in readdir("ConsumptionFinance/src")
    if endswith(file, ".cov")
        mv(joinpath("ConsumptionFinance/src", file), joinpath("ConsumptionFinance/coverage", file))
    end
end
for file in readdir("ConsumptionFinance/test")
    if endswith(file, ".cov")
        mv(joinpath("ConsumptionFinance/test", file), joinpath("ConsumptionFinance/coverage", file))
    end
end


using Coverage
covdir = "ConsumptionFinance/coverage"
lcovfile = "ConsumptionFinance/coverage"
cov = process_folder(covdir)
LCOV.writefile("ConsumptionFinance/coverage/coverage-lcov.info", cov)