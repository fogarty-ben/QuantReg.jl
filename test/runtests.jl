using Test, CSV, QuantReg, DelimitedFiles, JLD

const resultspath = joinpath(@__DIR__, "results/")

mtcars = CSV.read(joinpath(@__DIR__, "../data/mtcars.csv"))

mtcarscoef = readdlm(joinpath(resultspath, "./mtcars/coef.txt"), '\t', Float64, '\n')
mtcarsresid = readdlm(joinpath(resultspath, "./mtcars/resid.txt"), '\t', Float64, '\n')
mtcarsdual = readdlm(joinpath(resultspath, "./mtcars/dual.txt"), '\t', Float64, '\n')
mtcarsyhat = readdlm(joinpath(resultspath, "./mtcars/yhat.txt"), '\t', Float64, '\n')
mtcarsinf = load(joinpath(resultspath, "./mtcars/inf.jld"))["inference"]
@testset "mtcars" begin
    @testset "fit method=$fitmethod, τ=$tau" for fitmethod in ("br", "fn"),
                                                 tau in (0.1:0.1:0.9)
            mtcarsmodel = rq(@formula(mpg ~ disp + hp + cyl), mtcars; fitmethod=fitmethod, τ=tau)
            idx = Int(tau * 10)
            @test all(isapprox.(mtcarsmodel.fit.coef, mtcarscoef[:, idx]; atol=1e-9))
            @test all(isapprox.(mtcarsmodel.fit.resid, mtcarsresid[:, idx]; atol=1e-9))
            if (fitmethod == "br") & (tau == 0.9)
                @test all(isapprox.(mtcarsmodel.fit.dual, mtcarsdual; atol=1e-9))
            end
            @test all(isapprox.(mtcarsmodel.fit.yhat, mtcarsyhat[:, idx]; atol=1e-9))
    end
    @testset "inf τ=$τ, invers=$invers, α=$α, hs=$hs, iid=$iid, interp=$interp" for
        invers in (true, false), α in (0.01, 0.05, 0.10), hs in (true, false),
        iid in (true, false), interp in (true, false), τ in (0.25:0.25:0.75)
        local inf
        try
            mtcarsmodel = rq(@formula(mpg ~ disp + hp + cyl), mtcars; invers=invers, α=α,
                             hs=hs, iid=iid, interpolate=interp, τ=τ)
            if invers
                inf = hcat(mtcarsmodel.inf.lowerci, mtcarsmodel.inf.upperci)
            else
                inf = mtcarsmodel.inf.σ  
            end
        catch e
            inf = "ERROR"
        finally
            if inf == "ERROR"
                @test inf == mtcarsinf[τ][invers][α][hs][iid][interp]
            else
                # println(inf)
                # println(mtcarsinf[string(τ)][string(invers)][string(α)][string(hs)][string(iid)][string(interp)])
                @test all(isapprox.(inf, mtcarsinf[τ][invers][α][hs][iid][interp]; atol=1e-5))
            end
        end
    end
end

