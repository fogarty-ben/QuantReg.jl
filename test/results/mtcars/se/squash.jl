using DelimitedFiles, JLD

squashed = Dict()

for τ in (0.25:0.25:0.75)
    squashed[τ] = Dict()
    for invers in (true, false)
        squashed[τ][invers] = Dict()
        for α in (0.01, 0.05, 0.10)
            squashed[τ][invers][α] = Dict()
            for hs in (true, false)
                squashed[τ][invers][α][hs] = Dict()
                for iid in (true, false)
                    squashed[τ][invers][α][hs][iid] = Dict()
                    for interp in (true, false)
                        names = ["tau", τ, "invers", invers, "alpha", α, "hs", hs, "iid",
                                 iid, "interp", interp]
                        fn = fp = joinpath(@__DIR__, join(names, "_") * ".txt")
                        if isfile(joinpath(@__DIR__, fn))
                            squashed[τ][invers][α][hs][iid][interp] = readdlm(fp, '\t',
                                                                              Float64, '\n')
                        else
                            squashed[τ][invers][α][hs][iid][interp] = "ERROR"
                        end
                    end
                end
            end
        end
    end
end

save(joinpath(@__DIR__, "../inf.jld"), "inference", squashed)