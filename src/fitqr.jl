
# Move to other 
function fit(model::QuantRegModel)
    if !model.fitted
        fitfxn = nothing
        if model.method == "br"
            fitfxn = fitbr
        elseif model.method == "gurobi"
            fitfxn = fitgurobi
        end
    end
end

# Fit using Barrodale-Roberts from R quantreg
function fitbr
end

# Fit using Gurobi
function fitgurobi
end