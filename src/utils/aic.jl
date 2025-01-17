export AICc, AIC

function num_parameters(model::FBDconstant) return 3 end

function AICc(model::Model, tree::Root)
    n = 1.0
    k = num_parameters(model) 

    lnl = logL_root(model, tree)

    aic = 2*k - 2*lnl
    aicc = aic + (2*k^2 + 2*k) / (n-k-1)
    return(aicc)
end

function AIC(model::Model, tree::Root)
    n = 1.0
    k = num_parameters(model) 

    lnl = logL_root(model, tree)

    aic = 2*k - 2*lnl
    return(aic) 
end


