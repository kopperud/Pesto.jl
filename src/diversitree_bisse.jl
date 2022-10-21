export anc_state_prob_bisse

function anc_state_prob_bisse(treefile, datafile, model::SSEconstant)
    RCall.@rput treefile
    RCall.@rput datafile

    eta = model.η
    lambda = model.λ
    mu = model.μ

    RCall.@rput eta
    RCall.@rput lambda
    RCall.@rput mu

    RCall.R"""
    library(diversitree)
    library(ape)

    pars <- c(lambda, mu, eta, eta)
    phy <- read.nexus(treefile)
    tmp <- read.csv(datafile,header=TRUE, row.names=1)
    states <- tmp[,"state"]
    names(states) <- rownames(tmp)

    # calculate likelihood
    lik <- make.bisse(phy, states, control = list("backend", "deSolve"), strict=FALSE)
    rate <- pars[1]
    num_taxa <- length(states)
    #lnl <- lik(pars,root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE, intermediates=TRUE)
    lnl <- lik(pars,root=ROOT.FLAT, condition.surv=TRUE, intermediates=TRUE)

    D <- attr(lnl, "vals")[3:4]
    E <- attr(lnl, "vals")[1:2]
    lq <- attr(lnl, "intermediates")$lq

    freqs <- c(0.5, 0.5)
    nonextinct <- (1 - E)^2
    prob <- sum(freqs * D / (nonextinct * lambda)) 

    logL <- log(mean(D / (lambda * (1 - E)^2))) + sum(lq)
#    cat("logL =", logL, "\n")
    print(paste0("logL =", logL))

    # now infer marginal ancestral state reconstructions
    anc_states = asr.marginal(lik, pars)

    anc_states <- t(anc_states) 
    """
    RCall.@rget anc_states
    return(anc_states)
end
