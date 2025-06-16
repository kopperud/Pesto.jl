function check_if_converged(optres::Optim.MultivariateOptimizationResults)
    stp = optres.stopped_by  

    converged = stp.x_converged || stp.f_converged || stp.g_converged
end

