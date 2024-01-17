using LinearAlgebra
using Random
using HomotopyContinuation
using Combinatorics
using DelimitedFiles

"""
solve model_poly (filter out solutions that are not of model)
in variables x, for parameters pars, and n species
"""
function solvecount(x, pars, n, model, model_poly)
    #build system and solve it
    syst = System(model_poly(x, pars))
    res = HomotopyContinuation.solve(syst, stop_early_cb = findfeasible,
                compile = false,
                start_system = :total_degree,
                show_progress=true)
                    #certify and count all solutions
    sols = solutions(res)
    real_sols = solutions(res, only_real=true)
    cert_vec = certifysolutions(model, sols, pars)        
    cert_sols = sols[cert_vec, :]
    nsols = length(cert_sols)
    cert_vec = certifysolutions(model, real_sols, pars)        
    cert_real_sols = real_sols[cert_vec, :]
    nsolsreal = length(cert_real_sols)
    nsolsfeas = numberfeasible(cert_real_sols)
    #store results
    return [n nsols nsolsreal nsolsfeas]
end

function compute_pf(model, model_poly, nsim, nspp)

    return 
end