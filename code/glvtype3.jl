using LinearAlgebra
using HomotopyContinuation
using Random
using DelimitedFiles
using Distributions
using OptimizationOptimJL


function is_feasible(r::PathResult)
    return all(real(r.solution) .> 0)
end

function findfeasible(r::PathResult)
    #check if solution  is real
    if is_real(r)
        #check if its feasible
        if is_feasible(r)
            return true
        else
            return false
        end
    else
        return false
    end
end

"""
glv purely competitive model, sparse (polynomial version)
"""

struct Pars
    A::Array{Float64, 2}
    r::Array{Float64, 1}
    B::Array{Float64, 2}
end

"""
glv model with type II functional response
"""
function glvtype2(vars, pars)
    A = pars.A
    r = pars.r
    B = pars.B
    BDx = B*diagm(vars)
    Aeff = A ./ (1 .+ BDx)
    Aeff[diagind(Aeff)] .= diag(A)
    return r + Aeff*vars
end

"""
polynomial version of the above
"""
function glvtype2poly(vars, pars)
    A = pars.A
    r = pars.r
    B = pars.B
    #create matrix of products
    n = length(vars)
    P = Array{Expression}(undef, n, n)
    for i in 1:n
        denveci = 1 .+ B[i,:] .* vars
        prodi = prod(denveci)
        for j in 1:n
            if i==j
                P[i,i] = prodi
            else
                P[i,j] = prodi/denveci[j]
            end
        end
    end
    Aeff = A .* P
    #set diagonal to zero
    Aeff[diagind(Aeff)] .= 0.0
    polyeqs = diag(P) .* (r .+ diag(A) .* vars) + Aeff*vars
    return polyeqs
end

"""
from polynomial equations to long vector of polynomial coefficients
"""
function eqs2coeffs(polyeqs)
    n = length(polyeqs)
    coeffsvec = []
    for i in 1:n
        append!(coeffsvec,coefficients(polyeqs[i], vars))
    end
    return coeffsvec
end

"""
from vector of model parameters to vector of polynomial coefficients
"""
function pars2coeffs(x, B)
    n = size(B, 1)
    #apply new signs to model parameterws
    newA = reshape(x[1:(n*n)], (n, n))
    newr = x[(n*n+1):end]
    #compute coefficients of polynomial from new model parameters
    newparspoly = Pars(newA, newr, B)
    eqs = glvtype2poly(vars, newparspoly)
    return eqs2coeffs(eqs)
end

"""
par: parameters of the model A, B, r
x: signs of parameters Aij and ri
"""
function signdiff(x, par)
    B, ref_signs = par
    #get polynomial coefficents from model parameters
    coeffs = pars2coeffs(x, B)
    #get signs of coefficients
    signs = coeffs .> 0
    #compute number of sign differences
    return sum(ref_signs .!= signs)
end

"""
count number of positive solutions
"""
function numberfeasible(solution_vec)
    nsols = size(solution_vec, 1)
    nfeas=0
    for i in 1:nsols
        sol_i = solution_vec[i,:][1]
        if all(real(sol_i) .> 0)
            nfeas += 1
        end
    end
    return nfeas
end

function certifysolutions(model, sols, pars)
    nsol = size(sols, 1)
    cert_vec = []
    tol = 1e-6
    for i in 1:nsol
        res = model(sols[i], pars)
        if all(abs.(res) .< tol)
            push!(cert_vec, i)
        end
    end
    return cert_vec
end

function solvecount(x, pars, n, model, model_poly)
    #build system
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

"""
perform random search
"""
function randsearch(x, params, nattempts)
    nit = 0
    nx = length(x)
    ndiff = signdiff(x, params)
    ndiffnew = Inf
    while nit < nattempts
        xnew =  rand([-1;1], nx) .* abs.(x)
        #check number of different signs
        ndiffnew = signdiff(xnew, params)
        if ndiffnew < ndiff
            ndiff = ndiffnew
            x = xnew
        end
        nit += 1
    end
    return x
end

"""
perform directed search
"""
function dirsearch(x, params)
    #pick the best(s) sign flips
    #get a vector of scores
    #at each step, get rid of threads that are not best, 
    #and add those new arising that are best
    return 
end

####################################################################################
#code
####################################################################################
nmax = 6
seed = 1
rng = MersenneTwister(seed)
nsim = 2000

for sim in 1:nsim
    for n in 6:nmax
        #declare polynomial indeterminates as global variables
        @var x[1:n]
        global vars = x
        #declare model parameters
        A = randn(rng, (n,n))
        r = randn(rng, (n))
        B = randn(rng, (n, n))
        #construct reference system 
        refpars = Pars(A, r, B)
        refsyst = glvtype2poly(vars, refpars)
        #compute signs of reference system
        ref_signs = eqs2coeffs(refsyst) .> 0
        #build parameters of model to optimize
        _p = (abs.(B), ref_signs)
        #get parameters Aij and ri whose signs are to be optimized
        x0 = [collect(Iterators.flatten(A)); r] #flattening happens by stacking columns vertically
        #optimize
        xopt = randsearch(x0, _p, 20000)
        #form optimized parameters
        Aopt = reshape(xopt[1:(n*n)], (n, n))
        ropt = xopt[(n*n+1):end]
        optpars = Pars(Aopt, ropt, abs.(B))
        #solve and count solutions
        println("reference system")
        result_dn_ref = solvecount(vars, refpars, n, glvtype2, glvtype2poly)
        println("real system")
        result_dn = solvecount(vars, Pars(A, r, abs.(B)), n, glvtype2, glvtype2poly)
        println("optimized system")
        result_dn_opt = solvecount(vars, optpars, n, glvtype2, glvtype2poly)
        #save for this matrix of coefficients
        open("polynomialscreening/matching_coefficient_distributions/feas_results_ref.csv", "a") do io
            writedlm(io, result_dn_ref, ' ')
        end
        open("polynomialscreening/matching_coefficient_distributions/feas_results.csv", "a") do io
            writedlm(io, result_dn, ' ')
        end
        open("polynomialscreening/matching_coefficient_distributions/feas_results_opt.csv", "a") do io
            writedlm(io, result_dn_opt, ' ')
        end
    end
end

