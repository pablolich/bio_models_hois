include("flips.jl")
include("models.jl")

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

"""
check if sols are solutions to model
"""
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

"""
glv model with type II functional response
"""
function glvtype2(vars, pars)
    A = pars[1]
    r = pars[2]
    B = pars[3]
    BDx = B*diagm(vars)
    Aeff = A ./ (1 .+ BDx)
    Aeff[diagind(Aeff)] .= diag(A)
    return r + Aeff*vars
end 

"""
get a computationaly feasible kflips depth
"""
function getkflipsmax(nvar)
    ncombbool = [binomial(nvar, i) for i in 1:nvar] .< 1e6
    zeroinds = findall(ncombbool .== 0)
    if isempty(zeroinds)
        kflipsmax = length(ncombbool)
    else
        kflipsmax = minimum(zeroinds) - 1
    end
    return kflipsmax
end

#how deep in number of simultaneous flips to go
searchdepth = 0.25
#sample parameters
nmax = 6
#seed = parse(Int, ARGS[1])
seed = 2
rng = MersenneTwister(seed)
nsim = 2000

#initialize system by sampling parameters
#declare polynomial indeterminates as global variables
for sim in 1:nsim
    for n in 2:nmax 
        println("simulation: ", sim, " diversity: ", n)
        @var x[1:n]
        global vars = x
        #declare model parameters
        A = randn(rng, (n,n))
        r = randn(rng, (n))
        #form x0 from Aij and ri (whose signs are to be optimized)
        x0 = [collect(Iterators.flatten(A)); r] 
        #get number of variables that vary
        nvars = length(x0)
        #sample functional response coefficients (fixed)
        B = randn(rng, (n, n))
        B[diagind(B)] .= 0
        #create polynomial coefficients of reference systems (Bij are positive or negative)
        #refcoeffs = pars2coeffs(x0, B, n)
        #construct parameters of reference system 
        #pars = (refcoeffs, B, n)
        #allowed indices to flip
        #allowedinds = collect(1:length(x0))
        #determine number of kflips
        #kflipsmaxtheo = Int(ceil(n*(n+1)*searchdepth))
        #kflipsmaxcomp = getkflipsmax(nvars)
        #kflipsmax = minimum([kflipsmaxtheo, kflipsmaxcomp])
        #minimize with greedy algorithm
        #xoptgreedy = searchthroughkflips(allowedinds, kflipsmax, x0, pars)
        #form optimal parameters
        #Aopt = reshape(xoptgreedy[1:(n*n)], (n, n))
        #ropt = xoptgreedy[(n*n+1):end]
        #solve and count solutions
        println("reference system")
        result_dn_ref = solvecount(vars, (A, r, B), n, glvtype2, glvtype2poly)
        println("real system")
        result_dn = solvecount(vars, (A, r, abs.(B)), n, glvtype2, glvtype2poly)
        println("optimized system")
        #result_dn_opt = solvecount(vars, (Aopt, ropt, abs.(B)), n, glvtype2, glvtype2poly)
        #save for this matrix of coefficients
        open("../data/feas_results_ref.csv", "a") do io
            writedlm(io, result_dn_ref, ' ')
        end
        open("../data/feas_results_null.csv", "a") do io
            writedlm(io, result_dn, ' ')
        end
        #open("../data/feas_results_opt"*string(seed)*".csv", "a") do io
        #    writedlm(io, result_dn_opt, ' ')
        #end
    end
end