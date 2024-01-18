using LinearAlgebra
using Random
using HomotopyContinuation
using Combinatorics
using DelimitedFiles

"""
recast put ones in positions of vector not in allowedinds
vector: vector of signs of (possibly) smaller vector
allowedinds: components of the vector that are varied
n: original size of the vector
"""
function recast(vector, allowedinds, n)
    #get indices not in allowedinds
    allinds = collect(1:n)
    notallowedinds = allinds[(!in).(allinds, Ref(allowedinds))]
    #insert ones at fixed positions
    [insert!(vector, i, 1) for i in notallowedinds]
    return vector
end

"""
in a vector of ones, flip to -1 the inds2flip positions
n: number of model parameters (counting the fixed and the flippable ones)
inds2flip: positions of the parameters which sign is to be flipped
"""
function getvec4flip(n, inds2flip, allowedinds)
    vec = ones(Int, n)
    vec[inds2flip] .= -1
    return vec
end

"""
par: parameters of the model A, B, r
x: signs of parameters Aij and ri
"""
function signdiff(x, par)
    B, ref_coeffs = par
    #get polynomial coefficents from model parameters
    coeffs = pars2coeffs(x, B)
    #get signs of coefficients
    signs = coeffs .> 0
    ref_signs = ref_coeffs .> 0
    #compute number of sign differences
    return sum(ref_signs .!= signs)
end

"""
calculate distance between reference and current coefficients
weight signs mostly, but also numerical closeness. this will favor more similar 
coefficients in cases where the number of sign differeneces are the same.
"""
function distance(reference, current)
    signs_reference = sign.(reference)
    signs_current = sign.(current)
    return sum(signs_reference .!= signs_current) + 0.1 * norm(reference .- current)
end

"""
polynomial version of the model
"""
function glvtype2poly(vars, pars)
    A = pars[1]
    r = pars[2]
    B = pars[3]
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
x: vector of model parameters (which sign can change)
B: fixed parameters
"""
function pars2coeffs(x, B, n)
    #apply new signs to model parameterws
    newA = reshape(x[1:(n*n)], (n, n))
    newr = x[(n*n+1):end]
    #compute coefficients of polynomial from new model parameters
    newparspoly = (newA, newr, B)
    eqs = glvtype2poly(vars, newparspoly)
    return eqs2coeffs(eqs)
end

"""
from a binary vector to a vector of 1s and -1s
"""
function binary2signed(binaryvector)
    return binaryvector .+  binaryvector .- 1
end

"""
perform exhaustive search (only possible for n<5)
"""
function bruteforcesearch(x0, pars)
    #unpack parameters
    refcoeffs = pars[1]
    B = pars[2]
    n = pars[3]
    coeffs0 = pars2coeffs(x0, abs.(B), n)
    #calculate distance between reference and initial coefficients
    dist0 = distance(refcoeffs, coeffs0)
    xopt = deepcopy(x0)
    nx = length(x0)
    #iterate through all possible sign flips
    for i in 0:2^nx-1
        if rem(i,1e4) == 0
            println("Iteration: ", i)
        end
        signvec = binary2signed(digits(i, base=2, pad = nx))
        x = signvec .* abs.(x0)
        coeffs = pars2coeffs(x, abs.(B), n)
        #check number of different signs
        dist = distance(refcoeffs, coeffs)
        if dist < dist0
            dist0 = dist
            xopt = x
        end
    end
    return xopt
end

"""
greedy search
allowedinds: vector of allowed indices to flip
kflips: number of flips allowed at a time
n: number of species
"""
function greedysearch(allowedinds, kflips, x0, pars)
    #unpack parameters
    refcoeffs = pars[1]
    B = pars[2]
    n = pars[3]
    coeffs0 = pars2coeffs(x0, abs.(B), n)
    nvars = length(allowedinds)
    nvarori = length(x0)
    xopt = deepcopy(x0)
    #calculate distance between reference and initial coefficients
    dist0 = distance(refcoeffs, coeffs0)
    comb = 0
    notflipped = collect(1:n)
    #get all possible combinations of the allowed indices in groups of kflips
    mat2flip = collect(combinations(allowedinds, kflips))
    ncombs = length(mat2flip)
    for i in 1:ncombs
        if rem(i,1e4) == 0 && i > 0
            println("Iteration: ", i)
        end
        #turn to flipping vector
        flipvec = getvec4flip(nvarori, mat2flip[i], allowedinds)
        #flip corresponding indices of x0
        x = x0 .* flipvec
        #compute polynomial coefficients of flipped model with positive B
        newcoeffs = pars2coeffs(x, abs.(B), n)
        #get distance to reference model
        dist = distance(refcoeffs, newcoeffs)
        if dist < dist0
            #set new minimum distance
            dist0 = dist
            #set new optimal x
            xopt = x
            comb = i
            #set new vector of indices to flip (only the ones that have not been flipped yet)
            allowedinds = findall(flipvec .== 1)
        end
    end
    #after traversing through all kflips, there are 2 options; 
    if length(allowedinds) !== nvars
        #something has been flipped, call the function again, changing the allowedinds vector
        println("Combination number ", comb, " was succesful")
        println("fix advantageous flips, move down in recursion")
        return greedysearch(allowedinds, kflips, xopt, pars)
    else
        #2: nothing has been flipped, break the search
        println("no better answer was found")
        return (x0, dist0)
    end
end

"""
traverse through kflips
"""
function searchthroughkflips(allowedinds, kflipsmax, x0, pars)
    #number of variables
    nvar = length(x0)
    #initialize storing of best solution, and score per kflips
    xoptmat = zeros((kflipsmax, nvar))
    distvec = zeros(kflipsmax)
    for kflips in 1:kflipsmax
        println("trying with ", kflips, " (out of ", kflipsmax, ")", " at a time")
        #minimize
        xopt, dist = greedysearch(allowedinds, kflips, x0, pars)
        #store
        xoptmat[kflips,:] = xopt
        distvec[kflips] = dist
    end
    bestsolind = argmin(distvec)
    return xoptmat[bestsolind,:]
end

"""
main code
"""
function main()
    nmax = 4
    seed = 2
    rng = MersenneTwister(seed)
    nsim = 2000

    #initialize system by sampling parameters
    #declare polynomial indeterminates as global variables
    for n in 2:nmax
        for sim in 1:nsim
            @var x[1:n]
            global vars = x
            #declare model parameters
            A = randn(rng, (n,n))
            r = randn(rng, (n))
            #get parameters Aij and ri whose signs are to be optimized
            x0 = [collect(Iterators.flatten(A)); r] 
            nvars = length(x0)
            B = randn(rng, (n, n))
            refcoeffs = pars2coeffs(x0, B, n)
            #construct reference system 
            pars = (refcoeffs, B, n)
            #minimize with bruteforce
            xoptbrute = bruteforcesearch(x0, pars)
            xoptgreedy = ones(nvars)
            distvec = zeros(nvars)
            #allowed indices to flip
            allowedinds = collect(1:length(x0))
            #minimizie with greedy algorithm
            bestfound = false
            kflips = 1
            same = false
            while !bestfound
                xoptgreedy, dist = greedysearch(allowedinds, kflips, x0, pars)
                #check if answer is the same
                same = xoptgreedy == xoptbrute
                if same
                    println("found best answer for ", kflips, "flips")
                    bestfound = true
                end
                #save for this matrix of coefficients
                resultdn = [n sim kflips dist bestfound]
                open("comparedists_greedy_brute.csv", "a") do io
                    writedlm(io, resultdn, ' ')
                end
                kflips += 1
            end
        end
    end
end