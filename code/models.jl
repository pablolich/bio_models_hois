"""
glv type 2 dense
"""
function glvtype2dense(vars, pars)
    A = pars[1]
    r = pars[2]
    B = pars[3]
    S = diagm(1 ./ (1 .+ B*vars))
    Aeff = A*S
    Aeff[diagind(Aeff)] .= diag(A)
    return r - Aeff*vars
end

"""
glv type2 dense (test)
"""
function testmodel3spp(x, A, B)
  return 1 - A[1,1]*x[1] - A[1,2]*x[2]/(1+B[2,1]*x[1] + B[2,3]*x[3]) - 
           A[1,3]*x[3]/(1+B[3,2]*x[2] + B[3,1]*x[1])
end

"""
glv type 2 dense (polynomial version)
"""
function glvtype2densepoly(vars, pars)
    A = pars[1]
    r = pars[2]
    B = pars[3]
    #form rime rescaling products
    denvec = 1 .+ B*vars
    n = length(vars)
    prodtot = prod(denvec)
    #create matrix of products
    P = Array{Expression}(undef, n, n)
    for i in 1:n
        for j in 1:n
            if i==j
                P[i,i] = prodtot/denvec[i]
            else
                P[i,j] = prodtot/(denvec[i]*denvec[j])
            end
        end
    end
    Aeff = A .* P
    #set diagonal to zero
    Aeff[diagind(Aeff)] .= 0.0
    polyeqs = diag(P) .* (r .- diag(A) .* vars) - Aeff*vars
    return polyeqs
end
