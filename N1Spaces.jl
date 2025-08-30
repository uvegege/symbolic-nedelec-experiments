function Pspace(x::Symbolics.Arr{Num, 1}, degree)
    monos = generar_monomios(x, degree)
    monos = ordenar_grlex(monos, x)
    @variables cij[1:length(x), 1:length(monos)]
    Cij = Symbolics.scalarize(cij)
    pspace = [(dot(@views(Cij[i, :]), monos)) for i in 1:length(x)]
    return (; monos = monos, coefs = cij, space = pspace)
end


"""
https://www.dealii.org/reports/nedelec/nedelec.pdf

The dimension space is k in the case d = 2 and k(k+2) for d = 3

p_{ij} = sum_{j=1}^3 a_{ij} x_j, i = 1, 2, 3

The condition for p being in S^k is
    
## Example

@variables x[1:3]
ssp = Sspace(x, 1)

julia> ssp.space
3-element Vector{Tuple{Num, Num, Num}}:
 (0, x[3], -x[2])
 (x[2], -x[1], 0)
 (x[3], 0, -x[1])

@variables x[1:2]
ssp = Sspace(x, 3)

julia> ssp.space
3-element Vector{Tuple{Num, Num}}:
 (x[2]^3, -x[1]*(x[2]^2))
 (-(x[1]^2)*x[2], x[1]^3)
 (x[1]*(x[2]^2), -(x[1]^2)*x[2])


"""
function Sspace(x::Symbolics.Arr{Num, 1}, degree)
    # garantizo que p * x es de grado k + 1 de manera exacta y no incluye todos los grados
    # De esta manera el sistema no está sobredeterminado
    mp = generar_monomios_exacto(x, degree) 
    mp = ordenar_grlex(mp, x)
    dims = length(x)
    jmax = length(mp)
    @variables aij[1:dims, 1:jmax]
    Aij = Symbolics.scalarize(aij)

    monos_j = generar_monomios_exacto(x, degree+1)
    monos_j = ordenar_grlex(monos_j, x)
    
    p_i = [sum(Aij[ii,jj]*mp[jj] for jj in 1:jmax) for ii in 1:dims]
    eqs = dot(p_i, x)
    eqs = Symbolics.expand(eqs)
    
    relations = map(mj -> Symbolics.coeff(eqs, mj), monos_j)
    substitutions = []
    for r in relations
        ri = Symbolics.get_variables(r)
        if length(ri) > 1
            push!(substitutions, ri[1] => -sum(ri[2:end]))
        else
            push!(substitutions, ri[1] => 0)
        end
    end
    sdict = Dict(substitutions)
    aijvars = Symbolics.substitute(Aij, sdict)
    aijvars = vec(aijvars)
    aijvars = Symbolics.coeff.(aijvars)
    free_vars = filter(x-> !isempty(x), Symbolics.get_variables.(aijvars)) |> x->map(i->i[1], x)
    nfreedom = length(free_vars)
    spaces = Vector{NTuple{dims, Num}}(undef, nfreedom)
    for idx in 1:nfreedom
        free_vals = ntuple(i -> i == idx ? 1 : 0, nfreedom)
        p_x = copy(p_i)
        p_x = Symbolics.substitute.(p_x, Ref(sdict))
        dv = Dict(free_vars .=> free_vals)
        p_x = Symbolics.substitute(p_x, dv)
        spaces[idx] = tuple(p_x...)
    end
    @variables scoefs[1:nfreedom]
    return (;monos = x, coefs = scoefs, space = spaces)
end

function Rspace(degree, dimension)
    @variables u[1:dimension]
    pspace = Pspace(u, degree)
    sspace = Sspace(u, degree+1)

    V = copy(pspace.space)
    for (idx, sk) in enumerate(sspace.space)
        for idx2 in eachindex(V)
            V[idx2] += sspace.coefs[idx] * sk[idx2]
        end
    end
    
    allcoefs = vcat(vec(pspace.coefs), vec(sspace.coefs))
    return (; x = u, exprs = V, coefs = allcoefs, degree = degree)
end