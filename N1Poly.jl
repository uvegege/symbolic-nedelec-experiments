using DynamicPolynomials
using MultivariatePolynomials
# https://discourse.julialang.org/t/typeerror-when-using-symbolics-jl-with-polynomials-jl/57842/5 
# There is a problem with show() when the coefficients of the polynomial are `Num` (Symbolics.jl)

generar_monomios(x::Tuple, k::Int64) = monomials(x, 0:k)
generar_monomios(x::Vector, k::Int64) = monomials(x, 0:k)
generar_monomios_exacto(x::Tuple, k::Int64) = monomials(x, k)
generar_monomios_exacto(x::Vector, k::Int64) = monomials(x, k)

# Polynomial to Symbolics.jl
function polyeqtosym(eq)
    a = eq.a
    b = eq.x
    monos = b.vars
    exponents = b.Z
    n = length(monos)
    @variables u[1:n]
    result = Num(0)
    for (a_i, z) in zip(a, exponents)
        result += a_i * prod(i -> i[1]^i[2], zip(u,z))
    end
    return result
end


function Pspace(x, k)
    monos = generar_monomios(x, k)
    @variables cij[1:length(x), 1:length(monos)]
    Cij = Symbolics.scalarize(cij)
    pspace = [(dot(@views(Cij[i, :]), monos)) for i in 1:length(x)]
    return (; monos = monos, coefs = cij, space = pspace)
end


function Sspace(x, k)
    # Garantizar que p * x es de grado k + 1 de manera exacta y no incluye todos los grados
    # De esta manera el sistema no está sobredeterminado
    mp = generar_monomios_exacto(x, k) 
    dims = length(x)
    jmax = length(mp)
    @variables aij[1:dims, 1:jmax]
    Aij = Symbolics.scalarize(aij)

    #monos_j = generar_monomios_exacto(x, k+1)
    p_i = [sum(Aij[ii,jj]*mp[jj] for jj in 1:jmax) for ii in 1:dims]
    eqs = dot(p_i, x)

    relations = eqs.a
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
    for i in eachindex(p_i)
        p_i[i] = dot(Symbolics.substitute.(p_i[i].a, Ref(sdict)), p_i[i].x)
    end
    #spaces = Vector{NTuple{dims, Num}}(undef, nfreedom)
    spaces = Vector{Any}(undef, nfreedom)
    for idx in 1:nfreedom
        free_vals = ntuple(i -> i == idx ? 1 : 0, nfreedom)
        p_x = copy(p_i)
        for i in eachindex(p_x)
            p_x[i] =  dot(Symbolics.substitute.(p_x[i].a, Ref(Dict(free_vars .=> free_vals))),  p_i[i].x)
        end
        spaces[idx] = tuple(p_x...)
    end
    @variables scoefs[1:nfreedom]
    return (;monos = x, coefs = scoefs, space = spaces)
end


function Rspace_poly(k, dimension)
    if dimension == 2
        @polyvar x[1:2]
    else
        @polyvar x[1:3]
    end
    pspace = Pspace(x, k)
    sspace = Sspace(x, k+1)

    V = copy(pspace.space)
    for (idx, sk) in enumerate(sspace.space)
        for idx2 in eachindex(V)
            V[idx2] += sspace.coefs[idx] * sk[idx2]
        end
    end
    
    allcoefs = vcat(vec(pspace.coefs), vec(sspace.coefs))
    return (; x = x, exprs = V, coefs = allcoefs, degree = k)
end


# Edge Tangents
edge_tangents_poly(::Triangle) = (SA[1.0, 0.0], SA[0.0, 1.0], SA[-1.0, 1.0])
function edge_parametrization_poly(::Triangle, vars, s) 
    return ((vars[1], vars[2]) => (polynomial(s), 0),
     (vars[1], vars[2]) => (0, polynomial(s)), 
     (vars[1], vars[2]) => (polynomial(s), polynomial(1-s)))
end

face_parametrization_poly(vars, ::Triangle) = 1 - vars[1]

function symbolic_edge_moment_poly(v, t, L, pe, s)
    vl = dot(v, t)
    vl = MultivariatePolynomials.subs(vl, pe)
    intermedio = antidifferentiate(vl * L, s)
    return Symbolics.simplify(intermedio(s => 1) - intermedio(s => 0))
end

function edge_moments_NI_poly!(eqs, space, element = Triangle(), ptype = LegendreP())
    if dim(element) == 2
        @polyvar s
        s₀ = 2*s - 1 # Map to (-1, 1)
        V = space.exprs
        vars = space.x
        t = edge_tangents_poly(element)
        p = edge_parametrization_poly(element, vars, s)
        L = [legendrep(i, s₀) for i in 0:space.degree]
        for (pe, te) in zip(p, t)
            for le in L
                push!(eqs, symbolic_edge_moment_poly(V, te, le, pe, s))
            end
        end
        return nothing
    end
    @error "Not Implemented for $element and $ptype"
end

function face_moments_NI_poly!(eqs, space, element = Triangle(), ptype = DubinerP())
    if dim(element) == 2
        @polyvar s
        @polyvar t
        vars = space.x
        vl = space.exprs
        param = face_parametrization_poly((s,t), element)
        φ = MultivariatePolynomials.subs(vl, (vars[1] , vars[2]) => (t, s))
        vars = (s, t)
        M = generar_monomios(vars, space.degree-1)
        K = dubiner_coefs(space.degree-1)
        mons = map(k -> MultivariatePolynomials.polynomial(k,M), K)
        for m in mons
            push!(eqs, symbolic_face_moments_poly!(φ, (0,m), param, (s, t))) 
            push!(eqs, symbolic_face_moments_poly!(φ, (m,0), param, (s, t)))
        end
        return nothing
    end
    @error "Not Implemented for $element and $ptype"
end


function symbolic_face_moments_poly!(φ, L, pe, vars)
    s, t = vars
    intermedio = antidifferentiate(dot(φ, L), t)
    int_t = subs(intermedio, t => pe) - subs(intermedio, t => 0)
    int_s = antidifferentiate(int_t, s)
    return Symbolics.simplify(int_s(s => 1) - int_s(s => 0))
end

function N1kPoly(k, element = Triangle())

    space = Rspace_poly(k, dim(element));
    eqs = Num[];
    coeff = space.coefs;
    R = space.exprs;

    edge_moments_NI_poly!(eqs, space, element, LegendreP())
    face_moments_NI_poly!(eqs, space, element, DubinerP())
    
    A = sym_equation_to_matrix(eqs, coeff)

    Nj = map(eachindex(eqs)) do idx
        fv = SVector{length(eqs)}([i == idx ? 1 : 0 for i in eachindex(eqs)])
        soli = A\fv
        soli = roundcoefficients.(simplify_fractions(soli))
        dv = Dict(coeff .=> soli)
        map(eachindex(R)) do i
            Symbolics.simplify.(Symbolics.substitute(polyeqtosym(R[i]), dv))
        end
    end
    
    return Nj
end


