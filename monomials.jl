"""
Generate all the monomials of degree <= k
"""
function generar_monomios(x, k)
    dims = length(x)
    nmonos = binomial(dims+k,dims)
    monos = Vector{Num}(undef, nmonos)
    xi = Symbolics.scalarize(x)
    idx = 1
    for ijk in Iterators.product(ntuple(_ -> range(0, k), dims)...)
        sum(ijk) > k && continue
        m = 1
        for (xj, ex) in zip(xi, ijk)
            m *= xj^ex
        end
        monos[idx] = m
        idx += 1
        idx > nmonos && break
    end
    return ordenar_grlex(monos, x)
end

function generar_monomios_exacto(x, k)
    dims = length(x)
    nmonos = binomial(dims+k-1,k)
    monos = Vector{Num}(undef, nmonos)
    xi = Symbolics.scalarize(x)
    idx = 1
    for ijk in Iterators.product(ntuple(_ -> range(0, k), dims)...)
        if sum(ijk) === k
            m = 1
            for (xj, ex) in zip(xi, ijk)
                m *= xj^ex
            end
            monos[idx] = m
            idx += 1
            idx > nmonos && break
        end

    end
    return ordenar_grlex(monos, x)
end



# Dubiner basis

expsimp(x) = Symbolics.simplify(x; expand = true)
#expsimp(x) = Symbolics.expand(Symbolics.simplify(x))

function dubiner(s,t,m,n; rational = false)
    if rational == true
        return simplify(2^m * expsimp(jacobip(m, 0, 0, 2*s//(1-t)-1)) * expsimp((1-t)^m) * expsimp(jacobip(n, 2//1*m+1, 0, 2*t-1)); expand = true)
    else
        return simplify(2^m * expsimp(jacobip(m, 0, 0, 2*s/(1-t)-1)) * expsimp((1-t)^m) * expsimp(jacobip(n, 2*m+1, 0, 2*t-1)); expand = true)
    end
end

function dubiner_basis(vars, degree; rational = true)
    s, t = vars
    basis = []
    for n in 0:degree
        for m in 0:(degree - n)
            ϕ = dubiner(s, t, n, m; rational = rational) 
            ϕ = Symbolics.expand(ϕ)
            ϕ = Symbolics.simplify_fractions(ϕ)
            ϕ = Symbolics.simplify(ϕ)
            push!(basis, ϕ)
        end
    end
    return basis
end

function get_exponents(expr, vars)
    [Symbolics.degree(expr, v) for v in vars]
end

function grlex_less_than(e1, e2)
    deg1 = sum(e1)
    deg2 = sum(e2)
    if deg1 != deg2
        return deg1 < deg2
    else
        return e1 < e2  
    end
end

function ordenar_grlex(monos, vars)
    exp_pairs = [(get_exponents(m, vars), m) for m in monos]

    for i in 1:length(exp_pairs)-1
        for j in i+1:length(exp_pairs)
            if grlex_less_than(exp_pairs[j][1], exp_pairs[i][1])
                exp_pairs[i], exp_pairs[j] = exp_pairs[j], exp_pairs[i]
            end
        end
    end

    return [pair[2] for pair in exp_pairs]
end


function dubiner_coefs(k)
    @variables o[1:2]
    monos = generar_monomios(o, k)
    monos = ordenar_grlex(monos, o)

    dubbasis = dubiner_basis(o, k)

    return map(dubbasis) do db
        getindex.(Symbolics.factors.(Symbolics.terms(db)),1)
    end

end
