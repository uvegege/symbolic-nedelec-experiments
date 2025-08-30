struct Triangle end
dim(::Triangle) = 2

struct LegendreP end
struct DubinerP end 

function N1k(k, element = Triangle())

    space = Rspace(k, dim(element))
    eqs = Num[]
    coeff = space.coefs
    R = space.exprs
    edge_moments_NI!(eqs, space, element, LegendreP())
    face_moments_NI!(eqs, space, element, DubinerP())
    Nj = map(eachindex(eqs)) do idx
        free_vals = ntuple(i -> i == idx ? 1 : 0, length(eqs))
        soli = symbolic_linear_solve(eqs .~ free_vals, coeff)
        Symbolics.simplify.(Symbolics.substitute(R, Dict(coeff .=> soli)))
    end
    
    return Nj
end

function sym_equation_to_matrix(equations, coeff)
    A = zeros(Float64, length(coeff), length(coeff))
    for c in eachindex(coeff) 
        for e in eachindex(equations)
            A[e, c] = Symbolics.coeff(equations[e], coeff[c])
        end
    end
    return A
end

function N1k2(k, element = Triangle())

    space = Rspace(k, dim(element))
    eqs = Num[]
    coeff = space.coefs
    R = space.exprs

    edge_moments_NI!(eqs, space, element, LegendreP())
    face_moments_NI!(eqs, space, element, DubinerP())

    A = sym_equation_to_matrix(eqs, coeff)

    Nj = map(eachindex(coeff)) do idx
        fv = SVector{length(eqs)}([i == idx ? 1 : 0 for i in eachindex(eqs)])
        soli = A\fv
        soli = roundcoefficients.(simplify_fractions(soli))
        Symbolics.simplify.(Symbolics.substitute(R, Dict(coeff .=> soli)))
    end
    
    return Nj
end
