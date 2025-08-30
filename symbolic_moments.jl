edge_tangents(::Triangle) = (SA[1.0, 0.0], SA[0.0, 1.0], SA[-1.0, 1.0])

function edge_parametrization(::Triangle, vars, s) 
    return (Dict([vars[1] => s, vars[2] => 0]), Dict([vars[1] => 0, vars[2] => s]), Dict([vars[1] => s, vars[2] => 1-s]))
end

face_parametrization(vars, ::Triangle) = 1 - vars[1]


function symbolic_edge_moment(v, t, L, pe, s)
    vl = dot(v, t)
    vl = Symbolics.substitute(vl, pe)
    return integrate(vl * L, (s, 0, 1), symbolic = true, detailed = false)
end

function edge_moments_NI!(eqs, space, element = Triangle(), ptype = LegendreP())
    if dim(element) == 2
        @variables s
        s₀ = 2*s - 1 # Map to (-1, 1)
        V = space.exprs
        vars = space.x
        t = edge_tangents(element)
        p = edge_parametrization(element, vars, s)
        L = [Symbolics.expand(legendrep(i, s₀)) for i in 0:space.degree]
        for (pe, te) in zip(p, t)
            for le in L
                push!(eqs, symbolic_edge_moment(V, te, le, pe, s))
            end
        end
        return nothing
    end
    @error "Not Implemented for $element and $ptype"
end

expsimp(x) = Symbolics.simplify(Symbolics.expand(x))

polynomial_basis(vars, degree, ::DubinerP) = dubiner_basis(vars, degree)

function symbolic_face_moments!(φ, L, pe, vars)
    s, t = vars
    int_t = integrate(dot(φ, L), (t, 0, pe), symbolic = true, detailed = false)
    return simplify(integrate(int_t, (s, 0, 1), symbolic = true, detailed = false))
end

function face_moments_NI!(eqs, space, element = Triangle(), ptype = DubinerP())
    if dim(element) == 2
        @variables s, t
        vars = space.x
        vl = space.exprs
        param = face_parametrization((s,t), element)
        φ = Symbolics.substitute(vl, Dict(vars[1] => t, vars[2] => s))
        mons = polynomial_basis((s, t), space.degree-1, ptype)
        for m in mons
            push!(eqs, symbolic_face_moments!(φ, (0,m), param, (s, t))) 
            push!(eqs, symbolic_face_moments!(φ, (m,0), param, (s, t)))
        end
        return nothing
    end
    @error "Not Implemented for $element and $ptype"
end
