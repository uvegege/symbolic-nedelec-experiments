function integrate_face_moments!(eqs, φ, vars, degree, ::Triangle, ::DubinerP)
    degree <= 0 && return
    mons = dubiner_basis(vars, degree-1)
    s, t = vars
    for m in mons
        int_t = simplify(integrate(dot(φ, (0, m)), (t, 0, 1 - s), symbolic = true, detailed = false); expand = true)
        if all(<(1.0e-10), Symbolics.coeff.(roundcoefficients(int_t), generar_monomios(s, degree+1)))
            push!(eqs, 0)
        else
            push!(eqs, simplify(integrate(int_t, (s, 0, 1), symbolic = true, detailed = false); expand = true))
        end

        int_t = simplify(integrate(dot(φ, (m, 0)), (t, 0, 1 - s), symbolic = true, detailed = false); expand = true)
        if all(<(1.0e-10), Symbolics.coeff.(roundcoefficients(int_t), generar_monomios(s, degree+1)))
            push!(eqs, 0)
        else
            push!(eqs, simplify(integrate(int_t, (s, 0, 1), symbolic = true, detailed = false); expand = true))
        end
    end
    return nothing
end

function integrate_edge_moments!(eqs, vl, vars, degree, ::Triangle, ::LegendreP)
    s₀, s = vars
    for ᵢ in 0:degree
        push!(eqs, simplify(integrate(roundcoefficients(simplify(vl * legendrep(ᵢ, s₀); expand = true)), (s, 0, 1), symbolic = true, detailed = false)))
    end
    return nothing
end


function dofs_edges(Nj, vars, degree; ptype = LegendreP())
    @variables s
    s₀ = 2*s - 1 # Map to (-1, 1)
    t = edge_tangents(Triangle())
    p = edge_parametrization(Triangle(), vars, s)
    eqs = Num[]
    for (pe, te) in zip(p, t)
        vl = dot(Nj, te)
        vl = Symbolics.substitute(vl, pe)
        integrate_edge_moments!(eqs, vl, (s₀, s), degree, Triangle(), ptype)
    end
    return eqs
end


function dofs_faces(Nj, vars, degree; ptype = DubinerP())
    @variables s, t
    eqs = Num[]
    φ = Symbolics.substitute.(Nj, Ref(Dict(vars[1] => t, vars[2] => s))) 
    integrate_face_moments!(eqs, φ, (s, t), degree, Triangle(), ptype)
    return eqs
end


function verify(Nk, k)
    @variables u[1:2]
    p1 = reduce(hcat, dofs_edges.(Nk, Ref(u), k))'
    p2 = reduce(hcat, dofs_faces.(Nk, Ref(u), k))'
    Mk = simplify.(Symbolics.value.([p1  p2]))
    println("isdiag(M) = $(isdiag(Mk))")
    println("dof_i(Nj) = ")
    Mk = round.(Mk; digits = 8)
end


function curl_2d(F::Vector, vars)
    @assert length(F) == 2 && length(vars) == 2 
    Fx, Fy = F
    x, y = vars
    ∂_x = Differential(x)
    ∂_y = Differential(y)
    return expand_derivatives(∂_x(Fy)) - expand_derivatives(∂_y(Fx))
end



"""

## Example

```
    @variables u[1:2]
    N2 = N1k(2)
    verify_curl_degree(N2, u, 2)
```
"""
function verify_curl_degree(X, vars, k; detailed = true)
    if detailed == true
        xy = replace_u_xy.(X)
        @variables x, y
        i = 1
        for x_ in xy
            println("∇×N$i = ", curl_2d(x_, (x, y)))
            i+=1
        end
    end
    return all((Symbolics.degree(curl_2d(x, vars)) == k for x in X))
end



"""

    @variables u[1:2]
    k = 2
    Nk = N1k(k)
    test_edges(Nk, u, detailed = true)

"""
function test_edges(N, vars; detailed = true)

    @variables s
    t = edge_tangents(Triangle())
    p = edge_parametrization(Triangle(), vars, s)
    x, y = vars
    i = 1
    dims = length(vars)

    l = length(N)
    k = Int((-4 + sqrt(4^2 - 4*(3-l)*1))/(2*1)) + 1
    

    tangential_continue = true
    for (ti, pe) in zip(t, p)
        needed_degres = Set(0:k-1)
        #detailed && println("Current tangential continuity: $tangential_continue")
        detailed && println("--- Edge $i with tangent = ($(ti[1]), $(ti[2])) ---")
        for (idx,ni) in enumerate(N)
            ni_tang = Symbolics.simplify(Symbolics.substitute(dot(ni, ti), pe); expand = true)
            #ni_tang = Num(roundcoefficients(ni_tang, digits = 8))
            detailed && println("N$idx ⋅ ($(ti[1]), $(ti[2])) = $ni_tang")
            if i == ceil(Int, idx/k) 
                if isa(ni_tang.val, SymbolicUtils.Symbolic)
                    if (Symbolics.degree(ni_tang) in needed_degres)
                        pop!(needed_degres, Symbolics.degree(ni_tang))
                    else
                        tangential_continue *= false
                    end
                else
                    tangential_continue *= ni_tang.val == 1 ? true : false
                    tangential_continue && pop!(needed_degres, 0)
                end
            else
                if isa(ni_tang.val, SymbolicUtils.Symbolic) 
                    tangential_continue *= false
                else
                    tangential_continue *= Float64(ni_tang.val) == 0 ? true : false
                end
            end
        end
        i += 1
    end
    return tangential_continue
end