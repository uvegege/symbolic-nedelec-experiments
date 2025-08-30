using Latexify

function roundcoefficients(x; digits = 6)
    parts = Symbolics.terms(x)
    fac = Symbolics.factors.(parts)
    l = map(fac) do f
        v = isa(f[1].val, SymbolicUtils.Symbolic) ? f[1] : round(f[1].val, digits = digits)
        if length(f) > 1
            reduce(*, f[2:end]; init = v)
        else
            v
        end
    end
    Symbolics.simplify(sum(l))
end


function replace_u_xy(N)
    @variables u[1:2]
    @variables x, y
    ds = Dict(u[1] => x, u[2] => y)
    map(N) do Ni
        Symbolics.substitute.(Ni, Ref(ds))
    end
end

function create_function_file(name, Nsymbolic)

    @variables u[1:2]
    @variables x, y
    ds = Dict(u[1] => x, u[2] => y)
    Nsubs = map(Nsymbolic) do Ni
        Symbolics.substitute.(Ni, Ref(ds))
    end
    buff = IOBuffer()

    write(buff, "function ")
    write(buff, name)
    write(buff, "(element, u, i) \n")
    write(buff, "   x, y = u \n")

    for i in eachindex(Nsubs)
        write(buff, "   ")
        write(buff, string(i))
        write(buff, " == i && return @SVector[")
        #write(buff, join(string.(roundcoefficients.(Nsubs[i])), ", "))
        write(buff, join(string.(Nsubs[i]), ", "))
        write(buff, "]\n")
    end
    write(buff, "   return nothing\n")
    write(buff, "end")
    write("$name.jl", String(take!(buff)))

    return nothing
end

# https://defelement.org/elements/examples/triangle-nedelec1-legendre-1.html
function N1DefElement()
    @variables x, y

    #n0 = [y*(-4x - 4y + 3), x*(4x+4y-3)]
    #n1 = [sqrt(3)*y*(4x-4y+1)/3, sqrt(3)*x*(-4x+4y+1)/3]
    #n2 = [y*(1-4x), 4*x^2-5x+1]
    #n3 = [sqrt(3)*y * (4x+8y-5)/3, sqrt(3)*(-4x^2-8x*y + 7x +6y - 3)/3]
    #n4 = [4*y^2 - 5*y + 1, x*(1-4y)]
    #n5 = [sqrt(3)/3 * (-8*x*y +6x -4y^2+7y-3), sqrt(3)/3 * x * (8x+4y-5)]
    #n6 = [4*sqrt(2) * y * (-x-2y+2), 4*sqrt(2)*x*(x+2y-1)]
    #n7 = [4*sqrt(2)* y * (2x+y-1), 4*sqrt(2) * x * (-2x-y+2)]

    n0 = [y*(-4x - 4y + 3), x*(4x+4y-3)]
    n1 = [y*(4x-4y+1), x*(-4x+4y+1)]
    n2 = [y*(1-4x), 4*x^2-5x+1]
    n3 = [y * (4x+8y-5), (-4x^2-8x*y + 7x +6y - 3)]
    n4 = [4*y^2 - 5*y + 1, x*(1-4y)]
    n5 = [(-8*x*y +6x -4y^2+7y-3),  x * (8x+4y-5)]
    n6 = [4*y * (-x-2y+2), 4*x*(x+2y-1)]
    n7 = [4*y * (2x+y-1),  4*x * (-2x-y+2)]

    Nde = [n4, n5, n2, n3, n0, n1, n7, n6]

    map(Nde) do Ni
        Symbolics.expand.(Ni)
    end
end

function NTest()
    N1 = N1k(1)
    N1xy = replace_u_xy(N1)
    N1xyde = N1DefElement()
    for idx in eachindex(N1)
        println(" ----- Basis Function -----")
        for idx2 in eachindex(N1[idx])
            c = Symbolics.simplify(Symbolics.expand(roundcoefficients(N1xyde[idx][idx2])/roundcoefficients(N1xy[idx][idx2])))
            println("N$idx,$idx2 DefElement / N$idx,$idx2 =  $c")
        end
    end
end


function to_latex_string(N)
    for (i, Ni) in enumerate(replace_u_xy(N))
        nistring = replace(latexify(Ni), '\n' => "", "\\end{equation}" => "", "\\begin{equation}" => "")
        if i == 1 
            println("\\begin{array}{c} \\phi_{$i} =", nistring, " \\\\")
        else
            println("\\phi_{$i} =", nistring, " \\\\")
        end
        i == length(N) && println("\\end{array}")
    end
end
