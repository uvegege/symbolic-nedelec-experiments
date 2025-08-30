using Symbolics
using SymbolicNumericIntegration
using ClassicalOrthogonalPolynomials
using LinearAlgebra: dot, isdiag
using StaticArrays
using CairoMakie
#using GLMakie
using Latexify

include("./monomials.jl")
include("./N1kind.jl")
include("./N1Spaces.jl")
include("./symbolic_moments.jl")
include("./N1Poly.jl")
include("./outputs.jl")
include("./verifications.jl")
include("./vizN.jl")