module Trusty

using PaddedViews, Rotations, LinearAlgebra, Graphs, MetaGraphs, InvertedIndices, Unitful, SparseArrays, AngleBetweenVectors, Plots

export Section, Truss, deform, strain, volume, DOFComponent

unitdofs(dimensions) = PaddedView(0, fill(1, (1, 1)), (dimensions, dimensions))

"""
    Section

Section properties of a `Member`.

    Section(elasticity::Number, area::Number)

Construct `Section` object with `elasticity` and `area`.
"""
struct Section
    elasticity::Number
    area::Number
end

"""
	source::Int => section::Section => destination::Int

A truss member with a specified `section` connecting vertices `source` and `destination`.
"""
const Member = Pair{Int,Pair{Section,Int}}

"""
	AbstractVertexMatrix{T,D<:Unitful.Dimensions}

Matrix of vertices of element type `T`, possibly on unit dimensions `D`. Each column is a vertex's vector.

# Examples
```jldoctest
julia> V = let v₁ = [1; 1]
	       v₂ = [2; 2]
	       v₃ = [3; 3]
	
	       [v₁ v₂ v₃]
	   end
2×3 Matrix{Int64}:
 1  2  3
 1  2  3

julia> V isa AbstractVertexMatrix
true

julia> using Unitful

julia> U = let v₁ = [1; 1]u"inch"
	       v₂ = [2; 2]u"cm"
	       v₃ = [3; 3]u"m"
	
	       [v₁ v₂ v₃] # all vectors in dimension 𝐋 (length)
	   end
2×3 Matrix{Quantity{Rational{Int64}, 𝐋 , Unitful.FreeUnits{(m,), 𝐋 , nothing}}}:
 127//5000 m  1//50 m  3//1 m
 127//5000 m  1//50 m  3//1 m

julia> U isa AbstractVertexMatrix
true
```
"""
const AbstractVertexMatrix{T,D} = M where {M<:AbstractMatrix{<:Union{T,Quantity{T,D,U}}}} where {T<:Real,D,U}

function trussgraph(memberpairs::Vector{Member}, vertexlocations::AbstractVertexMatrix)
    edgesandsections = Dict(Edge(a => b) => s for (a, (s, b)) in memberpairs)
    edges = keys(edgesandsections)

    graph = MetaGraph(Graph(collect(edges)))

    for (edge, section) in edgesandsections
        set_prop!(graph, edge, :section, section)
    end

    for (vertex, location) in enumerate(eachcol(vertexlocations))
        set_prop!(graph, vertex, :location, location)
    end

    graph
end

"""
    DOFComponent

Key for referring to a component of a DOF along some dimension.

    DOFComponent(vertex::Int, dimension::Int)

Construct `DOFComponent` representing freedom along `dimension` at `vertex`.
"""
struct DOFComponent
    vertex::Int
    dimension::Int
end

abstract type Structure end

const ConstraintSet = Tuple{Vararg{Vector{Int}}}

"""
	Truss{N}

Truss in `N` dimensions.

	Truss{N}(edges::Vector{Member}, vertices::AbstractVertexMatrix, constraints::ConstraintSet) where N

Construct a `Truss` in `N` dimensions from `Member`s, vertices, and constraints.

See also [`Member`](@ref), [`AbstractVertexMatrix`](@ref), and [`ConstraintSet`](@ref).

	Truss(edges::Vector{Member}, vertices::AbstractVertexMatrix, constraints::ConstraintSet)

Construct a `Truss`, setting the number of dimensions `N` at runtime instead of compile time.

# Examples
```jldoctest
julia> Truss(1:2 .=> Section(1, 1) .=> 2:3,
	     [0 1 2
	      0 1 0],
	     ([1, 3], [1, 3]))
```
"""
struct Truss{N} <: Structure
    topology::MetaGraph
    constraints::ConstraintSet

    function Truss{N}(edges::Vector{Member}, vertices::AbstractVertexMatrix, constraints::ConstraintSet) where {N}
        @assert N isa Int

        new{N}(trussgraph(edges, vertices), constraints)
    end

    function Truss(edges::Vector{Member}, vertices::AbstractVertexMatrix, constraints::ConstraintSet)
        Truss{size(vertices)[1]}(edges, vertices, constraints)
    end
end

section((; topology)::Truss, edge::Edge)::Section = get_prop(topology, edge, :section)

location((; topology)::Truss, vertex::Int)::Vector = get_prop(topology, vertex, :location)

componentindex(ndims, vertex::Int, dimension::Int) = ndims * (vertex - 1) + dimension

edgevec(truss::Truss, (; src, dst)::Edge) = location(truss, dst) - location(truss, src)

const AbstractStiffnessMatrix{T} = AbstractMatrix{T} where T<:Number

"""
	aligndofs(ndims::Val{N}, dofs::AbstractMatrix, vector::Vector)

Rotate a member's DOF matrix `dofs` in `N` dimensions to align to some `vector` (in practice, the vector between the member's two vertices).
"""
function aligndofs(ndims::Val{2}, dofs::AbstractMatrix, vector::Vector)
    rotmatrix = RotMatrix(atan(reverse(vector)...))

    rotmatrix * dofs * rotmatrix'
end

function aligndofs(ndims::Val{3}, dofs::AbstractMatrix, vector::Vector)
    u = unit(upreferred(first(vector)))
    unitless = ustrip.(u, vector)

    rotmatrix = Matrix(rotation_between([1, 0, 0], unitless))

    rotmatrix * dofs * rotmatrix'
end

function stiffness(truss::Truss{N}, edge::Edge)::SparseMatrixCSC where {N}
    (source, destination) = (edge.src, edge.dst)
    components = N * nv(truss.topology)
    matrixsize = (components, components)

    (; elasticity, area) = section(truss, edge)
    vec = edgevec(truss, edge)

    vertexstiffness =
        aligndofs(Val(N), unitdofs(N), vec) * elasticity * area / norm(vec)

    sourceindex = componentindex(N, source, 1)
    destindex = componentindex(N, destination, 1)

    soft = zero(eltype(vertexstiffness))

    selfstiffness =
        PaddedView(soft, vertexstiffness, matrixsize, (sourceindex, sourceindex)) +
        PaddedView(soft, vertexstiffness, matrixsize, (destindex, destindex))

    neighborstiffness =
        PaddedView(soft, vertexstiffness, matrixsize, (sourceindex, destindex)) +
        PaddedView(soft, vertexstiffness, matrixsize, (destindex, sourceindex))

    selfstiffness - neighborstiffness
end

stiffness(truss::Truss) = sum(stiffness(truss, edge) for edge in edges(truss.topology))

constrainedindices(ndims::Int, constraints::ConstraintSet) =
    [componentindex(ndims, vertex, dimension)
     for (dimension, vertices) in enumerate(constraints) for vertex in vertices]

function constrain(ndims::Int, matrix::AbstractStiffnessMatrix, constraints::ConstraintSet)
    freeindices = Not(constrainedindices(ndims, constraints)...)

    freestiffness = matrix[freeindices, freeindices]

    vertices = 1:(size(matrix)[1]÷ndims)

    alldofs =
        [DOFComponent(vertex, dimension) for vertex in vertices for dimension in 1:ndims]

    freedofs = alldofs[freeindices]

    (freestiffness, freedofs)
end

struct StiffnessProblem{S<:Structure,K<:AbstractStiffnessMatrix}
    structure::S
    stiffness::K
    freedofs::Vector{DOFComponent}
end

Base.convert(::Type{StiffnessProblem}, truss::Truss{N}) where {N} =
    StiffnessProblem(truss, constrain(N, stiffness(truss), truss.constraints)...)

Base.convert(::Type{Truss}, problem::StiffnessProblem{<:Truss}) = problem.structure

"""
	SparseLoadMap{F<:Number}

A `Dict` of mappings between vertices (`Int`) and force vectors on those vertices with elements of type `F`.

# Examples
```jldoctest
julia> f = Dict(1 => [1; 0], 4 => [2, 3])
Dict{Int64, Vector{Int64}} with 2 entries:
  4 => [2, 3]
  1 => [1, 0]

julia> f isa SparseLoadMap
true
```
"""
const SparseLoadMap{F<:Number} = Dict{Int,<:Union{Vector{F},Tuple{Vararg{F}}}}

function loadvec(ndims::Int, dofs::Vector{DOFComponent}, loads::L) where {L<:SparseLoadMap}
    noload = zeros(eltype(first(values(loads))), ndims)

    [get(loads, vertex, noload)[dimension] for (; vertex, dimension) in dofs]
end

rawdisplacements(stiffness::AbstractMatrix{<:Real}, loadvector::Vector) = stiffness \ loadvector

function rawdisplacements(stiffness::AbstractStiffnessMatrix{U}, loadvector::Vector{V}) where {U,V}
    u = unit(U)
    v = unit(V)

    inv(u) * rawdisplacements(ustrip.(u, stiffness), ustrip.(v, loadvector)) * v
end

const SparseDisplacementMap = Dict{DOFComponent,<:Number}

function dofdisplacements((; stiffness, freedofs)::StiffnessProblem{Truss{N}}, loads::SparseLoadMap)::SparseDisplacementMap where {N}
    displacements = rawdisplacements(stiffness, loadvec(N, freedofs, loads))

    Dict(zip(freedofs, displacements))
end

dofdisplacements(truss::Truss, loads::SparseLoadMap) = dofdisplacements(convert(StiffnessProblem, truss), loads)

edgepairs(truss::Truss) =
    [edge.src => section(truss, edge) => edge.dst for edge in edges(truss.topology)]

vertexlocations(truss::Truss{N}) where {N} =
    [location(truss, vertex)[dimension] for dimension in 1:N, vertex in vertices(truss.topology)]

function deform(truss::Truss{N}, displacements::SparseDisplacementMap{D}) where {N,D}
    (; topology, constraints) = truss

    static = zero(D)

    vertexdisplacements =
        [get(displacements, DOFComponent(vertex, dimension), static)
         for dimension in 1:N, vertex in vertices(topology)]

    Truss{N}(edgepairs(truss), vertexlocations(truss) + vertexdisplacements, constraints)
end

deform((; structure)::StiffnessProblem, displacements::SparseDisplacementMap) =
    deform(structure, displacements)

deform(problem::P, displacements::D) where {P<:Union{StiffnessProblem,Truss},D<:Dict{Int,<:Any}} =
    deform(problem, Dict(k => convert(Vector, v) for (k, v) in displacements))

"""
	deform(truss::Truss{N}, load::SparseLoadMap) where {N}

Return a truss whose vertex locations are those of `truss` after deforming subject to some `load` scenario.

See also [`SparseLoadMap`](@ref).
"""
deform(t::T, load::SparseLoadMap) where {T<:Union{Truss,StiffnessProblem}} =
    deform(t, dofdisplacements(t, load))

Unitful.ustrip(u::Unitful.Units, truss::Truss) =
    Truss(edgepairs(truss), [ustrip.(u, point) for point in vertexlocations(truss)], truss.constraints)

edgevecs(truss::Truss) = edgevec.(truss, edges(truss.topology))

"""
	strain(original::Truss, deformed::Truss)

Return a `Dict` mapping a truss' `Member`s to the amount of strain in those members between the truss's `original` shape and its `deformed` shape.
"""
strain(original::Truss, deformed::Truss) =
    edgepairs(original) .=>
        NoUnits.(norm.(edgevecs(deformed)) ./ norm.(edgevecs(original))) .- 1

"""
	volume(truss::Truss, member::Member)

Return the physical volume of `member` in `truss`.
"""
function volume(truss::Truss, (a, ((; area), b))::Member)
    length = norm(edgevec(truss, Edge(a, b)))

    length * area
end

Broadcast.broadcastable(section::Section) = Ref(section)

Broadcast.broadcastable(truss::Truss) = Ref(truss)

"See [`plottruss`](@ref)."
function plottruss!(p::Plots.Plot, truss::Truss; vertexlabels=true,
    textcolor="black",
    width=5,
    label=nothing,
    edgecolor=:auto,
    kw...)

    vertexlocs = collect.(zip(location.(truss,
        vertices(truss.topology))...))

    edges = [collect.(zip(location.(truss, [a, b])...))
             for (a, (_, b)) in edgepairs(truss)]

    with(label=label, kw...) do

        for (edge, color) in Pair.(edges, edgecolor)
            plot!(p, edge..., width=width, linecolor=color)
        end

        scatter!(p, vertexlocs...,
            series_annotations=
            (vertexlabels in (false, nothing)
             ? nothing
             : text.(vertices(truss.topology),
                :bottom,
                (vertexlabels == true
                 ? :auto
                 : vertexlabels),
                color=textcolor)))
    end

    p
end

"""
	plottruss(truss::Truss; kw...)
	plottruss!(p::Plot, truss::Truss; kw...)

Plot a `truss`, optionally atop an existing plot `p`.

Only supports `truss`es with unitless vertex locations. If you use units via `Unitful` to define some truss `t`, decide which units `u` to plot in and use `ustrip(u, t)`.

# Arguments
- Trusty-specific:
  - `vertexlabels=true`: whether to plot vertex labels, or the labels' size
  - `textcolor=:black`: the color of text like vertex labels
  - `width=5`: the width of members
  - `edgecolor=:auto`: the color(s) of edges
- Others:
  - Any keyword argument supported by Plots
"""
plottruss(truss::Truss; kw...) = plottruss!(plot(; kw...), truss; kw...)

end # module
