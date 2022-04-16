module Trusty

using PaddedViews, Rotations, LinearAlgebra, Graphs, MetaGraphs, InvertedIndices, Unitful

export Section, PlanarTruss, dofdisplacements, deform, strain

unitdofs(dimensions) = PaddedView(0, fill(1, (1, 1)), (dimensions, dimensions))

struct Section
    elasticity::Number
    area::Number
end

Member = Pair{Int,Pair{Section,Int}}

Geometry{N<:Number} = Vector{<:Vector{N}}

function trussgraph(memberpairs::Vector{Member}, vertexlocations::Geometry)
    edgesandsections = Dict(Edge(a => b) => s for (a, (s, b)) in memberpairs)
    edges = keys(edgesandsections)

    graph = MetaGraph(Graph(collect(edges)))

    for (edge, section) in edgesandsections
        set_prop!(graph, edge, :section, section)
    end

    for (vertex, location) in enumerate(vertexlocations)
        set_prop!(graph, vertex, :location, location)
    end

    graph
end

struct DOFComponent
    vertex::Int
    dimension::Int
end

section(graph::MetaGraph, edge::Edge)::Section = get_prop(graph, edge, :section)

location(graph::MetaGraph, vertex::Int)::Vector = get_prop(graph, vertex, :location)

componentindex(vertex::Int, dimension::Int) = 2(vertex - 1) + dimension

edgevec(graph::MetaGraph, (; src, dst)::Edge) = location(graph, dst) - location(graph, src)

function stiffness(graph::MetaGraph, edge::Edge)
    (source, destination) = (edge.src, edge.dst)

    doublevertices = 2nv(graph)
    matrixsize = (doublevertices, doublevertices)

    (; elasticity, area) = section(graph, edge)
    vec = edgevec(graph, edge)
    rotmatrix = RotMatrix(atan(reverse(vec)...))
    vertexstiffness =
        rotmatrix * unitdofs(2) * rotmatrix' * elasticity * area / norm(vec)

    sourceindex = componentindex(source, 1)
    destindex = componentindex(destination, 1)

    soft = zero(eltype(vertexstiffness))

    selfstiffness =
        PaddedView(soft, vertexstiffness, matrixsize, (sourceindex, sourceindex)) +
        PaddedView(soft, vertexstiffness, matrixsize, (destindex, destindex))

    neighborstiffness =
        PaddedView(soft, vertexstiffness, matrixsize, (sourceindex, destindex)) +
        PaddedView(soft, vertexstiffness, matrixsize, (destindex, sourceindex))

    selfstiffness - neighborstiffness
end

stiffness(truss::MetaGraph) = sum(stiffness(truss, edge) for edge in edges(truss))

ConstraintSet = NTuple{2,Vector{Int}}

constrainedindices(constraints::ConstraintSet) =
    [componentindex(vertex, dimension)
     for (dimension, vertices) in enumerate(constraints) for vertex in vertices]

function constrain(matrix::Matrix{<:Number}, constraints::ConstraintSet)
    freeindices = Not(constrainedindices(constraints)...)

    freestiffness = matrix[freeindices, freeindices]

    vertices = 1:(size(matrix)[1]÷2)

    alldofs =
        [DOFComponent(vertex, dimension) for vertex in vertices for dimension in 1:2]

    freedofs = alldofs[freeindices]

    (freestiffness, freedofs)
end

abstract type Structure end

struct PlanarTruss <: Structure
    topology::MetaGraph
    constraints::ConstraintSet

    PlanarTruss(edges::Vector{Member}, geometry::Geometry, constraints::ConstraintSet) =
        new(trussgraph(edges, geometry), constraints)
end

struct StiffnessProblem{S<:Structure,M<:Number}
    structure::S
    stiffness::Matrix{M}
    freedofs::Vector{DOFComponent}
end

Base.convert(::Type{StiffnessProblem}, truss::PlanarTruss) =
    StiffnessProblem(truss, constrain(stiffness(truss.topology), truss.constraints)...)

SparseLoadMap{F<:Number} = Dict{Int,<:Union{Vector{F},Tuple{Vararg{F}}}}

function loadvec(dofs, loads::L) where {L<:SparseLoadMap}
    noload = zeros(eltype(first(values(loads))), 2)

    [get(loads, vertex, noload)[dimension] for (; vertex, dimension) in dofs]
end

rawdisplacements(stiffness::Matrix{<:Real}, loadvector::Vector{<:Real}) = stiffness \ loadvector

function rawdisplacements(stiffness::Matrix{U}, loadvector::Vector{V}) where {U<:Number,V<:Number}

    u = unit(U)
    v = unit(V)

    inv(u) * rawdisplacements(ustrip.(u, stiffness), ustrip.(v, loadvector)) * v
end

SparseDisplacementMap = Dict{DOFComponent,<:Number}

function dofdisplacements((; stiffness, freedofs)::StiffnessProblem, loads::L)::SparseDisplacementMap where {L<:SparseLoadMap}
    displacements = rawdisplacements(stiffness, loadvec(freedofs, loads))

    Dict(zip(freedofs, displacements))
end

dofdisplacements(truss::PlanarTruss, loads::SparseLoadMap) = dofdisplacements(convert(StiffnessProblem, truss), loads)

edgepairs((; topology)::PlanarTruss) =
    [edge.src => section(topology, edge) => edge.dst for edge in edges(topology)]

vertexlocations((; topology)::PlanarTruss) =
    [location(topology, vertex) for vertex in vertices(topology)]

function deform(truss::PlanarTruss,
    displacements::SparseDisplacementMap{D}) where {D}

    (; topology, constraints) = truss

    static = zero(D)

    vertexdisplacements =
        [[get(displacements, DOFComponent(vertex, dimension), static)
          for dimension in 1:2] for vertex in vertices(topology)]

    PlanarTruss(edgepairs(truss), vertexlocations(truss) .+ vertexdisplacements, constraints)
end

deform((; structure)::StiffnessProblem, displacements::SparseDisplacementMap) =
    deform(structure, displacements)

deform(problem::P, displacements::D) where {P<:Union{StiffnessProblem,PlanarTruss},D<:Dict{Int,<:Any}} =
    deform(problem, Dict(k => convert(Vector, v) for (k, v) in displacements))

deform(t::T, load::SparseLoadMap) where {T<:Union{PlanarTruss,StiffnessProblem}} =
    deform(t, dofdisplacements(t, load))

Unitful.ustrip(u::Unitful.Units, truss::PlanarTruss) =
    PlanarTruss(edgepairs(truss), [ustrip.(u, point) for point in vertexlocations(truss)], truss.constraints)

edgevecs(graph::MetaGraph) = edgevec.(Ref(graph), edges(graph))

edgevecs(truss::PlanarTruss) = edgevecs(truss.topology)

strain(original::T, deformed::T) where {T<:Union{PlanarTruss,MetaGraph}} =
    NoUnits.(norm.(edgevecs(deformed)) ./ norm.(edgevecs(original))) .- 1

end # module
