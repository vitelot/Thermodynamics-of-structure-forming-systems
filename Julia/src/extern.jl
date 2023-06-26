"""
 Structures are defined
"""

using CSV, DataFrames, Plots;

Double = Float64;

# structures
mutable struct Atom
    id::Int; # identification number
    σ::Int;  # spin (±1, for the moment)
    # isfree::Bool; # true if this atom is not combined
end

mutable struct Molecule
    # natoms::Int; # number of atoms composing the molecule
    atoms::Vector{Atom}; # the atoms composing the molecule
end

mutable struct Box # the simulation box
    M::Dict{Int,Int}; # magnetisation: number of up and down spin
    atoms::Set{Atom};
    molecules::Set{Molecule};
end

