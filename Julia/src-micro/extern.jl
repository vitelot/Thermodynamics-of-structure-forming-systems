"""
 Structures are defined
"""

using CSV, DataFrames, Plots, Statistics;
import Base.copy;

###################################
###       Global variables      ###
# ProgramVersion = v"0.4.2";
Opt = Dict{String,Any}(); # options read from par.ini

Double = Float64;
###################################

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
    N::Int; # total nr of atoms
    T::Double; # temperature
    J::Double; # spin-spin coupling
    H::Double; # external magnetic field
    M::Dict{Int,Int}; # magnetisation: number of up and down spin
    atoms::Set{Atom};
    molecules::Set{Molecule};
end

