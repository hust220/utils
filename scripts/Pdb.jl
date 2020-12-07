# Manipulation of PDB and mol2 files

using Printf

struct PDB end
struct Mol2 end

mutable struct Atom
    name::String
    ind::Int
    bfactor::Float64
    element::String
    charge::Float64
    coords::Vector{Float64}
end
Atom() = Atom("", -1, 0, "", 0, [])
Atom(name::AbstractString) = Atom(name, -1, 0, "", 0, [])
Atom(name::AbstractString, ind::Int) = Atom(name, ind, 0, "", 0, [])
Atom(name::AbstractString, x::Float64, y::Float64, z::Float64) = Atom(name, -1, 0, "", 0, [x, y, z])
Atom(name::AbstractString, ind::Int, x::Float64, y::Float64, z::Float64) = Atom(name, ind, 0, "", 0, [x, y, z])
Base.getindex(atom::Atom, i::Int) = atom.coords[i]
function Base.setindex!(atom::Atom, i::Int, n::Float64)
    atom.coords[i] = n
end
Base.iterate(atom::Atom, state=1) = state > 3 ? nothing : (atom.coords[state], state+1)
Base.length(atom::Atom) = length(atom.coords)
function distance(a1::Atom, a2::Atom)
    sum = 0.0
    for i in 1:3
        d = a1[i] - a2[i]
        sum += d*d
    end
    return sqrt(sum)
end

mutable struct Residue
    name::String
    ind::Int
    atoms::Vector{Atom}
end
Residue() = Residue("", -1, [])
Base.getindex(residue::Residue, i::Int) = residue.atoms[i]
Base.getindex(residue::Residue, atomname::String) = filter(atom->atom.name==atomname, residue.atoms)[1]
function hasatom(residue::Residue, atomname::String)
    return !isempty(filter(atom->atom.name==atomname, residue.atoms))
end
function Base.setindex!(residue::Residue, i::Int, atom::Atom)
    residue.atoms[i] = atom
end
function Base.iterate(residue::Residue, state=1)
    if state > length(residue.atoms) 
        return nothing
    else
        return (residue.atoms[state], state+1)
    end
end
Base.length(residue::Residue) = length(residue.atoms)

mutable struct Chain
    name::String
    residues::Vector{Residue}
end
Chain() = Chain("", [])
Base.getindex(chain::Chain, i::Int) = chain.residues[i]
function Base.setindex!(chain::Chain, i::Int, residue::Residue)
    chain.residues[i] = residue
end
Base.iterate(chain::Chain, state=1) = state > length(chain.residues) ? nothing : (chain.residues[state], state+1)
Base.length(chain::Chain) = length(chain.residues)

mutable struct Model
    chains::Vector{Chain}
end
Model() = Model([])
Base.getindex(model::Model, i::Int) = model.chains[i]
function Base.setindex!(model::Model, i::Int, chain::Chain)
    model.chains[i] = chain
end
Base.iterate(model::Model, state=1) = state > length(model.chains) ? nothing : (model.chains[state], state+1)
Base.length(model::Model) = length(model.chains)
function residues(model::Model)
    rs = []
    for chain in model
        for res in chain
            push!(rs, res)
        end
    end
    return rs
end

mutable struct Structure
    models::Vector{Model}
end
Structure() = Structure([])
Base.getindex(structure::Structure, i::Int) = structure.models[i]
function Base.setindex!(structure::Structure, i::Int, model::Model)
    structure.models[i] = model
end
Base.iterate(structure::Structure, state=1) = state > length(structure.models) ? nothing : (structure.models[state], state+1)
Base.length(structure::Structure) = length(structure.models)

function read(filename::AbstractString, ::Type{Mol2})
    residue = Residue()
    flag = false
    for line in eachline(filename)
        if strip(line) == "@<TRIPOS>ATOM"
            flag = true
        elseif startswith(line, "@<TRIPOS>")
            flag = false
        else
            if flag
                v = split(strip(line), r"\s+")
                #                println(v)
                atom_id = parse(Int, v[1])
                atom_name = v[2]
                x = parse(Float64, v[3])
                y = parse(Float64, v[4])
                z = parse(Float64, v[5])
                atom_type = v[6]
                atom = Atom(atom_name, atom_id, x, y, z)
                push!(residue.atoms, atom)
            end
        end
    end
    return residue
end

function read(filename::AbstractString, ::Type{PDB})
    residue = Residue()
    chain = Chain()
    model = Model()
    structure = Structure()

    function addresidue!()
        if length(residue.atoms) > 0
            push!(chain.residues, residue)
            residue = Residue()
        end
    end

    function addchain!()
        addresidue!()
        if length(chain.residues) > 0
            push!(model.chains, chain)
            chain = Chain()
        end
    end

    function addmodel!()
        addchain!()
        if length(model.chains) > 0
            push!(structure.models, model)
            model = Model()
        end
    end

    res_name = ""
    chain_name = ""
    res_num = -1

    for line in eachline(filename)
        line = rstrip(line)
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            res_name_old = res_name;
            chain_name_old = chain_name;
            res_num_old = res_num;

            atom_num = parse(Int, line[7:11]);
            atom_name = strip(line[13:16]);
            alt_loc = line[17:17];
            res_name = strip(line[18:20]);
            chain_name = strip(line[21:22]);
            res_num = parse(Int, line[23:26]);
            insertion = line[27:27];
            x = parse(Float64, line[31:38]);
            y = parse(Float64, line[39:46]);
            z = parse(Float64, line[47:54]);

            if (length(line) > 60) occupancy = parse(Float64, line[55:60]) end
            if (length(line) > 66) bfactor = parse(Float64, line[61:66]) end
            if (length(line) > 78) element = strip(line[77:78]) end
            if (length(line) > 80) charge = strip(line[79:80]) end

            atom = Atom(atom_name, atom_num, x, y, z)

            if res_name_old != res_name || res_num_old != res_num
                addresidue!()
                if chain_name_old != chain_name
                    addchain!()
                end
            end
            push!(residue.atoms, atom)
            residue.name = res_name
            residue.ind = res_num
            chain.name = chain_name
        elseif startswith(line, "MODEL")
            addmodel!()
        elseif startswith(line, "ENDMDL")
            addmodel!()
        elseif line == "TER"
            addchain!()
        elseif line == "END"
            break
        end
    end

    addmodel!()
    #    println(structure)
    return structure
end

function center(residue::Residue)
    return reduce((c, a)->map(i->c[i]+a[i], 1:3), residue; init=zeros(3)) / length(residue)
end
