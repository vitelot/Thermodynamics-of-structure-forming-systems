function initializeBox(N::Int)
    B = Box(Atom[], Molecule[]);

    for i=1:N
        push!(B.atoms, Atom(i,1,true));
    end
    return B;
end

function energy(B::Box)
    Es = Eh = 0.0;
    Atoms = B.atoms;
    n = length(Atoms);

    if H > 0.0 # well, it doesn't really matter ... it's so fast anyway
        for i in eachindex(Atoms)
            si = Atoms[i].σ;
            # Eh += si;
            for j in 1:i-1
                Es += si * Atoms[j].σ;
            end
        end
    else
        for i in eachindex(Atoms)
            si = Atoms[i].σ;
            Eh += si;
            for j in 1:i-1
                Es += si * Atoms[j].σ;
            end
        end
    end

    return -Eh * H - Es * 2.0 * J /(n-1);
end

function spinFlipDeltaEnergy(B::Box, atom::Atom)
    E = 0.0;
    n = length(B.atoms);
    for a in B.atoms
        E += a.σ;
    end 
    E -= atom.σ;

    E = (4.0 * J * E/(n-1) + 2.0 * H) * atom.σ;
    return E;
end

function montecarlo(B::Box)
    Elist = Real[];
    Atoms = B.atoms;

    for _ in 1:10000
        a = rand(Atoms);
        ΔE = spinFlipDeltaEnergy(B,a);
        if ΔE < 0
            a.σ = -a.σ;
        else
            p = rand();
            if p<exp(-ΔE/T)
                a.σ = -a.σ;
            end
        end
        push!(Elist, energy(B));
    end
    return Elist;
end

