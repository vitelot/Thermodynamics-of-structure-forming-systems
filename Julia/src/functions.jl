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

    if H > 0.0
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




