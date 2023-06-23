function initializeBox(N::Int)
    Atoms = Atom[];
    for i=1:N
        push!(Atoms, Atom(i,rand([-1,1])));
    end
    allspins = [x.σ for x in Atoms];
    nup = count(x->x.σ>0, Atoms);
    ndn = count(x->x.σ<0, Atoms);

    B = Box(Dict(-1=>ndn, 1=>nup),Set(Atoms), Set{Molecule}());
    
    return B;
end

function energy(B::Box)
    Es = Eh = 0.0;
    Atoms = collect(B.atoms);
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
    Avrg  = Real[];

    en = energy(B);
    sumEnergy = en;

    for i in 1:Steps # do at max Step steps
        
        en += spinFlipMove(B);
        sumEnergy += en; 
        # push!(Elist, en);
        if i%avrgStep == 0
            avrg = sumEnergy/i;
            # avrg = sum(Elist[end-avrgStep+1:end])/avrgStep;
            push!(Avrg, avrg);
        end
    end
    return Avrg;
end

"""
    checks whether a spin flip occurs; does it and returns the ΔE 
"""
function spinFlipMove(B::Box)
    a = rand(B.atoms);
    ΔE = spinFlipDeltaEnergy(B,a);
    if ΔE < 0
        B.M[a.σ] -= 1; # decrease number of old spins
        a.σ = -a.σ;
        B.M[a.σ] += 1; # and increase the new one
        return ΔE; # spin flip accepted
    else
        p = rand();
        if p<exp(-ΔE/T)
            B.M[a.σ] -= 1;
            a.σ = -a.σ;
            B.M[a.σ] += 1;
            return ΔE; # spin flip accepted
        end
    end
    return 0.0; # flipping discarded
end

function moleculeMove(B::Box)
    if rand() < 0.5
        moleculeSplit(B);
    else
        moleculeJoin(B);
    end
end

"""
    checks whether two atoms join; does it and returns the ΔE 
"""
function moleculeJoin(B::Box)
    Atoms = B.atoms;
    if length(Atoms) < 2
        @info "No atoms available to join";
        return 0.0;
    end

    # select two random different atoms
    a1 = rand(Atoms);
    a2 = rand(Atoms);
    while true
        a2.id == a1.id || break;
        a2 = rand(Atoms);
    end
    # removing an atom corresponds to half a spin flip
    ΔE = 0.5*(spinFlipDeltaEnergy(B,a1)+spinFlipDeltaEnergy(B,a2));
    
    # current spin population
    nup = B.M[1];
    ndn = B.M[-1];
    # spin population after removal
    nup_tilde = nup - (a1.σ + 1)÷2 - (a2.σ + 1)÷2;
    ndn_tilde = ndn - (1 - a1.σ)÷2 - (1 - a2.σ)÷2;
    # current nr of molecules
    nmol = length(B.molecules);

    prob = factorialRatio(nup,nup_tilde)*
           factorialRatio(ndn,ndn_tilde)/
           (2*nmol+2)*
           exp(-ΔE/T);

    if rand() < prob
        # accept the joining
        # remove the two atoms from box1
        pop!(Atoms, a1);
        pop!(Atoms, a2);
        mol = Molecule([a1,a2]);
        # add the new molecule into box2
        push!(B.molecules, mol);
        B.M[1] = nup_tilde;
        B.M[-1]= ndn_tilde;
        return ΔE;
    else
        # do not join and do nothing
        return 0.0;
    end
end

"""
calculates the ratio n!/m!
"""
function factorialRatio(n::Int, m::Int)::Real
    n==m && return 1.0;
    if n==0 || n==1 
        return 1.0/factorial(m);
    end
    if m==0 || m==1 
        return factorial(n);
    end
    f=1.0;
    if n>m
        for i=m+1:n
            f *= i;
        end
        return f;
    else
        for i = n+1:m
            f /= i
        end
        return f;
    end
end