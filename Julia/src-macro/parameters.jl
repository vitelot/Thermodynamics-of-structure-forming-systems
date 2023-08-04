function loadParameters(file::String)
    
    if !isfile(file)
        createIniFile(file)
    end

    for line in eachline(file)
        occursin(r"^#", line) && continue # ignore lines beginning with #
        df = split(line, r"\s+")
        length(df) < 2 && continue # ignore empty lines
        key = df[1] ; val = df[2]
        ####################################################################
        if(key=="Version")                    Opt[key] = val
        ####################################################################
        elseif(key=="Ntot")                   Opt[key] = parse(Int, val)
        elseif(key=="Hfield")                 Opt[key] = parse(Double, val)
        elseif(key=="Jcoupling")              Opt[key] = parse(Double, val)
        elseif(key=="MolEnFormation")         Opt[key] = parse(Double, val)
        elseif(key=="initialSpinUpFraction")  Opt[key] = parse(Double, val)
        elseif(key=="fractionMinNrFreeAtoms") Opt[key] = parse(Double, val)
        ####################################################################
        elseif(key=="thermalisationSteps")  Opt[key] = round(Int, parse(Float64, val))
        elseif(key=="Steps")                Opt[key] = round(Int, parse(Float64, val))
        elseif(key=="avrgStep")             Opt[key] = round(Int, parse(Float64, val))
        ####################################################################
        elseif(key=="Tmin")                 Opt[key] = parse(Double, val)
        elseif(key=="Tstep")                Opt[key] = parse(Double, val)
        elseif(key=="Tmax")                 Opt[key] = parse(Double, val)
        ####################################################################
        elseif(key=="verbosity")            Opt[key] = parse(Int, val)
        elseif(key=="onlyAccepted")         Opt[key] = parse(Bool, val)
        elseif(key=="saveResults")          Opt[key] = parse(Bool, val)
        elseif(key=="doPlot")               Opt[key] = parse(Bool, val)
        elseif(key=="goParallel")           Opt[key] = parse(Bool, val)
        ####################################################################
        else @warn "WARNING: input parameter $key does not exist";
        end
    end

    if VersionNumber(Opt["Version"]) != ProgramVersion
        println("""
                The par.ini file corresponds to an older version of the program.
                Delete it and rerun the simulation to create a new one.
                Version found: $(Opt["Version"]) --- Current version: $ProgramVersion
                """);
        exit(1);
    end

    ts = Opt["thermalisationSteps"];
    Opt["thermalisationSteps"] = ifelse(ts>0, ts, Opt["Ntot"]*Opt["Ntot"]÷2); # run these nr of steps at the beginning discarding the energy

    st = Opt["Steps"];
    Opt["Steps"]              = ifelse(st>0, st, 100*Opt["thermalisationSteps"]); # total nr of steps

    av = Opt["avrgStep"];
    Opt["avrgStep"]           = ifelse(av>0, av, floor(Int, sqrt(Opt["Steps"]))); # used to estimate the average energy

    @info("Parameters loaded, starting the program.")
end

function createIniFile(file::String)

    INI = open(file, "w")
        print(INI,
"""
#key                    value
#############################
Version     $ProgramVersion   # Program's version
#############################
Ntot        200     # total nr of particles
Hfield      0.0     # external magnetic field
Jcoupling   1.0     # spin-spin coupling
MolEnFormation 0.0  # molecule energy formation; can be of both signs
#############################
initialSpinUpFraction  0.5    # initialize the system with this spin up fraction
fractionMinNrFreeAtoms 0.0    # free atoms cannot drop under this fraction
#############################
thermalisationSteps      0    # if <= 0 use the default: Ntot^2/2; (scientific notation allowed, e.g., 1e6)
Steps                    0    # steps used for each temperature; if <= 0 use the default: 100*thermalisationSteps; (scientific notation allowed, e.g., 1e6)
avrgStep                 0    # steps interval used to calculate averages; if <= 0 use the default: sqrt(Steps); (scientific notation allowed, e.g., 1e6)
#############################
Tmin        0.0     # minimum temperature
Tstep       0.2     # temperature increasing step (can be negative in case Tmin>Tmax)
Tmax        2.0     # maximum temperature
#############################
verbosity    0       # level of info output
onlyAccepted 0       # step counter increases only if the step is accepted
saveResults  0       # save result into file results.csv
doPlot       0       # plot energy, magnetization and nr of molecules as a function of temperature
goParallel   0       # try multi-threading
#############################
"""
)
    close(INI);
    println("Parameter file \"$file\" was missing and a default one was created.\nPlease edit it and rerun.");
    exit();

end
