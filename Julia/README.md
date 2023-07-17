[figure]: ./multiplicities.png "Transition multiplicities"

# Montecarlo simulation of structure-forming systems

## 1. Description

We organize the Julia scripts in two folders:
- src-macro
- src-micro

In the first one, the system is described by means of macro variables such as magnetization and number of molecules.
In the second one, we follow every atom singularly. This variant is useful when multi-particle generalization are needed. 

## 2. General instructions to run the Julia scripts

### 2.1 Preparation
You have to do the following steps only once.

- **First of all**, you need to install the Julia language binaries. Download Julia from the official source repository: [Julia download](https://julialang.org/downloads/).
- **Second**, enter the `Julia` folder.
- **Third**, run the command `julia --project`. This opens the Julia REPL (read-eval-print loop) and loads the environment.
- **Fourth**, enter the package manager by typing `]`.
- **Fifth**, execute the command `instantiate`. This will load the packages needed by the simulation.
- **Sixth**, exit the package manager with a backspace and exit julia with `^d`.

### 2.2 Running the simulation sequentially
Let's suppose you would like to run the macro simulation.

First, enter the src-macro folder
`cd src-macro`
then execute once
`julia --project=.. main.jl`
Executing the previous command for the first time generates the parameter file `par.ini`, which is used to define the parameters of the model.
This is how a `par.ini` file looks like:
> #key                    value <br>
> ############################# <br>
Ntot        100     # total nr of particles <br>
Hfield      0.0     # external magnetic field <br>
Jcoupling   1.0     # spin-spin coupling <br>
> #############################
initialSpinUpFraction  1.0    # initialize the system with this spin up fraction <br>
fractionMinNrFreeAtoms 0.0    # free atoms cannot drop under this fraction <br>
> ############################# <br>
Tmin        0.0     # minimum temperature <br>
Tstep       0.01     # temperature increasing step (can be negative in case Tmin>Tmax) <br>
Tmax        0.5     # maximum temperature <br>
> ############################# <br>
saveResults 1       # save result into file results.csv <br>
doPlot      0       # plot energy, magnetization and nr of <br> molecules as a function of temperature <br>
goParallel  0       # try multi-threading <br>
>############################# <br>

The next executions of `julia --project=.. main.jl` will read the previous file and run the simulation in the temperature range specified.

**Important:** In the sequential runs, *the initial configuration is taken from the previous temperature state.* 
**This allows us to observe hysteresis cycles.**

### 2.3 Running the multi-threading simulation
First, you need to set the `goParallel` option in the `par.ini` to `1`.

Second, execute the command
`julia -t x --project=.. main.jl`
with `x` set to the number of cores you would like to use.
**Important:** In the in the multi-threading runs, *the initial configuration is always set according to the `initialSpinUpFraction` defined in the `par.ini` file.*
Therefore, **Hysteresis cannot be reproduced.**

## 3. Multiplicity table
In the macro-simulation we use the following transition multiplicities:
![Transition multiplicities][figure]
