# StepWise-Sampling
Stepwise Transition is an umbrella sampling hamiltonian replica-exchange that consists of an 8step Transition. Al steps together to take the system from state A to state B. here, for transmembrane MHP1, the sampling uses the umbrella windows from the outward-facing open (OF) state to the Outward-facing occluded (OC) state.

The configuration of all the codes are shown by coding.png
NAMD2 runs the Umbrella sampling simulations along with replica exchange of the targeted system.
Along with NAMD2, three other files sequentially are called after each other: 1) sim.conf 2) replica.namd 3) sim_base.conf
replica.namd, reads all the initial coordinates form the initial folder too.
sim_bace.conf, is uses to run the Main_Subcodes/Main.tcl, Main.tcl reads all the initial values of each step transition. Therefore, it decides the system with what window and condition needs to be run based on each replica.  
The initial conditions for each steps are in the folder input.
