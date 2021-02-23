# StepWise-Sampling
Stepwise Transition is an all atom umbrella sampling hamiltonian replica-exchange simulation on transmembrane protein MHP1 that consists of 8 transition steps. All steps cause the system go through a transition from state A to state B. Here, the target is a transmembrane protein name as MHP1. The sampling uses the umbrella windows from the outward-facing open (OF) state to the Outward-facing occluded (OC) state.

1) The schematics of all the codes are shown by the picture below.
2) NAMD2 runs the Umbrella sampling simulations along with replica exchange for the targeted system.
3) Along with NAMD2, three other files sequentially are called after each other: 1) sim.conf 2) replica.namd 3) sim_base.conf.
4) replica.namd, reads all the initial coordinates form the initial folder too.
5) sim_bace.conf which is the configuration file, runs the Main_Subcodes/Main.tcl. The Main.tcl file reads all the initial values of each step transition. Therefore, based on each replica, it decides the system with what window and condition needs to be run.  
6) inside the Main_Subcodes folder, there are 8 other codes as 8 different transition steps which have all the essential commands that measures the condition of each transition. (See G358 for more detail).
7) initial conditions for each step transition are in the folder input.


![alt text](https://github.com/hamedmeshkin/StepWise-Transition-Model/blob/main/schematic%20diagram%20of%20the%20codes.png?raw=true)
