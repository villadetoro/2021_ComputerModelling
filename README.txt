PROJECT: SOLAR SYSTEM SIMULATION (s1867522 - Marta-Villa de Toro Sánchez)

To initiate the program the command line takes 6 arguments:

0: SimulationSS.py 
1: file with initial conditions for the bodies
2: name of the file to write the trajectory (XYZ format)
3: name of the file to write the observables computed
4: name of the file to write the energies computed
5: parameters for the simulation

Please use: python3 SimulationSS.py initials.txt tajectory.xyz observables.txt energy.txt parameters.txt

The parameters used are:
numstep: 200000
dt: 1
G: 1.48818517e-34

Time unit: days
Distance unit: AU
Mass unit: kg
Gravitational constant (in consistent units): 1.48818517e-34