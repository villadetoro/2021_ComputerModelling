
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan  23 13:34:01 2021

Solar Syatem simulation, N-body system interacting through Newtonian gravity

The code simulate the motion of the main bodies of the solar system, starting 
from realistic initial conditions. It also provides some observables (aapsides, 
period orbits and energies)

@author: Marta-Villa de Toro Sánchez s1867522


"""


import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D
from scipy import signal



def force_gravitational(G, p3d_list, separation_array, moduli_array):
    """
    Computes the gravitational intereactions between the particles of the system
                
    :param G: gravitational constant (float)
    :param p3d_list: list in which each item is a P3D instance
    :param separation:array: [N,N,3] array with the pair separations 
    :param moduli_array: [N,N] array with the pair distances between particles
        
    :return force_array:  array containing gravitational forces 
    """
    N = len(p3d_list)
    force_i_array = np.zeros((N, N, 3))    
    
    
    for i in range(N):
        for j in range(i+1,N):
            #by specifying i+1 in the range only half of the matrix is computed, as the whole matrix is symmetric 
            #calculate the force:
            force_i_array[i,j] = -G * p3d_list[i].mass * p3d_list[j].mass * separation_array[i, j]/(moduli_array[i, j])**3
            force_i_array[j,i] = -force_i_array[i, j]
            
    force_array = np.sum(force_i_array, axis = 1)

    return force_array


def potential_gravitational(G, p3d_list, moduli_array):
    """
    Computes the total potential energy (due to the gravitational interaction)
                
    :param G: gravitational constant (float)
    :param p3d_list: list in which each item is a P3D instance
    :param moduli_array: [N,N] array with the pair distances between particles
        
    :return potential_energy: total potential energy of the system
    """

    N = len(p3d_list)    
    #potential_array = np.zeros([N,N])
    potential_energy = 0.0
    
    for i in range(N):
        for j in range(i+1,N):
            #calculate the potential energy in each interaction 
            potential_energy -= p3d_list[i].mass*p3d_list[j].mass*G/moduli_array[i,j]
              
    return potential_energy



def orbital_period(sun_distances, p3d_list, time_list, observables, moon_earth, p3d_list_moon):
    """
    Computes the orbital periods for all bodies except the sun, and the Moon orbit
    around the Earth. The values obtained are written into a given output file.
                
    :param sun_distances: [numstep, N-2] array containig the distances between
    the Sun and all the other bodies (except the Moon)
    :param p3d_list: list in which each item is a P3D instance
    :param time_list: list contaninng all the time steps
    :param observables: outfile to write the period values in 
    :param moon_eart: [numstep, N] array containg the distances between the Moon and the Earth
    :p3d_list_moon: list in which each item is a P3D instance (without the moon)
        

    """
    #for every body except the Moon
    
    #create a list to store all the period values
    period_list = []
  
    #the range is with the -2 because the Sun and the Moon are not taken in account
    for i in range(len(p3d_list)-2):
        
        #max_values used to calculate the prominance
        max_value = np.zeros([len(p3d_list)-2])
        max_value[i] = sun_distances[:,i].max()

        #prominance taken in account to avoid small peaks in Neptune orbit
        peaks = signal.find_peaks(sun_distances[:,i], prominence = max_value[i]/100)[0]
   
        period = 0
        
        for j in range(len(peaks)-1):
            #add all the orbits done by each planet (taking the time from the time_list)
            period += time_list[peaks[j+1]] - time_list[peaks[j]]

        #in case that more than 1 orbit was completed by a planet the value will be averaged
        if len(peaks) == 0:
            period_list.append("Not enough time to complete one period, please try again with a larger numstep/dt.")
        if len(peaks) == 1:
            period_list.append(period)
        else:
            period = period/(len(peaks)-1)
            period_list.append(period)
    
    #write the period values obtained into the output file
    for i in range(1, len(p3d_list)-1):      
        observables.write("{0:s} Orbit Period: {1:s}\n".format(p3d_list_moon[i].label,
                          str(period_list[i-1])))
    
    #procedure to find the period of the Moon around the Earth
    
    #moon_earth only contains the distances between the Moon and the Earth
    moon_peak = signal.find_peaks(moon_earth[:,0])[0]
    moon_periods = 0
    
    for y in range(len(moon_peak)-1):
        moon_periods += time_list[moon_peak[y+1]] - time_list[moon_peak[y]]
    
    #average in case that mora than one orbit was completed
    if len(moon_peak) == 0:
        period_list.append("Not enough time to complete one period, please try again with a larger Numstep.")
    if len(moon_peak) == 1:
        period_list.append(moon_periods)
    else:
        moon_periods = moon_periods/(len(moon_peak)-1)
        period_list.append(moon_periods)
    
    #write the period value into the output file
    observables.write("Moon Orbit Period: {0:8f}\n".format(moon_periods))
    
    observables.write("\n")



def apsides(sun_distances, p3d_list, observables, moon_earth, p3d_list_moon):
    """
    Computes the orbit apsides for all bodies except the sun, and the Moon apsides
    for its orbit around the Earth. The values obtained are written into a given output file.
                
    :param sun_distances: [numstep, N-2] array containig the distances between
                          the Sun and all the other bodies (except the Moon)
    :param p3d_list: list in which each item is a P3D instance
    :param observables: outfile to write the period values in 
    :param moon_eart: [numstep, N] array containg the distances between the Moon and the Earth
    :p3d_list_moon: list in which each item is a P3D instance (without the moon)
        

    """
    
    #for every body except the Moon
    
    apside_max = [] #list to store the Aphelia
    apside_max = sun_distances.max(axis=0)
    apside_min = [] #list to store the Periphelia
    apside_min = sun_distances.min(axis=0)    
        
    
    #write the apside values obtained into the ourput file
    for i in range(len(p3d_list)-2):
        observables.write("{0:s} Aphelion: {1:8f}, Periphelion: {2:8f}\n".format(p3d_list_moon[i+1].label,
                          apside_max[i],
                          apside_min[i]))

    #for the Moon apsides
    
    moon_max = moon_earth.max()
    moon_min = moon_earth.min()
        
    observables.write("Moon Aphelion: {0:8f}, Periphelion: {1:8f}\n".format(moon_max, moon_min))




def main():
    
    """
    After reading the input files from the command line and making use of the
    Particle3D module, it carries out a Velocity Verlet Time Integration, calls 
    the methods to compute the period orbits an the apsides. Plots the energies 
    into graphs to help to compare them.

    """
  
    # Read name of output files and input files from command line
    #example: python3 simulation.py initials.txt tajectory.xyz observables.txt energy.txt parameters.txt
    
    if len(sys.argv)!=6:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file>")
        quit()
    else:
        outfile_name = sys.argv[2]  #to write the trajectories (XYZ-format)
        input_file = sys.argv[1]  #input of the initial conditions of the particles
        observables_file = sys.argv[3]  #to write the apsides and orbital
        energy_file = sys.argv[4]  #to crite the energy values
        parameters = sys.argv[5]  #input of the time integration parameters

    #read parameters for the time integration from input file
    parameters = open(parameters, 'r')
    for line in parameters.readlines():
        dt = float(line.split(",")[0])
        numstep = int(line.split(",")[1])
        G = float(line.split(",")[2])
        
    parameters.close()
        
    
    # Open output files
    outfile = open(outfile_name, "w")
    observables = open(observables_file, "w")
    energy_output = open(energy_file, "w")
    
    
    #create the list with the Particle3D instances 
    p3d_list = Particle3D.read_particle(input_file)
    N = len(p3d_list)

    #correct for the centre of mass
    total_mass, com_vel = Particle3D.com_velocity(p3d_list)
    for i in range(N):   
        p3d_list[i].velocity -= com_vel
    
    #compute separation and moduli arrays by using the Particle3D module
    separation_array, moduli_array = Particle3D.pair_separations(p3d_list)

    
    time = 0.0   #starting time of the time integration
    
    kinetic_energy = Particle3D.sys_kinetic(p3d_list)  #computes total kinetic energy of the system
    potential_energy = potential_gravitational(G, p3d_list, moduli_array)  #computes total potential energy of the system
    total_energy = kinetic_energy + potential_energy    #total energy of the system
    
    #write total energies values into the output file for energies (initial value)
    energy_output.write("{0:f} {1:8f} {2:8f} {3:8f}\n".format(time,
                        total_energy,
                        potential_energy,
                        kinetic_energy))
    
    #create necessary lists 
    time_list = [time]
    energy_list = [total_energy]
    potential_list = [potential_energy]
    kinetic_list = [kinetic_energy]  
          
    #create xero arrays used inside the time integration loop
    sun_distances = np.zeros([numstep, N-2])
    moon_earth = np.zeros([numstep,1])
    

    # Start the time integration loop
    for i in range(numstep):
        #start writing into the output file for the trajectories
        outfile.write("{0}\n".format(N))
        outfile.write("point={0}\n".format(i))
        
       
                    
           
        #calculate pair wise forces between bodies
        force_array = force_gravitational(G, p3d_list, separation_array, moduli_array)
        Particle3D.update_pos_2nd_list(p3d_list, dt, force_array) # Update particles position
        
        #update positions and distances arrays
        separation_array, moduli_array = Particle3D.pair_separations(p3d_list)
            
                
        # Update the pair-wise forces
        force_array_updated = force_gravitational(G, p3d_list, separation_array, moduli_array)
        force_total = 0.5*(force_array + force_array_updated)
        # Update particles velocities by averaging current and new forces
        Particle3D.update_vel_list(p3d_list, dt, force_total)
    
        # Re-define pair-wise forces value
        force_array = force_array_updated
             
        # Increase time
        time += dt
            
        #update energies
        potential_energy = potential_gravitational(G, p3d_list, moduli_array)
        kinetic_energy = Particle3D.sys_kinetic(p3d_list)
        total_energy = potential_energy + kinetic_energy
        
        # Output energy information
        energy_output.write("{0:f} {1:8f} {2:8f} {3:8f}\n".format(time,
                        total_energy,
                        potential_energy,
                        kinetic_energy))
        
        # Append information to data lists
        time_list.append(time)
        energy_list.append(total_energy)
        potential_list.append(potential_energy)
        kinetic_list.append(kinetic_energy)


        #takes the row of the separation between the sun and the other particles and takes off the sun-sun separation colum
        #fills moon-eart distances array
        for g in range(N):
            
            if p3d_list[g].label.lower() == "sun":
                sun_index = g
          
                for p in range(N): 
                    if p3d_list[p].label.lower() == "moon":
                        moon_index = p
                        
                        #create a list without the moon 
                        p3d_list_moon = p3d_list[:moon_index] + p3d_list[moon_index+1:]
                        
                        moduli_apsides = np.delete(moduli_array, moon_index, 1)                   
                        moduli_apsides = moduli_apsides[sun_index, sun_index+1:]
                        #sun_distances contains the distances between the Sun and the other bodies except the Moon
                        sun_distances[i,:] = moduli_apsides  
                 
                        for h in range(N):
                            if p3d_list[h].label.lower() == "earth":
                                earth_index = h
                                moduli_apsides = moduli_array[earth_index, moon_index]
                                #moon_earth contains the distances between the Moon and the Earth
                                moon_earth[i] = moduli_apsides

        for j in range(N):
            #write into outfile for trajectories, for each of the bodies in the simulation
            outfile.write("{0:s} {1:8f} {2:8f} {3:8f}\n".format(p3d_list[j].label,
                          p3d_list[j].position[0],
                          p3d_list[j].position[1],
                          p3d_list[j].position[2]))
        

    # Post-simulation:
    
    #call methods to compute the apsides and the orbit periods
    orbital_period(sun_distances, p3d_list, time_list, observables, moon_earth, p3d_list_moon)
    apsides(sun_distances, p3d_list, observables, moon_earth, p3d_list_moon)

    
    # Close files
    outfile.close()
    energy_output.close()
    observables.close()

    # Plot total energies to screen
    pyplot.title('Energy vs Time')
    pyplot.xlabel('Time (days)')
    pyplot.ylabel('Energy')
    pyplot.plot(time_list, energy_list, "r")
    pyplot.plot(time_list, potential_list, "b")
    pyplot.plot(time_list, kinetic_list, "g")
    pyplot.legend(("Total Energy", "Potential Energy", "Kinetic Energy"))
    pyplot.show()
    
    #plot total energy to screen
    pyplot.title('Total Energy vs Time')
    pyplot.xlabel('Time (days)')
    pyplot.ylabel('total Energy')
    pyplot.plot(time_list, energy_list)
    pyplot.show()     

    #plot potential energy to screen
    pyplot.title('Potential Energy vs Time')
    pyplot.xlabel('Time (days)')
    pyplot.ylabel('Potential Energy')
    pyplot.plot(time_list, potential_list)
    pyplot.show()    
    
    #plot potential energy to screen
    pyplot.title('Kinetic Energy vs Time')
    pyplot.xlabel('Time (days)')
    pyplot.ylabel(' Kinetic Energy')
    pyplot.plot(time_list, kinetic_list)
    pyplot.show()  
    
    
# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
    

    
    
