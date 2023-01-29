
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 19:24:23 2020

Particle3D, a class to describe point particles in 3D space

An instance describes a particle in Euclidean 3D space: 
velocity and position are [3] arrays

Includes time integrator methods (for a single particle and for a system,
 contained in a list) + methods to calculate the kinetic energy and 
momentum of a particle + methods to calculate the total kinetic energy, total 
mass and Center of Mass velocity of a system of particles contained in a list

@author: Marta-Villa de Toro Sánchez s1867522

"""
import math
import numpy as np


class Particle3D(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label: name of the particle
    mass: mass of the particle
    pos: position of the particle
    vel: velocity of the particle

        Methods:
    __init__
    __str__
    kinetic_e  - computes the kinetic energy
    momentum - computes the linear momentum
    update_pos_1st - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel - updates the velocity

        Static Methods:
    new_particle - initializes a P3D instance from a file handle
    sys_kinetic - computes total K.E. of a p3d list
    com_velocity - computes total mass and CoM velocity of a p3d list
    """

    def __init__(self,label, mass, position, velocity):
        """
        Initialises a particle in 3D space

        :param label: String w/ the name of the particle
        :param mass: float, mass of the particle
        :param position: [3] float array w/ position
        :param velocity: [3] float array w/ velocity
        """
        #the variables to define each self come from the input, handled in new_particle static method
        self.label = label        #string 
        self.mass = mass            #float
        self.position = position    #array of 3 floats
        self.velocity = velocity   #array of 3 floats
        self.label = label
        


    @staticmethod
    def read_particle(inputfile):
        """
        Initialises a Particle3D instance given an input file handle.
        
        The input file should contain one per planet in the following format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>
        
        :param inputfile: Readable file handle in the above format

        :return Particle3D instance
        """
        p3d_list = [] #create a list to store all the particles
        
        inputfile = open(inputfile, "r")   #input file from main module  
        
        for line in inputfile.readlines(): 
        #read lines of the input file,each line is a particle
            
            label = line.split(",")[0]
            mass = float(line.split(",")[1])
            position = np.array(line.split(",")[2:5], float)
            velocity = np.array(line.split(",")[5:8], float)
            particle = Particle3D(label, mass, position, velocity)

            p3d_list.append(particle)

        inputfile.close() 
      
        #to acces to a singular particle use p3d_list(i)    
        
        return p3d_list #call Particle3D __innit__ method

    def __str__(self):
        """
        XYZ-compliant string. The format is
        <label>    <x>  <y>  <z>
        """        
        
        return str(self.label) +" "+ str(self.position[0]) +" "+ str(self.position[1]) +" "+ str(self.position[2])


    @staticmethod
    def pair_separations(p3d_list):
        """
        Creates arrays with the pair separations between the particles and the pair distances
        
        :param p3d_list: list in which each item is a P3D instance
        
        :returns separation_array: [N,N,3] array with the pair separations 
        :returns moduli_array: [N,N] array with the pair distances between particles
        
        """
    
        N = len(p3d_list)
        separation_array = np.zeros([N, N, 3])
        moduli_array = np.zeros([N, N])
        
        
        for i in range(N):       
            for j in range(i+1,N):
            #by specifying i+1 in the range only half of the matrix is computed, as the whole matrix is symmetric 
                separation_array[i,j] = p3d_list[i].position - p3d_list[j].position
                moduli_array[i,j] = np.linalg.norm(separation_array[i,j])

        return separation_array, moduli_array
        
    def kinetic_e(self):
        """
        Returns the kinetic energy of a Particle3D instance

        :return ke: float, 1/2 m v**2
        """
    
        return 0.5*(self.mass*(np.linalg.norm(self.velocity)**2))


    def momentum(self):
        """
        Returns the linear momentum of a Particle3D instance
        :return p: (3) float np.array, m*v
        """
        
        return self.mass*self.velocity


    def update_pos_1st(self, dt):
        """
        1st order position update

        :param dt: timestep
        """
        self.position += dt*self.velocity
        


    def update_pos_2nd(self, dt, force):
        """
        2nd order position update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        """
        self.position += dt*self.velocity + dt**2*force/(2*self.mass)
        


    def update_vel(self, dt, force):
        """
        Velocity update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        """
        self.velocity += dt*force/self.mass
        
     
    @staticmethod   
    def update_pos_2nd_list(p3d_list, dt, force_array):
        """
        Computes the updated position (2nd) for the whole system
        
        :param p3d_list: list in which each item is a P3D instance
        :param dt: timestep
        :param force_array: [N,N,3] array with the pair wise forces between all particles
        
        """
        
        N = len(p3d_list)
        for i in range(N):    
            #making use of the method update_pos_2nd
            p3d_list[i].update_pos_2nd(dt, force_array[i])
            
    @staticmethod        
    def update_vel_list(p3d_list, dt, force_array):
        """
        Computes the updated velocity for the whole system
        
        :param p3d_list: list in which each item is a P3D instance
        :param dt: timestep
        :param force_array: [N,N,3] array with the pair wise forces between all particles
        
        """
    
        N = len(p3d_list)
        for i in range(N): 
            #making use of the method update_vel
            p3d_list[i].update_vel(dt, force_array[i])
        


    @staticmethod
    def sys_kinetic(p3d_list):
        """
        Returns the kinetic energy of the whole system

        :param p3d_list: list in which each item is a P3D instance

        :return sys_ke: \sum 1/2 m_i v_i^2 
        """

        sys_ke = sum([particle.kinetic_e() for particle in p3d_list])
        
        return sys_ke
    


    @staticmethod
    def com_velocity( p3d_list ):
        """
        Computes the total mass and CoM velocity of a list of P3D's

        :param p3d_list: list in which each item is a P3D instance
        :return total_mass: The total mass of the system 
        :return com_vel: Centre-of-mass velocity
        """
        total_momentum = sum([particle.momentum() for particle in p3d_list])
    
        total_mass = sum([particle.mass for particle in p3d_list])
        
        com_vel = total_momentum/total_mass
    
        return total_mass, com_vel
