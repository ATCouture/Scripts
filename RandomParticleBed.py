import random as rd
import numpy as np
import fortranformat as ff
import math
import matplotlib.pyplot as plt

# Version 1.0
# Code requests parameters from user and outputs a points.dat file describing a randomly distributed
# particle bed in the following format:
# Number_of_Particles
# Particle Material | Particle Center x-coordinate | Particle Center y-coordinate | Particle Center z-coordinate | Particle Diameter
# Code currently supports 3 Dimensional Cartesian, Cylindrical, and Spherical Coordinate Input and Outputs in Cartesian Coordinates

fout = 'points.dat'
writer = ff.FortranRecordWriter('(A16, 4E23.16)')
writer_head = ff.FortranRecordWriter('(I16)')
error = 0

print('Enter desired particle volume fraction: ')
particles_vol_fract = float(input())
print('Particle Volume Fraction = ',particles_vol_fract*100,'%')

print('Enter particle diameter:')
particles_dia = float(input())
print('Particles diameter: ',particles_dia)
particles_vol = math.pi/6*particles_dia**3

print('Enter particle material:')
part_material = input()
print('Particle Material: ',part_material)

print('Type coordinates (cart=1, cyld=2, sph=3)')
coord = int(input())

if coord == 1:
    print('Cartesian Coordinate Frame (3-dimensional)')
    print('Enter minimum x-coordinate boundary:')
    x_min = float(input()) + particles_dia/2
    print('Enter maximum x-coordinate boundary:')
    x_max = float(input()) - particles_dia/2
    print('Enter minimum y-coordinate boundary:')
    y_min = float(input()) + particles_dia/2
    print('Enter maximum y-coordinate boundary:')
    y_max = float(input()) - particles_dia/2
    print('Enter minimum z-coordinate boundary:')
    z_min = float(input()) + particles_dia/2
    print('Enter maximum z-coordinate boundary:')
    z_max = float(input()) - particles_dia/2

    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min

    total_vol = (x_range + particles_dia)*(y_range + particles_dia)*(z_range + particles_dia)

    particles_num = int(particles_vol_fract*total_vol/particles_vol)
    print('Number of particles: '+ str(particles_num))
    part_coord = np.zeros((particles_num,3))
    i = 0
    while i < particles_num:
        # Initial partical location
        part_coord[i,0] = rd.random()*x_range+x_min
        part_coord[i,1] = rd.random()*y_range+y_min
        part_coord[i,2] = rd.random()*z_range+z_min
        overlap = 0
        if i > 0:
            for j in range(i): #loops through all nonzero positioned particles
                delta = ((part_coord[i,0] - part_coord[j,0])**2 + (part_coord[i,1] - part_coord[j,1])**2 + (part_coord[i,2] - part_coord[j,2])**2)**(1/2)
                if delta <= particles_dia and i != j: # When true, particle is overlapped with another existing particle
                    overlap = 1 # This skips i = i + 1 and repeats while loop for this particle
                    break # This ensures that code doesn't look for more overlaps after initial overlap is found
        if overlap == 0: # This means that the new particle is in a unique location
            i = i + 1
        percent_done = i/particles_num*100
        print('Percent of particles placed: ' + str(percent_done), end='\r')
elif coord == 2:
    print('Cylindrical Coordinate Frame')

    print('Enter minimum radial-coordinate boundary:')
    rad_min = float(input()) + particles_dia/2
    print('Enter maximum radial-coordinate boundary:')
    rad_max = float(input()) - particles_dia/2
    print('Assuing the azimuth is full 360 degrees.')
    print('Enter minimum z-coordinate boundary:')
    z_min = float(input()) + particles_dia/2
    print('Enter maximum z-coordinate boundary:')
    z_max = float(input()) - particles_dia/2
    rad_range = rad_max - rad_min
    az_range = 2*math.pi
    z_range = z_max - z_min
    total_vol = math.pi*((rad_max + particles_dia/2)**2 - (rad_min - particles_dia/2)**2)*(z_range + particles_dia)

    particles_num = int(particles_vol_fract*total_vol/particles_vol)
    print('Number of particles: '+ str(particles_num))
    part_coord = np.zeros((particles_num,3))
    i = 0
    while i < particles_num:
        # Initial partical location
        az_angle = az_range*rd.random()
        part_coord[i,0] = (rd.random()*rad_range + rad_min)*math.cos(az_angle)
        part_coord[i,1] = (rd.random()*rad_range + rad_min)*math.sin(az_angle)
        part_coord[i,2] = rd.random()*z_range+z_min
        overlap = 0
        if i > 0:
            for j in range(i): #loops through all nonzero positioned particles
                delta = ((part_coord[i,0] - part_coord[j,0])**2 + (part_coord[i,1] - part_coord[j,1])**2 + (part_coord[i,2] - part_coord[j,2])**2)**(1/2)
                if delta <= particles_dia and i != j: # When true, particle is overlapped with another existing particle
                    overlap = 1 # This skips i = i + 1 and repeats while loop for this particle
                    break # This ensures that code doesn't look for more overlaps after initial overlap is found
        if overlap == 0: # This means that the new particle is in a unique location
            i = i + 1
        percent_done = i/particles_num*100
        print('Percent of particles placed: ' + str(percent_done), end='\r')

elif coord == 3:
    print('Spherical Coordinate Frame')
    print('Enter minimum radial-coordinate boundary:')
    rad_min = float(input()) + particles_dia/2
    print('Enter maximum radial-coordinate boundary:')
    rad_max = float(input()) - particles_dia/2
    print('Assuing the azimuth and polar angles are full 360 degrees.')

    rad_range = rad_max - rad_min
    az_range = 2*math.pi
    psi_range = 2*math.pi
    total_vol = 4/3*math.pi*((rad_max + particles_dia/2)**3 - (rad_min - particles_dia/2)**3)

    particles_num = int(particles_vol_fract*total_vol/particles_vol)
    print('Number of particles: '+ str(particles_num))
    part_coord = np.zeros((particles_num,3))
    i = 0
    while i < particles_num:
        # Initial partical location
        az_angle = az_range*rd.random()
        psi_angle = az_range*rd.random()
        part_coord[i,0] = (rd.random()*rad_range + rad_min)*math.cos(rd.random()*az_range)*math.sin(rd.random()*psi_range)
        part_coord[i,1] = (rd.random()*rad_range + rad_min)*math.sin(rd.random()*az_range)*math.sin(rd.random()*psi_range)
        part_coord[i,2] = (rd.random()*rad_range + rad_min)*math.cos(rd.random()*psi_range)
        overlap = 0
        if i > 0:
            for j in range(i): #loops through all nonzero positioned particles
                delta = ((part_coord[i,0] - part_coord[j,0])**2 + (part_coord[i,1] - part_coord[j,1])**2 + (part_coord[i,2] - part_coord[j,2])**2)**(1/2)
                if delta <= particles_dia and i != j: # When true, particle is overlapped with another existing particle
                    overlap = 1 # This skips i = i + 1 and repeats while loop for this particle
                    break # This ensures that code doesn't look for more overlaps after initial overlap is found
        if overlap == 0: # This means that the new particle is in a unique location
            i = i + 1
        percent_done = i/particles_num*100
        print('Percent of particles placed: ' + str(percent_done), end='\r')
else:
    print('Error: incorect coordinate input')
    error == 1


size = np.ones(len(part_coord[:,1]))*2
plt.scatter(part_coord[:,0],part_coord[:,1], s=size)
plt.show()


if error != 1:
    with open(fout,'w') as f:
        f.write(writer_head.write([particles_num])+'\n')
        for k in range(particles_num):
            f.write(writer.write([part_material,part_coord[k,0],part_coord[k,1],part_coord[k,2],particles_dia]) + '\n')
print('Done!                                                        ')
