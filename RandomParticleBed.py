import random as rd
import numpy as np
import fortranformat as ff
import math
fout = 'points.dat'
writer = ff.FortranRecordWriter('(A16, 4E23.16)')
writer_head = ff.FortranRecordWriter('(I16)')

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

    total_vol = x_range*y_range*z_range
    #x_range = total_vol**(1/3) # for equivalent volume of barrel
    #y_range = total_vol**(1/3) # for equivalent volume of barrel
    #z_range = total_vol**(1/3) # for equivalent volume of barrel

elif coord == 2:
    print('Cylindrical Coordinate Frame')

    print('Enter minimum radial-coordinate boundary:')
    rad_min = float(input()) + particles_dia/2
    print('Enter maximum radial-coordinate boundary:')
    rad_max = float(input()) - particles_dia/2
    print('Enter minimum azimuth angle in degrees:')
    az_min = input()

elif coord == 3:
    print('Spherical Coordinate Frame')

else:
    print('Error: incorect coordinate input')

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

final_string = writer_head.write([particles_num]) + '\n'
#for k in range(particles_num):
#    final_string += writer.write([part_material,part_coord[k,0],part_coord[k,1],part_coord[k,2],particles_dia]) + '\n'
with open(fout,'w') as f:
    f.write(writer_head.write([particles_num])+'\n')
    for k in range(particles_num):
        f.write(writer.write([part_material,part_coord[k,0],part_coord[k,1],part_coord[k,2],particles_dia]) + '\n')
print('Done!                                                        ')
