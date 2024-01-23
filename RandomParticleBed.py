import random as rd
import numpy as np
import fortranformat as ff
#fout = 'points.dat'
#writer = ff.FortranRecordWriter('(A16, 4E23.16)')
#writer_head = ff.FortranRecordWriter('(I16)')

# print('Enter desired particle volume fraction: ')
particles_vol_fract = 0.0001 # float(input())
print(particles_vol_fract)

# print('Enter particle diameter?')
particles_dia = 0.000025 # float(input())
print('Particles diameter: ' + str(particles_dia))
particles_vol = 3.14/6*particles_dia**3

# print('Type coordinates (cart=1, cyld=2, sph=3)')
coord = 1 #int(input())

if coord == 1:
    print('Cartesian Coordinate Frame')
    x_min = 0 #float(input())
    x_max = .001 #float(input())
    y_min = 0 #float(input())
    y_max = .001 #float(input())
    z_min = 0 #float(input())
    z_max = .001 #float(input())
    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min

    total_vol = 3.14*0.006**3 #x_range*y_range*z_range
    x_range = total_vol**(1/3) # for equivalent volume of barrel
    y_range = total_vol**(1/3) # for equivalent volume of barrel
    z_range = total_vol**(1/3) # for equivalent volume of barrel

elif coord == 2:
    print('Cylindrical Coordinate Frame')

elif coord == 3:
    print('Spherical Coordinate Frame')

else:
    print('Error: incorect coordinate input')

# print('Enter number of particles?')
particles_num = int(particles_vol_fract*total_vol/particles_vol) # 1000 # int(input()) Change to input later. maybe vol fraciton?
print('Number of particles: '+ str(particles_num))

part_x = np.zeros(particles_num)
part_y = np.zeros(particles_num)
part_z = np.zeros(particles_num)
count1 = 0
i = 0
while i < particles_num:
    # Initial partical location
    # insert function for random
    part_x[i] = rd.random()*x_range+x_min
    part_y[i] = rd.random()*y_range+y_min
    part_z[i] = rd.random()*z_range+z_min
    overlap = 0
    if i > 0:
        for j in range(i):
            delta = ((part_x[i] - part_x[j])**2 + (part_y[i] - part_y[j])**2 + (part_z[i] - part_z[j])**2)**(1/2)
            if delta <= particles_dia and i != j:
                overlap = 1
                count1 = count1 + 1
                break
    if overlap == 0:
        i = i + 1
    percent_done = i/particles_num*100
    print(str(percent_done), end='\r')

#count2 = 0
#for k in range(particles_num):
#    if k != 1:
#        for l in range(particles_num):
#            delta = ((part_x[k] - part_x[l])**2 + (part_y[k] - part_y[l])**2 + (part_z[k] - part_z[l])**2)**(1/2)
#            if delta < particles_dia + particles_dia/1000 and k != l:
#                count2 = count2 + 1
#                print('Overlap at: ',k,l)

print('Number of duplicates:',count1)
print('Number of particles:', particles_num)
print('Particle volume fraction percent:', particles_vol_fract*100)
print('Done!')
