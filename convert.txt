# read in a parXXX.vtu file and return a text file with particle positions only
# points.dat format: 
# Number of points
# Material | x location | y location | z location | diameter
import numpy as np
import pyvista as pv
import fortranformat as ff

fin = 'ConvertToPoints.vtu'
fout = 'points.dat'

# Additional column data entries
print('Enter Particle Diameter: ')
particle_diameter = float(input())
print('Enter Material of Particle: ')
particle_material = input()

# Set for barrel case when barrel length is 40mm and particle bed length is 6mm
x_offset = 0.0
y_offset = 0.0
z_offset = 0.034

data = pv.read('ConvertToPoints.vtu')
points = data.points

# sort in x direction
points = points[points[:, 0].argsort()]

# add x offset
points[:,0] = points[:,0] + x_offset

# add y offset
points[:,1] = points[:,1] + y_offset

# add z offset
points[:,2] = points[:,2] + z_offset


writer = ff.FortranRecordWriter('(A16, 4E23.16)')
writer_head = ff.FortranRecordWriter('(I16)')

n = data.n_points
s0 = writer_head.write([n]) + '\n'
s1 = ''
s2 = ''
s3 = ''
s4 = ''
s5 = ''
s6 = ''
s7 = ''
s8 = ''
s9 = ''
s10 = ''
s11 = ''
print('Number of Particles: ',n-1)
print('x,y,z offsets are hardcoded')
print('x-offset: ',x_offset,'(m)')
print('y-offset: ',y_offset,'(m)')
print('z-offset: ',z_offset,'(m)')
for j in range(0,11):
    a = int(n/10)*j
    b = int(n/10)*(j+1)
    if j == 10:
        b = n
    if j == 0:
        print('Precent complete: 0% (Updates every 10%)', end='\r')
    for i in range(a,b):   
        if j == 0:
            s1 += writer.write([particle_material, points[i,0], points[i,1], points[i,2], particle_diameter]) + '\n'
        elif j == 1:
            s2 += writer.write([particle_material, points[i,0], points[i,1], points[i,2], particle_diameter]) + '\n'    
        elif j == 2:
            s3 += writer.write([particle_material, points[i,0], points[i,1], points[i,2], particle_diameter]) + '\n'
        elif j == 3:
            s4 += writer.write([particle_material, points[i,0], points[i,1], points[i,2], particle_diameter]) + '\n'
        elif j == 4:
            s5 += writer.write([particle_material, points[i,0], points[i,1], points[i,2], particle_diameter]) + '\n'
        elif j == 5:
            s6 += writer.write([particle_material, points[i,0], points[i,1], points[i,2], particle_diameter]) + '\n'
        elif j == 6:
            s7 += writer.write([particle_material, points[i,0], points[i,1], points[i,2], particle_diameter]) + '\n'
        elif j == 7:
            s8 += writer.write([particle_material, points[i,0], points[i,1], points[i,2], particle_diameter]) + '\n'
        elif j == 8:
            s9 += writer.write([particle_material, points[i,0], points[i,1], points[i,2], particle_diameter]) + '\n'
        elif j == 9:
            s10 += writer.write([particle_material, points[i,0], points[i,1], points[i,2], particle_diameter]) + '\n'
        elif j == 10:
            s11 += writer.write([particle_material, points[i,0], points[i,1], points[i,2], particle_diameter]) + '\n'
    if j <9:
        print('Precent complete: ', (j+1)*10,'%                      ', end='\r')
final = s0+s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11
with open(fout, 'w') as f:
    f.write(final)
print('Done                                                             ')
