sat_h = nc.variables['goes_imager_projection'].perspective_point_height
sat_lon = nc.variables['goes_imager_projection'].longitude_of_projection_origin
sat_sweep = nc.variables['goes_imager_projection'].sweep_angle_axis
X = nc.variables['x'][:] * sat_h
Y = nc.variables['y'][:] * sat_h
XX, YY = np.meshgrid(X, Y)
fk1 = f.variables['planck_fk1'][0]
fk2 = f.variables['planck_fk2'][0]
bc1 = f.variables['planck_bc1'][0]
bc2 = f.variables['planck_bc2'][0]

print("fk1 = ",fk1)
print("fk2 = ",fk2)
print("bc1 = ",bc1)
print("bc2 = ",bc2)
print(" ")


# Here's the scale factor and offset that are referred to - 
# I am not using them at all

scalefactor = f.variables['Rad'].scale_factor
add_offset = f.variables['Rad'].add_offset
print("scale factor is ",f.variables['Rad'].scale_factor)
print("add_offset   is ",f.variables['Rad'].add_offset)

# Reaed in Radiance values into data_var

data_var = f.variables['Rad'][:]

#Compute brightness temp from radiances
data_var = (fk2 / ( np.log((fk1 / data_var) + 1 )) - bc1) / bc2
