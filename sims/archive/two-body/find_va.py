import numpy as np

va_array = np.loadtxt("outparray.csv", delimiter=',')
mag = 0
for i in range(len(va_array[:,0])):
    line_m = np.sqrt(va_array[i,3]**2 + va_array[i,4]**2)
    if line_m > mag:
        mag = line_m
        x_coord = va_array[i,1]
        y_coord = va_array[i,2]
        ra = np.sqrt(x_coord**2 + y_coord**2)

print("A maximum velocity of "+str(mag)+" was reached at a distance "+str(ra)+" with components "+str(x_coord)+" and "+str(y_coord))



