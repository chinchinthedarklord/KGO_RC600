import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from lmfit import Model

def linear(x, A):
	result =  A * x
	return result

linear_model = Model(linear)

dim1 = dim2 = 2048
filter_image = 'Ic'

#BIAS
#print('Enter number of biases:')
#number_flats = int(input())
#------------------------------------------ 
number_bias = 32
#------------------------------------------ 
data_bias_mean = np.zeros((dim1, dim2))
for i in range(number_bias):
	name_bias = 'BIAS-1MHz-76mcs-PreampX4-' + str(i+1).zfill(4) + '.fit'
	file_bias = fits.open('data/bias/' + name_bias)
	data_bias = file_bias[0].data #2D array
	data_bias = data_bias.astype('float64')
	file_bias.close()
	data_bias_mean += data_bias
data_bias_mean = data_bias_mean / float(number_bias)

#FLAT
#print('Enter number of flats:')
#number_flats = int(input())
#------------------------------------------ 
number_flats = 12 
#------------------------------------------ 
signal = []
noise = []
for i in range(number_flats):
	name_flat = 'FLAT-1MHz-76mcs-PreampX4-'+ str(i+1).zfill(4) + filter_image + '.fit'
	file_flat = fits.open('data/flat/' + name_flat)
	data_flat = file_flat[0].data - data_bias_mean #2D array
	average_value = np.sum(data_flat[dim1/ 4: dim1 * 3 / 4, dim2/ 4: dim2 * 3 / 4]) * 4 / float(dim1 * dim2)
	noise_single = np.var(data_flat[dim1/ 4: dim1 * 3 / 4, dim2/ 4: dim2 * 3 / 4], dtype=np.float64)
	signal.append(average_value)
	noise.append(noise_single)

#fitting to find GAIN
result = linear_model.fit(noise, x = signal, A = 1.0)
GAIN = result.params.get('A').value
dGAIN = result.params.get('A').stderr
print('GAIN for filter ' + filter_image + ' = ' + str(GAIN) + ' +/- ' + str(dGAIN))

#plotting	
plt.plot(signal, noise, 'ko', label = 'data')
plt.plot(signal, result.best_fit, 'r', label = 'best fit')
plt.title('GAIN for filter ' + filter_image)
plt.xlabel('signal (ADU)')
plt.ylabel(r'$\sigma^2$' + ' (ADU)')
plt.grid()
plt.legend()
plt.show()