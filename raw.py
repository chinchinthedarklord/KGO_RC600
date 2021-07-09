import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os

def cleaning(number_image, filter_image, graphing, saving):
	#downloading image
	name_image = 'IR02143-1MHz-76mcs-PreampX4-' + str(number_image).zfill(4) + str(filter_image) + '.fit'
	file_image = fits.open('data/2021-11-10/' + name_image)
	data_image = file_image[0].data #2D array
	data_image = data_image.astype('float64')
	dim1, dim2 = len(data_image), len(data_image[0,:]) 
	exposure_time = file_image[0].header['EXPTIME'] #exposure time
	altitude = file_image[0].header['OBJCTALT'] #altitude
	file_image.close()

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
		data_image = data_image.astype('float64')
		file_bias.close()
		data_bias_mean += data_bias
	data_bias_mean = data_bias_mean / float(number_bias)	

	#DARK
	name_dark = 'DARK-1MHz-76mcs-PreampX4-' + str(number_image).zfill(4) + '-T010.fit'
	file_dark = fits.open('data/dark/' + name_dark)
	data_dark = file_dark[0].data #2D array
	data_dark = data_dark.astype('float64')
	dark_time = file_dark[0].header['EXPTIME'] #exposure time for dark
	file_dark.close()

	#FLAT
	#print('Enter number of flats:')
	#number_flats = int(input())
	#------------------------------------------ 
	number_flats = 12 
	#------------------------------------------ 
	data_flat_mean = np.zeros((dim1, dim2))
	for i in range(number_flats):
		name_flat = 'FLAT-1MHz-76mcs-PreampX4-'+ str(i+1).zfill(4) + filter_image + '.fit'
		file_flat = fits.open('data/flat/' + name_flat)
		data_flat = file_flat[0].data #2D array
		average_value = np.sum(data_flat[dim1/ 4: dim1 * 3 / 4, dim2/ 4: dim2 * 3 / 4]) * 4 / float(dim1 * dim2)
		data_flat = data_flat / average_value
		data_flat_mean += data_flat		
	data_flat_mean = data_flat_mean / float(number_flats)	

	#getting real image
	data_real = np.zeros((dim1, dim2))
	data_real = ((data_image - data_bias_mean) - (data_dark - data_bias_mean) * exposure_time / dark_time) / data_flat_mean

	#showing image
	if(graphing == 1):
		vmin, vmax = np.percentile(data_real, (5, 99.5))
		plt.imshow(data_real, cmap='gray', vmin = vmin, vmax = vmax)
		plt.colorbar()
		plt.title('real_image')
		plt.show()
	
	#saving fits
	if(saving == 1):
		hdu = fits.PrimaryHDU(data_real)
		hdul = fits.HDUList([hdu])
		hdul.writeto('data/clean_images/clean' + str(number_image).zfill(4) + str(filter_image) + '.fits')
	return 0

#cleaning(number_image, filter_image, graphing = 1, saving = 1)
for filename in os.listdir('data/2021-11-10'):
	number_index1 = filename.find('00')
	number_index2 = filename.find('.fit')
	number_image = int(filename[number_index1+2 : number_index1+4])
	filter_image = filename[number_index1+4: number_index2]
	print(number_image, filter_image)
	cleaning(number_image, filter_image, graphing = 0, saving = 1)