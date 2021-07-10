import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import glob

def photometry(number_image, filter_image, number_star, GAIN, graphing):
	#downloading raw image
	name_image = 'IR02143-1MHz-76mcs-PreampX4-' + str(number_image).zfill(4) + str(filter_image) + '.fit'
	file_image = fits.open('data/2021-11-10/' + name_image)
	altitude = float(file_image[0].header['OBJCTALT']) #altitude
	exposure_time = float(file_image[0].header['EXPTIME']) #exposure time
	file_image.close()
	
	#downloading clean image
	name_image = 'clean' + str(number_image).zfill(4) + str(filter_image) + '.fits'
	file_image = fits.open('data/clean_images/' + name_image)
	data_image = file_image[0].data #2D array
	dim1, dim2 = len(data_image), len(data_image[0,:]) 
	
	ADU = 0
	#downloading coordinates
	name_coordinates = 'coordinates.txt'
	file_coordinates = open(name_coordinates)
	cnt= 0
	for line in file_coordinates:
		dots = [int(number) for number in line.strip().split(', ')]
		cnt+=1
		if(cnt == int(number_image)):
			if (number_star == 'st2'):
				coordinates = dots[0], dots[1]
			if (number_star == 'st3'):
				coordinates = dots[2], dots[3]
			break
	
	#calculating ADU
	frame_size = 15 #photutils???
	x_center, y_center = coordinates
	star = data_image[y_center - frame_size : y_center + frame_size, x_center - frame_size: x_center + frame_size]
	sky = np.median(star) #substracting sky background
	
	#image of a star
	if (graphing == 1):
		plt.imshow(star - sky, cmap='gray')
		plt.colorbar()
		plt.title('star')
		plt.show()
		
	ADU = np.sum(star - sky) / exposure_time #exposure varies with image number?
	z = 90.0 - altitude
	result = ADU * GAIN, z
	return result

def links_filter(filter_image):
	names = glob.glob('data/clean_images/*' + filter_image + '.fits')
	numbers = [int(name[23:27]) for name in names]
	return numbers

#testing	
#print(photometry(number_image = 1, filter_image = 'B', number_star = 'st3', GAIN = 1, graphing = 1))
#print(links_filter('B'))

filters = ['U', 'B', 'V', 'Rc', 'Ic']
stars = ['st2', 'st3']
GAIN = 1.0

for filter_test in filters:
	for star_test in stars:
		numbers_image = links_filter(filter_test)
		N = len(numbers_image)
		cnt = 0
		z = []
		electron_count = []
		
		#getting z, e-
		for number_image in numbers_image:
			electron_count_single, z_single = photometry(number_image, filter_test, star_test, GAIN, graphing = 0)
			electron_count.append(electron_count_single)
			z.append(z_single * np.pi / 180.0)
			cnt += 1
			#print(str(cnt) + ' / ' + str(N) + ' done')
	
		#saving graphs
		plt.plot(1/np.cos(np.array(z)), np.log(np.array(electron_count)), 'ko')
		plt.title('object ' + star_test + ' in filter ' + filter_test + ' (GAIN = 1)')
		plt.xlabel(r'$\sec(z)$')
		plt.ylabel(r'$\ln (e^{-} / t_{exp})$')
		plt.savefig('data/extinction_graphs/' + star_test + '_' + filter_test + '.png')
		print('done ' + filter_test + ' ' + star_test)
		plt.clf()


