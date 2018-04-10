import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from astropy.io import fits

#define the kernels
sobel_x = np.array([[1,0,-1],[2,0,-2],[1,0,-1]])
sobel_y = np.array([[1,2,1],[0,0,0],[-1,-2,-1]])


#load the test image
def rawFits(filename):
    f = fits.open(filename)
    h = f[0].header
    d = f[0].data
    f.close()
    return h, d

def image_stack(im_list,verbose=False):
	#function to stack a list of images, producing a higher signal to noise ratio
	total = np.copy(im_list[0])
	try:
		for i in im_list[1:]:
			total += i
	except:
		print('ERROR: You did not provide enough images to stack!')
	if verbose:
		fig = plt.figure()
		plt.imshow(total, cmap='gray',origin='lower',interpolate='none')
		plt.show()
		plt.close()
	return total

def brightness(im, verbose=False):
	#function to measure the mean log intensity of an image for cloud cover measurements
	NIGHT = 0.
	DAY = 1.
	brightness = (np.log10(np.nanmean(im)) - NIGHT) / DAY
	return 


def sobel(im,verbose=False):
	#----------sobel edge detector--------------
	#convolve the image with each kernel
	gx = signal.convolve2d(im, sobel_x,mode='same',boundary='symm')
	gy = signal.convolve2d(im, sobel_y,mode='same',boundary='symm')

	#compute the magnitude G = sqrt(Gx^2+Gy^2) and direction Theta = atan(Gy/Gx) of the edges
	G = np.nan_to_num(np.sqrt(np.power(gx,2.) + np.power(gy,2.)))
	Theta = np.arctan2(gy,gx)

	#keep a threshold value of edges???
	new = np.copy(G)
	new[G > 10] = 10
	num_edges = len(new[new > 0.5])
	print(num_edges/(640*480.) * 100, '% Edges')

	if verbose:
		#show the detected edges
		fig,ax = plt.subplots(2,2)
		ax[0,0].imshow(np.rot90(im),origin='lower',cmap='gray')
		ax[0,0].set_title('Original Image')
		ax[1,0].imshow(np.rot90(new),origin='lower')
		ax[1,0].set_title('Edges (Gradient Magnitude)')
		ax[1,1].imshow(np.rot90(Theta),origin='lower')
		ax[1,1].set_title('Edges (Gradient Direction)')
		ax[0,1].quiver(gx[::5,::5],gy[::5,::5])
		ax[0,0].axis('off')
		ax[1,0].axis('off')
		ax[1,1].axis('off')
		ax[0,1].axis('off')
		plt.show()

if __name__ == '__main__':
	h, TEST = rawFits('./assets/tesCurrentImage.FIT')
	h,TEST = rawFits('./assets/samenight1.FIT')
	im = -1*np.log(TEST)
	sobel(im,verbose=True)