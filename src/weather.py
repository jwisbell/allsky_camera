import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from astropy.io import fits
from matplotlib import gridspec

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
		print 'ERROR: You did not provide enough images to stack!'
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
	brightness = (np.log10(np.nanmean(-1*im)) - NIGHT) / DAY
	return brightness


def sobel(im,verbose=False,save=False):
	#----------sobel edge detector--------------
	#convolve the image with each kernel
	gx = signal.convolve2d(im, sobel_x,mode='same',boundary='symm')
	gy = signal.convolve2d(im, sobel_y,mode='same',boundary='symm')

	#compute the magnitude G = sqrt(Gx^2+Gy^2) and direction Theta = atan(Gy/Gx) of the edges
	G = np.nan_to_num(np.sqrt(np.power(gx,2.) + np.power(gy,2.)))
	Theta = np.arctan2(gy,gx)

	#keep a threshold value of edges???
	new = np.copy(G)
	new[G > 0.5] = 10
	num_edges = len(new[new > 0.5])
	print num_edges/(640*480.) * 100, '% Edges'

	if verbose:
		#show the detected edges
		fig=plt.figure(figsize=(15,3))#,ax = plt.subplots(1,3,sharey=True)
		gs = gridspec.GridSpec(1, 3) #height_ratios=[3,1]) 
		ax = [None, None, None]
		ax[0] = plt.subplot(gs[0])
		ax[1] = plt.subplot(gs[1])
		ax[2] = plt.subplot(gs[2])
		ax[0].imshow(np.rot90(im),origin='lower',cmap='gray')
		ax[0].set_title('Original Image')
		ax[1].imshow(np.rot90(new)*-1,origin='lower',cmap='gray')
		ax[1].set_title('Edges (Gradient Magnitude)')
		ax[2].imshow(np.rot90(Theta),origin='lower',cmap='viridis')
		ax[2].set_title('Edges (Gradient Direction)')
		#ax[0,1].quiver(gx[::5,::5],gy[::5,::5])
		ax[0].axis('off')
		ax[1].axis('off')
		ax[2].axis('off')
		#ax[0,1].axis('off')
		if save:
			plt.savefig('sobel_cloudy.pdf'%(im),bbinches='tight')
		plt.show()
		plt.close()
		fig = plt.figure()
		plt.imshow(np.rot90(new)*-1,origin='lower',cmap='gray')
		plt.title('Edges -- cloudy (%.2f)'%(num_edges/(640*480.) * 100))
		plt.savefig('edges_cloudy.png',bbinches='tight')
		plt.axis('off')
		plt.close()

if __name__ == '__main__':
	h, TEST = rawFits('./assets/tesCurrentImage.FIT')
	cloudy = './assets/cloudy.FIT'
	clear = './assets/clear.FIT'
	moon = './assets/moon.FIT'
	h,TEST = rawFits('./assets/031854.FIT')
	im = -1*np.log(TEST)
	sobel(im,verbose=True,save=False)
	print brightness(im)