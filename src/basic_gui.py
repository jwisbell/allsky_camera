import sys
import time
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
from subprocess import call
from argparse import ArgumentParser
from astropy.io import fits


def rawFits(filename):
    f = fits.open(filename)
    h = f[0].header
    d = f[0].data
    f.close
    return h, d

h, TEST = rawFits('./assets/tesCurrentImage.FIT')
im = -1*np.log(TEST)

class ImageViewer(QtGui.QWidget):
	 def __init__(self, parent=None):
		super(ImageViewer,self).__init__(parent)
		#### Create Gui Elements ###########
		self.mainbox = QtGui.QWidget()
		self.setLayout(QtGui.QGridLayout())
		self.setCursor(QtCore.Qt.CrossCursor)

		#create the plotting canvas and define its size
		self.canvas = pg.GraphicsLayoutWidget()
		self.layout().addWidget(self.canvas)
		self.view = self.canvas.addViewBox()
		self.view.setAspectLocked(True)
		self.view.setRange(xRange=[0,480],yRange=[0,640],padding=-1)
		self.view.setMouseEnabled(x=False,y=False)
		self.view.setBackgroundColor((.5,.5,.5,1.))
		#self.xsize = float(xspan); self.ysize = float(yspan)

		#define the geometry of the entire widget
		self.setGeometry(0,0,480,640)  
		self.canvas.nextRow()
		#self.setContentsMargins(lmargin,tmargin,rmargin,bmargin)

		#intitialize the image data
		self.img = pg.ImageItem(border=None)
		self.img.setZValue(-100) #set image below overlays
		self.img.setImage(im)#autoLevels=True)
		self.view.addItem(self.img)




"""Main widget that holds and organizes all other widgets (Surface, AboutScreen, UIowaScreen, Settings, WelcomeScreen) and 
all other threads (GravityThread). Called at the start of the program."""
class Display(QtGui.QWidget):
    def __init__(self, parent=None):
		super(Display,self).__init__(parent)
		#### Create Gui Elements ###########
		self.mainbox = QtGui.QWidget()
		self.setLayout(QtGui.QGridLayout())
		self.setCursor(QtCore.Qt.CrossCursor)
		self.setMouseTracking(True)
		self.setGeometry(0,0,640,640)
		self.setFocus()

		self.imviewer = ImageViewer(self)
		self.imviewer.move(640-480, 0)

		#dummy buttons
		self.stack_btn = QtGui.QPushButton('Stack',self)
		self.toggle_btn = QtGui.QPushButton('Overlay',self)
		self.animate_btn = QtGui.QPushButton('Animate',self)
		self.test_btn = QtGui.QPushButton('Button',self)

		self.stack_btn.move(25,100)
		self.toggle_btn.move(25,150)
		self.animate_btn.move(25,200)
		self.test_btn.move(25,250)



#boilerplate code, initialize the Display widget and start the program
if __name__ == '__main__':
    #args = parser.parse_args()
    app = QtGui.QApplication(sys.argv)
    thisapp = Display()
    thisapp.show()
    sys.exit(app.exec_())