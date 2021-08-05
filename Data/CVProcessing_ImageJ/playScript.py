
import os
import datetime
import csv
import glob
import string
from ij import IJ, ImagePlus
from ij.measure import ResultsTable
from ij.process import ImageProcessor, FloatProcessor
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.plugin.filter import BackgroundSubtracter
from ij.plugin import MontageMaker, Thresholder


# process all cellavista images
# parentFolder contains entire cellavista experiment

def main():

	allparents = ['/media/matthew/DataDisk/CellavistaExperiments_Images/20160422_MDAMB231_Dox_6hr/', '/media/matthew/DataDisk/CellavistaExperiments_Images/20160422_MDAMB231_Dox_12hr/', '/media/matthew/DataDisk/CellavistaExperiments_Images/20160422_MDAMB231_Dox_24hr/']
	# be careful with 231 data, also have FUCCI channel
	for parentFolder in allparents:
		expName_start = parentFolder.rfind('/',0,len(parentFolder)-1)
		print(parentFolder)
		#outputPath = '/home/matthew/Documents/Vanderbilt/Gs/CVProcessing_ImageJ/' + parentFolder[expName_start+1:len(parentFolder)-1] + '.csv'
		outputPath = '/home/matthew/Documents/Vanderbilt/Gs/CVProcessing_ImageJ/dummy.csv'
		processExperiment(parentFolder, outputPath)


def processExperiment (parentFolder, outputCSV): 

	# go through parentFolder directory to collect all experimental timepoints
	allTPs = []
	for d in os.listdir(parentFolder):
		if os.path.isdir(parentFolder + d):
			allTPs.append(int(d))
	allTPs.sort()
	allTPs = [1,2,3,4,5]
	
	WellLUT = ['A','B','C','D','E','F','G','H']
	allWells = []
	allWellNames = []
	for r in range(2,3):
		for c in range(11,12):
			allWellNames.append('%s%02d' %(WellLUT[r-1], c)) 
			allWells.append("-R%02d-C%02d" % (r, c))
	
	# need to edit out corners!
	# order of images for montage
	#theOrder = [3 2 9
	#	4 1 8
	#	5 6 7]
	#theOrder = [13 12 11 10 25,... 
	#    14 3 2 9 24,...
	#    15 4 1 8 23,...
	#    16 5 6 7 22,...
	#    17 18 19 20 21]
	#theOrder = [1000 12 11 10 1000,... 
	#    13 3 2 9 21,...
	#    14 4 1 8 20,...
	#    15 5 6 7 19,...
	#    1000 16 17 18 1000]
	#theOrder = [1000, 12, 11, 10, 1000, 13, 3, 2, 9, 21, 14, 4, 1, 8, 20, 15, 5, 6, 7, 19, 1000, 16, 17, 18, 1000]
	# collected 3 channels for the 231 data, sample every third
	theOrder = [1000, 34, 31, 28, 1000, 37, 7, 4, 25, 61, 40, 10, 1, 22, 58, 43, 13, 16, 19, 55, 1000, 46, 49, 52, 1000]
	#theOrder = [1, 12, 11, 10, 1, 13, 3, 2, 9, 21, 14, 4, 1, 8, 20, 15, 5, 6, 7, 19, 1, 16, 17, 18, 1]
	
	oneImageOrder = [1000, 12, 11, 10, 1000, 13, 3, 2, 9, 21, 14, 4, 1, 8, 20, 15, 5, 6, 7, 19, 1000, 16, 17, 18, 1000]
	twoImageOrder = [1000, 23, 21, 19, 1000, 25, 5, 3, 17, 41, 27, 7, 1, 15, 39, 29, 9, 11, 13, 37, 1000, 31, 33, 35, 1000]
	threeImageOrder = [1000, 34, 31, 28, 1000, 37, 7, 4, 25, 61, 40, 10, 1, 22, 58, 43, 13, 16, 19, 55, 1000, 46, 49, 52, 1000]
	
	# create csv writer to output variables
	f = open(outputCSV, 'wb')
	writer = csv.writer(f)
	writer.writerow(['Well', 'CellCount', 'MeanCellArea'])
	
	# particle analyzer filter settings
	#options = PA.DISPLAY_SUMMARY
	options = 0
	rt = ResultsTable()
	minSize = 50
	maxSize = 1000
	minCirc = .2
	maxCirc = 1
	p = PA(options, PA.AREA + PA.KURTOSIS, rt, minSize, maxSize, minCirc, maxCirc)
	p.setFontSize(12)
	p.setHideOutputImage(True)
	
	# montage settings
	columns = 5
	rows = 5
	scale = 1.0
	first = 1
	last = 25
	inc = 1
	borderWidth = 0
	labels = False
	
	# rolling ball subtraction settings
	rb_radius = 50
	rb_createBackground = False
	rb_lightBackground = False
	rb_useParaboloid = True
	rb_presmooth = False
	rb_correctCorners = False
	
	haveBackground = False
	
	old_tpi = -1000
	for tpi in allTPs:
	
		tpi_str = str(tpi)
		CurrentDay = parentFolder + tpi_str
		current_tpi = tpi
		# cycle through each well
		for wellIter, currentWell in enumerate(allWells):
			wellFiles = glob.glob(CurrentDay + '/*' + currentWell + '*')
			wellFiles.sort()
			newStack = 1
			
			# if in a new folder, print out the time in the csv
			if current_tpi != old_tpi:
				dateString = wellFiles[0]
				dateClip_end = dateString.find('-')
				dateClip_start = dateString.rfind('/')+1
				dateString = dateString[dateClip_start:dateClip_end]
				if tpi == 1:
					tp0 = datetime.datetime.strptime(dateString,'%Y%m%d%H%M%S')
				old_tpi = current_tpi
				ctp = datetime.datetime.strptime(dateString,'%Y%m%d%H%M%S')
				timeDiff = ctp-tp0
				timeDiff_hours = timeDiff.days * 24 + timeDiff.seconds/float(3600)
				imageTime = round(timeDiff_hours,4)
				writer.writerow([tpi, '', ''])
				writer.writerow(['Timespan(h):', imageTime, ''])
			
			if len(wellFiles) == 63:
				theOrder = threeImageOrder
			elif len(wellFiles) == 42:
				theOrder = twoImageOrder
			else:
				theOrder = oneImageOrder
				
	
			# collect all images for each well
			for cOrder in theOrder:
				try:			
					WellIm = wellFiles[cOrder-1]
					imp = IJ.openImage(WellIm)
				except IndexError:
					if haveBackground:
						imp = backgroundImage
					else:
						itList = 0		
						while True:
							try:
								WellIm = wellFiles[itList]
								imp = IJ.openImage(WellIm)
								BackgroundSubtracter().rollingBallBackground(imp.getProcessor(), rb_radius, True, rb_lightBackground, rb_useParaboloid, True, rb_correctCorners)
								haveBackground = True
								backgroundImage = imp
								break
							except IndexError:
								itList += 1
			        if newStack == 1:
					newStack += 1
					myStack = imp.createEmptyStack()
				ip = imp.getProcessor()
				myStack.addSlice(ip)
			ipStack = ImagePlus('ipStack',myStack)
	
			# create montage for each well, apply median filter, and subtract background
			montagedImage = MontageMaker().makeMontage2(ipStack, columns, rows, scale, first, last, inc, borderWidth, labels)
			montagedImage.getProcessor().medianFilter()
			montagedImage.show()

			#IJ.run(montagedImage, 'Subtract Background...', 'rolling = 50 stack')
			#BackgroundSubtracter().rollingBallBackground(montagedImage.getProcessor(), rb_radius, rb_createBackground, rb_lightBackground, rb_useParaboloid, rb_presmooth, rb_correctCorners)
	
			# threshold montaged image and run watershed segmentation
			#montagedImage.getProcessor().setThreshold(10, 1000,ImageProcessor.NO_LUT_UPDATE)
			#IJ.run(montagedImage, 'Convert to Mask', '')
			#IJ.run(montagedImage, 'Watershed', '')
			#p.analyze(montagedImage)
			
			#cellAreas = rt.getColumn(0)
			#print([allWellNames[wellIter], rt.size(), round(sum(cellAreas)/len(cellAreas),2)])
			#writer.writerow([allWellNames[wellIter], rt.size(), round(sum(cellAreas)/len(cellAreas),2)])
			#rt.reset()
					
	f.close()
	return True
	
main()
