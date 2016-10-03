from ij.gui import GenericDialog
from ij.plugin import HyperStackConverter
from script.imglib.analysis import DoGPeaks
from script.imglib.color import Red
from script.imglib.algorithm import Scale2D
from script.imglib.math import Compute
from script.imglib import ImgLib
from ij3d import Image3DUniverse
from javax.vecmath import Color3f, Point3f
from ij.process import StackStatistics as IS
from ij.io import DirectoryChooser
from Utilities import Counter3D
from ij.plugin.frame import RoiManager
import csv
import os
import fnmatch
from ij import VirtualStack, IJ, CompositeImage
from ij.process import ColorProcessor
from ij.io import DirectoryChooser
from ij.gui import YesNoCancelDialog
from mpicbg.imglib.image import ImagePlusAdapter
from mpicbg.imglib.algorithm.fft import PhaseCorrelation
from javax.vecmath import Point3i
from java.io import File, FilenameFilter
from shutil import rmtree
from loci.plugins import BF
from ij import WindowManager as WM


#the drift correction of stacks from a time series was adapted from the code by Robert Bryson-Richardson and Albert Cardona (see below)
#------------------------------------------------------------------------------------------------------------------------------------
# Robert Bryson-Richardson and Albert Cardona 2010-10-08 at Estoril, Portugal
# EMBO Developmental Imaging course by Gabriel Martins
#
# Register time frames (stacks) to each other using Stitching_3D library
# to compute translations only, in all 3 spatial axes.
# Operates on a virtual stack.




# imp stands for ij.ImagePlus instance


def compute_stitch(imp1, imp2):
	""" Compute a Point3i that expressed the translation of imp2 relative to imp1."""
	phc = PhaseCorrelation(ImagePlusAdapter.wrap(imp1), ImagePlusAdapter.wrap(imp2), 5, True)
	phc.process()
	return Point3i(phc.getShift().getPosition())


def extract_frame(imp, frame, channel):
	""" From a VirtualStack that is a hyperstack, contained in imp,
	extract the timepoint frame as an ImageStack, and return it.
	It will do so only for the given channel. """
	stack = imp.getStack() # multi-time point virtual stack
	vs = ImageStack(imp.width, imp.height, None)
	for s in range(1, imp.getNSlices()+1):
		i = imp.getStackIndex(channel, s, frame)
		vs.addSlice(str(s), stack.getProcessor(i))
	return vs

def compute_frame_translations(imp, channel):
	""" imp contains a hyper virtual stack, and we want to compute
	the X,Y,Z translation between every time point in it
	using the given preferred channel. """
	t1_vs = extract_frame(imp, 1, channel)
	shifts = []
	# store the first shift: between t1 and t2
	shifts.append(Point3i(0, 0, 0))
	# append the rest:
	IJ.showProgress(0)
	i = 1
	for t in range(2, imp.getNFrames()+1):
		t2_vs = extract_frame(imp, t, channel)
		shift = compute_stitch(ImagePlus("1", t1_vs), ImagePlus("2", t2_vs))
		shifts.append(shift)
		t1_vs = t2_vs
		IJ.showProgress(i / float(imp.getNFrames()))
		i += 1
	IJ.showProgress(1)
	return shifts

def concatenate_shifts(shifts):
	""" Take the shifts, which are relative to the previous shift,
	and sum them up so that all of them are relative to the first."""
	# the first shift is 0,0,0
	for i in range(2, len(shifts)): # we start at the third
		s0 = shifts[i-1]
		s1 = shifts[i]
		s1.x += s0.x
		s1.y += s0.y
		s1.z += s0.z
	return shifts

def compute_min_max(shifts):
	""" Find out the top left up corner, and the right bottom down corner,
	namely the bounds of the new virtual stack to create.
	Expects absolute shifts. """
	minx = Integer.MAX_VALUE
	miny = Integer.MAX_VALUE
	minz = Integer.MAX_VALUE
	maxx = -Integer.MAX_VALUE
	maxy = -Integer.MAX_VALUE
	maxz = -Integer.MAX_VALUE
	for shift in shifts:
		minx = min(minx, shift.x)
		miny = min(miny, shift.y)
		minz = min(minz, shift.z)
		maxx = max(maxx, shift.x)
		maxy = max(maxy, shift.y)
		maxz = max(maxz, shift.z)

	return minx, miny, minz, maxx, maxy, maxz

def zero_pad(num, digits):
	""" for 34, 4 --> '0034' """
	str_num = str(num)
	while (len(str_num) < digits):
		str_num = '0' + str_num
	return str_num

def create_registered_hyperstack(imp, target_folder, channel):
	""" Takes the imp, which contains a virtual hyper stack,
	and determines the x,y,z drift for each pair of time points,
	using the preferred given channel,
	and output one image for each slide into the target folder."""
	shifts = compute_frame_translations(imp, channel)
	# Make shifts relative to 0,0,0 of the original imp:
	shifts = concatenate_shifts(shifts)
	print "shifts concatenated:"
	for s in shifts:
		print s.x, s.y, s.z
	# Compute bounds of the new volume,
	# which accounts for all translations:
	minx, miny, minz, maxx, maxy, maxz = compute_min_max(shifts)
	# Make shifts relative to new canvas dimensions
	# so that the min values become 0,0,0
	for shift in shifts:
		shift.x -= minx
		shift.y -= miny
		shift.z -= minz
	print "shifts relative to new dimensions:"
	for s in shifts:
		print s.x, s.y, s.z
	# new canvas dimensions:
	width = imp.width + maxx - minx
	height = maxy - miny + imp.height
	slices = maxz - minz + imp.getNSlices()

	print "New dimensions:", width, height, slices
	# Count number of digits of each dimension, to output zero-padded numbers:
	slice_digits = len(str(slices))
	frame_digits = len(str(imp.getNFrames()))
	channel_digits = len(str(imp.getNChannels()))
	# List to accumulate all created names:
	names = []
	# Prepare empty slice to pad in Z when necessary
	empty = imp.getProcessor().createProcessor(width, height)
	# if it's RGB, fill the empty slice with blackness
	if isinstance(empty, ColorProcessor):
		empty.setValue(0)
		empty.fill()
	# Write all slices to files:
	stack = imp.getStack()
	for frame in range(1, imp.getNFrames()+1):
		shift = shifts[frame-1]
		fr = "t" + zero_pad(frame, frame_digits)
 		# Pad with mpty slices before reaching the first slice
		for s in range(shift.z):
			ss = "_z" + zero_pad(s + 1, slice_digits) # slices start at 1
 			for ch in range(1, imp.getNChannels()+1):
 				name = fr + ss + "_c" + zero_pad(ch, channel_digits) +".tif"
				names.append(name)
				FileSaver(ImagePlus("", empty)).saveAsTiff(target_folder + "/" + name)
	# Add all proper slices
		for s in range(1, imp.getNSlices()+1):
			ss = "_z" + zero_pad(s + shift.z, slice_digits)
			for ch in range(1, imp.getNChannels()+1):
				ip = stack.getProcessor(imp.getStackIndex(ch, s, frame))
				ip2 = ip.createProcessor(width, height) # potentially larger
				ip2.insert(ip, shift.x, shift.y)
				name = fr + ss + "_c" + zero_pad(ch, channel_digits) +".tif"
				names.append(name)
				FileSaver(ImagePlus("", ip2)).saveAsTiff(target_folder + "/" + name)
	# Pad the end
		for s in range(shift.z + imp.getNSlices(), slices):
			ss = "_z" + zero_pad(s + 1, slice_digits)
			for ch in range(1, imp.getNChannels()+1):
				name = fr + ss + "_c" + zero_pad(ch, channel_digits) +".tif"
				names.append(name)
				FileSaver(ImagePlus("", empty)).saveAsTiff(target_folder + "/" + name)

	# Create virtual hyper stack with the result
	vs = ImageStack(width, height, None)
	for name in names:
		vs.addSlice(IJ.openImage(target_folder+"/"+name).getProcessor())
	vs_imp = ImagePlus("registered time points", vs)
	vs_imp.setDimensions(imp.getNChannels(), len(names) / (imp.getNChannels() * imp.getNFrames()), imp.getNFrames())
	vs_imp.setOpenAsHyperStack(True)
	IJ.log("\nHyperstack dimensions: time frames:" + str(vs_imp.getNFrames()) + ", slices: " + str(vs_imp.getNSlices()) + ", channels: " + str(vs_imp.getNChannels()))
	if 1 == vs_imp.getNSlices():
		return vs_imp
	# Else, as composite
	mode = CompositeImage.COLOR;
	if isinstance(imp, CompositeImage):
		mode = imp.getMode()
	else:
		return vs_imp
	return CompositeImage(vs_imp, mode)

class Filter(FilenameFilter):
	def accept(self, folder, name):
		return not File(folder.getAbsolutePath() + "/" + name).isHidden()

def validate(target_folder):
	f = File(target_folder)
	if len(File(target_folder).list(Filter())) > 0:
		yn = YesNoCancelDialog(IJ.getInstance(), "Warning!", "Target folder is not empty! May overwrite files! Continue?")
		if yn.yesPressed():
			return True
		else:
      			return False
  	return True

def run(imp,folder,channel):

	if imp is None:
		return
	if not imp.isHyperStack():
 		print "Not a hyper stack!"
		return
	if 1 == imp.getNFrames():
		print "There is only one time frame!"
		return
	if 1 == imp.getNSlices():
		print "To register slices of a stack, use 'Register Virtual Stack Slices'"
		return
 # dc = DirectoryChooser("Choose target folder")
 # target_folder = dc.getDirectory()
	target_folder = folder+'slices/'
	os.mkdir(target_folder)
	if target_folder is None:
		return # user canceled the dialog
	if not validate(target_folder):
		return
	#gd = GenericDialog("Options")
	#channels = []
	#for ch in range(1, imp.getNChannels()+1 ):
	#	channels.append(str(ch))
	#gd.addChoice("channel:", channels, channels[0])
	#gd.showDialog()
	#if gd.wasCanceled():
	#	return
	#channel = gd.getNextChoiceIndex() + 1  # zero-based
	vs_imp = create_registered_hyperstack(imp, target_folder, channel)
	return vs_imp

# -----------------------------------------------------------------------------------------------------------------------------



def remove_dark_timepoints(imp, channel, threshold):
	nc = imp.getNChannels();
	nt = imp.getNFrames();
	ns = imp.getNSlices();
	indices = [];

	ndeleted = 0;
 	for t in range(1, imp.getNFrames() + 1):
 		mean = 0;
		for z in range(1, imp.getNSlices() + 1):
			idx = imp.getStackIndex(channel, z, t);
			ip = imp.getStack().getProcessor(idx);
			m = ip.getStatistics().max;
			if m > mean:
				mean = m
		if mean < threshold:
			ndeleted += 1
			black.append(t)
			for z in range(1, imp.getNSlices() + 1):
				for c in range(1, imp.getNChannels() + 1):
					indices.append(imp.getStackIndex(c, z, t));


	indices.sort();
	indices.reverse();


	for idx in indices:
		print("here")
		imp.getStack().deleteSlice(idx);

	print(nc);
	print(ns);
	print(nt);
	print(ndeleted);
	print(black)
	imp.setDimensions(nc, ns, nt - ndeleted);

def extract_time_points(imp, t0, dt):
	nc = imp.getNChannels()
	nz = imp.getNSlices()
	nt = imp.getNFrames()
	new = ImageStack(imp.width, imp.height)
	for t in range(t0,t0 + dt):
		for z in range (1,nz+1):
			for c in range(1,nc+1):
				idx=imp.getStackIndex(c,z,t);
				slice = imp.getStack().getProcessor(idx)
				cp = slice.duplicate()
				new.addSlice(imp.getStack().getSliceLabel(idx),cp)

	new_imp = ImagePlus("test" + imp.title, new)
	new_imp.setDimensions(nc,nz,dt)
	new_imp.setOpenAsHyperStack(True)
	return new_imp

def extract_channel(imp, t0, dt,ch):
	nc = imp.getNChannels()
	nz = imp.getNSlices()
	nt = imp.getNFrames()
	new = ImageStack(imp.width, imp.height)
	for t in range(t0,t0 + dt):
		for z in range (1,nz+1):
			for c in range(1,nc+1):
				if c == int(ch):
					idx=imp.getStackIndex(c,z,t);
					slice = imp.getStack().getProcessor(idx)
					cp = slice.duplicate()
					new.addSlice(imp.getStack().getSliceLabel(idx),cp)

	new_imp = ImagePlus("test" + imp.title, new)
	new_imp.setDimensions(1,nz,dt)
	new_imp.setOpenAsHyperStack(True)
	return new_imp

def print_means(imp, channel):
	indices = [];
	for t in range(1, imp.getNFrames() + 1):
		mean = 0;
		for z in range(1, imp.getNSlices() + 1):
			idx = imp.getStackIndex(channel, z, t);
			ip = imp.getStack().getProcessor(idx);
			m = ip.getStatistics().max;
			if m > mean:
				mean = m
		print(mean);

def print_histogram(imp, channel):
	indices = [];
	for t in range(1, imp.getNFrames() + 1):
		mean = 0;
		for z in range(1, imp.getNSlices() + 1):
			idx = imp.getStackIndex(channel, z, t);
			ip = imp.getStack().getProcessor(idx);
			m = ip.getHistogram();
		print(m[0]);

def acc(lst):
	cnt = 0
	res = []

	for i in lst:
		cnt += i
		res.append(cnt)

	return res

def analyze_cells(imp, size, i, zmin):

	test = []
	cdf = []

	#ip = imp.getProcessor()
	#print "grabbing image..."
	#test = ip.getHistogram()
	#test =  StackStatistics(imp).getHistogram()
	#print "calculating stack statistics..."
	#total = sum(test)
	#print "calculate threshold"
	#cdf = map(lambda x: x/float(total), acc(test))
	#thresh =  min(filter(lambda x: cdf[x] >= quantile, xrange(1,len(cdf)) ))

	max_int = StackStatistics(imp).max
	cal= imp.getCalibration()
	scale2D = cal.pixelWidth / cal.pixelDepth
	sigma = (size / cal.pixelWidth) * scale2D
	iso = Compute.inFloats(Scale2D(ImgLib.wrap(imp),scale2D))
	peaks = DoGPeaks(iso, sigma, sigma * 0.5, thresh, 1)
	print "FOund", len(peaks), "peaks"

	ps = []
	file = open(folder+str(i).zfill(4)+'_test_out.csv','w')
	exporter = csv.writer(file)

	for peak in peaks:
		if peak[2]>=zmin:
			print "raw",peak
			p = Point3f(peak)
			p.scale(cal.pixelWidth * 1/scale2D)
			print "sf", cal.pixelWidth * 1/scale2D
			print "scaled", p
			ps.append(p)
			t = ()
			exporter.writerow([p.x, p.y, p.z])

	file.close()

	if vis:
		iso = Compute.inFloats(Scale2D(Red(ImgLib.wrap(imp)),scale2D))
		univ = Image3DUniverse(512,512)
		univ.addIcospheres(ps,Color3f(1,0,0), 2, size/2, "Cells").setLocked(True)
		univ.addOrthoslice(imp).setLocked(True)
		univ.show()

def extract_structures(imp, folder, ost_ch, vasc_ch):
	n = imp.getNFrames()
	name = imp.getTitle()
	if struct_exp:
		imp2 = extract_channel(imp,1,n,ost_ch)
		imp2.show()
		IJ.run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		IJ.run("Z Project...", "start=1 stop="+str(n)+" projection=[Median] all");
		imp2 = IJ.getImage()
		IJ.run(imp2,"Auto Threshold", "method=Moments white stack use_stack_histogram");
		IJ.run(imp2,"Invert", "stack");
		IJ.run("Exact Euclidean Distance Transform (3D)");
		imp2 = IJ.getImage()
		FileSaver(ImagePlus("export",imp2.getStack())).saveAsTiffStack(folder+exp_name+"_osteoblasts.tif")
		imp2.close()

	if vasc_exp:
		imp2 = extract_channel(imp,1,n,vasc_ch)
		imp2.show()
		IJ.run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		IJ.run("Z Project...", "start=1 stop="+str(n)+" projection=[Median] all");
		imp2 = IJ.getImage()
		IJ.run(imp2,"Auto Threshold", "method=Moments white stack use_stack_histogram");
		IJ.run(imp2,"Invert", "stack");
		IJ.run("Exact Signed Euclidean Distance Transform (3D)");
		imp2 = IJ.getImage()
		FileSaver(ImagePlus("export",imp2.getStack())).saveAsTiffStack(folder+exp_name+"_vasculature.tif")
		imp2.close()


def extract_osteoblast_structures(imp,folder,channel):
	n = imp.getNFrames()
	imp2 = extract_channel(imp,1,n,channel)
	imp.close()
	imp2.show()

	IJ.run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	#IJ.run("Z Project...", "start=1 stop="+str(n)+" projection=[Max Intensity] all");
	IJ.run("Z Project...", "start=1 stop="+str(n)+" projection=[Median] all");

	imp2 = IJ.getImage()
	#IJ.run(imp2,"Auto Threshold", "method=MaxEntropy white stack use_stack_histogram");
	IJ.run(imp2,"Auto Threshold", "method=Moments white stack use_stack_histogram");
	IJ.run(imp2,"Invert", "stack");
	IJ.run("Exact Euclidean Distance Transform (3D)");
	imp2 = IJ.getImage()
	FileSaver(ImagePlus("export",imp2.getStack())).saveAsTiffStack(folder+exp_name+"_structures.tif")
	imp2.close()

#TODO implement vasculature detection such that it works with the same opened registered image as osteoblast detection !
def extract_vasculature_structures(imp,folder,channel):
	n = imp.getNFrames()
	imp2 = extract_channel(imp,1,n,channel)
	imp.close()
	imp2.show()

	IJ.run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	#IJ.run("Z Project...", "start=1 stop="+str(n)+" projection=[Max Intensity] all");
	IJ.run("Z Project...", "start=1 stop="+str(n)+" projection=[Median] all");

	imp2 = IJ.getImage()
	#IJ.run(imp2,"Auto Threshold", "method=MaxEntropy white stack use_stack_histogram");
	IJ.run(imp2,"Auto Threshold", "method=Moments white stack use_stack_histogram");
	IJ.run(imp2,"Invert", "stack");
	IJ.run("Exact Euclidean Distance Transform (3D)");
	imp2 = IJ.getImage()
	FileSaver(ImagePlus("export",imp2.getStack())).saveAsTiffStack(folder+exp_name+"_structures.tif")
	imp2.close()

def detect_objects_3D(imp, size, i, zmin, zmax):

	test = []
	cdf = []

	IJ.run("3D OC Options", "  nb_of_obj._voxels median_gray_value centroid dots_size=5 font_size=10 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=none");

	#ip = imp.getProcessor()
	#print "grabbing image..."
	#test = ip.getHistogram()
	#test =  StackStatistics(imp).getHistogram()
	#print "calculating stack statistics..."
	#total = sum(test)
	#print "calculate threshold"
	#cdf = map(lambda x: x/float(total), acc(test))
	#thresh =  min(filter(lambda x: cdf[x] >= quantile, xrange(1,len(cdf)) ))



	max_int = StackStatistics(imp).max

	cal= imp.getCalibration()
	scale2D = cal.pixelWidth / cal.pixelDepth
	sigma = (size / cal.pixelWidth) * scale2D


	iso = Compute.inFloats(Scale2D(ImgLib.wrap(imp),scale2D))
	#thresh = 0
	thresh = round(max_int/25.5)
	peaks = DoGPeaks(iso, sigma, sigma * 0.5, thresh, 1)
	ps = []
 #	print max_int
#	for p in peaks:
#		if p[2]>zmin and p[2]<zmax:
#			p2 = Point3f(p)
#			p2.scale(cal.pixelWidth * 1/scale2D)
#			ps.append(p2)


	#ob3d = Counter3D(imp, thresh)
	#peaks = ob3d.getCentroidList()
	print "FOund", len(peaks), "peaks"

	#ps = []
	file = open(folder+exp_name+'/'+'cells/'+str(i).zfill(4)+'_cells.csv','w')
	exporter = csv.writer(file)

#TODO: check consistency of coordinates !

	for peak in peaks:
		print peak[2]
		if peak[2]>=zmin and peak[2]<=zmax:
			#print "raw",peak
			p = Point3f(peak)
			p.scale(cal.pixelWidth * 1./scale2D)
			#print "sf", cal.pixelWidth * 1/scale2D
			#print "scaled", p
			ps.append(p)
			t = ()
			exporter.writerow([p.x, p.y, p.z])

	file.close()

#	if vis:
#		iso = Compute.inFloats(Scale2D(Red(ImgLib.wrap(imp)),scale2D))
#		univ = Image3DUniverse(512,512)
#		univ.addIcospheres(ps,Color3f(1,0,0), 2, size/2, "Cells").setLocked(True)
#		univ.addOrthoslice(imp).setLocked(True)
#		univ.show()

#TODO select channel for vasculature and osteoblasts

def get_opts():
	gd = GenericDialog("setup plugin")
	#channels = map(lambda x: str(x), range(1,imp.getNChannels()+1))
	channels = ['1','2','3'];
	modes = ['batch mode', 'manual']
	gd.addChoice("processing mode",modes,modes[1])
	gd.addChoice("choose channel for pre-registration: ", channels, channels[0])
	gd.addChoice("choose channel for registration: ", channels, channels[2])
	gd.addChoice("channel for cell detection: ", channels, channels[1])
	gd.addChoice("channel for osteoblast detection: ", channels, channels[2])
	gd.addChoice("channel for vasculature detection: ", channels, channels[0])
#	gd.addNumericField("quantile for cell segmentation", 0.9995, 36)
	gd.addNumericField("rough cell size in microns", 15., 2)
#	gd.addSlider("minimal z depth (voxels) for detected cells", 1, imp.getNSlices()+1, 8)
	gd.addSlider("minimal z depth (voxels) for detected cells", 1,36, 5)
	gd.addSlider("minimal z depth (voxels) for detected cells", 1,36, 25)

	gd.addCheckbox("delete black frames (recommended)", True)
	gd.addChoice("channel to estimate black frames: ", channels, channels[2])
	gd.addCheckbox("rough registration of time series", True)
	gd.addCheckbox("descriptor-based registration of time series", True)
	gd.addCheckbox("detect cells", False)
	gd.addCheckbox("export osteoblast structures", False)
	gd.addCheckbox("export vasculature structures", False)
	gd.addCheckbox("save registered movie", True)
	gd.addCheckbox("delete temporary files", True)
	#gd.addCheckbox("show 3d vis", False)

	gd.showDialog()
	proc_mode = gd.getNextChoice()
	ch_pre = gd.getNextChoice()
	reg = gd.getNextChoice()
	cc = gd.getNextChoice()
	ost_ch = gd.getNextChoice()
	vasc_ch = gd.getNextChoice()
	#quantile = gd.getNextNumber()
	size = gd.getNextNumber()
	minz = gd.getNextNumber()
	zmax = gd.getNextNumber()
	del_black = gd.getNextBoolean();
	black_ch = gd.getNextChoice()
	preregflag = gd.getNextBoolean();
	regflag = gd.getNextBoolean()
	detect_cells = gd.getNextBoolean()
	struct_exp = gd.getNextBoolean()
	vasc_exp = gd.getNextBoolean()
	save_res = gd.getNextBoolean()
	del_tmp = gd.getNextBoolean()
	#vis = gd.getNextBoolean()

	return proc_mode, ch_pre, reg, cc, ost_ch, vasc_ch, size, minz, zmax, del_black, black_ch, preregflag, regflag, detect_cells, struct_exp, vasc_exp, save_res, del_tmp

def proc_movie(imp, exp_name, olength, pwidth, pheight, pdepth):
	if opts is not None:
		reg = opts
		ch0 = int(reg)

	if del_black:
		IJ.log("removing dark stacks")
		remove_dark_timepoints(imp, int(black_ch), 150);
		file = open(folder+exp_name+'/'+'black_frames.csv','w')
		exporter = csv.writer(file)
		exporter.writerow(black)
		file.close()

	if preregflag:
		IJ.log("pre-registration of movie")
		imp_pre_reg = run(imp,folder,int(ch_pre))
		imp_pre_reg.show()
		imp.close()
		if save_res:
			imp_tmp=ImagePlus("export",imp_pre_reg.getStack())
			imp_tmp.setDimensions(imp_pre_reg.getNChannels(), imp_pre_reg.getNSlices(),imp_pre_reg.getNFrames())
			imp_tmp.setOpenAsHyperStack(True)
			FileSaver(imp_tmp).saveAsTiffStack(folder+exp_name+'/'+exp_name+"_prereg.tif")
			imp_tmp.close()
		if del_tmp:
			rmtree(folder+'slices/')
		IJ.log("done...")
		imp = imp_pre_reg


	if regflag:
		IJ.log("descriptor-based registration")
		IJ.run(imp,"Descriptor-based series registration (2d/3d + t)" ,"series_of_images=["+imp.getTitle()+"] brightness_of=Strong approximate_size=[3 px] type_of_detections=[Maxima only] transformation_model=[Translation (3d)] number_of_neighbors=3 redundancy=1 significance=3 allowed_error_for_ransac=5 global_optimization=[All-to-all matching with range ('reasonable' global optimization)] range=5 choose_registration_channel="+str(ch0) +" image=[Fuse and display]")
		IJ.log("done...")
		imp_reg = WM.getImage("[XYZCT] registered "+imp.getTitle())
	#	imp = WM.getImage("registered time points")
		imp.close()
		if save_res:
				imp_tmp=ImagePlus("export",imp_reg.getStack())
				imp_tmp.setDimensions(imp_reg.getNChannels(), imp_reg.getNSlices(),imp_reg.getNFrames())
				imp_tmp.setOpenAsHyperStack(True)
				FileSaver(imp_tmp).saveAsTiffStack(folder+exp_name+'/'+exp_name+"_registered.tif")
				imp_tmp.close()

		imp = imp_reg

#TODO: fine-tune and test cell detection, feed size constraints and different options...
#TODO: test, how DoG filter works on the data

	if detect_cells:
		IJ.log("detecting cells...")

		dict ={}
		index_lst = range(1,olength+1)
		count_del = 0

		for d in black:
			del index_lst[d-1-count_del]
			count_del+=1

		for el in range(1,olength+1):
			if el in black:
				dict[el] = -1
			else:
				dict[el] = index_lst.index(el)

		if not os.path.exists(folder+exp_name+'/'+'cells/'):
			os.mkdir(folder+exp_name+'/'+'cells/')

		for i in xrange(1,olength+1):
			if dict[i] > -1:
				tmp = extract_channel(imp,dict[i],1,cc)
				cal = tmp.getCalibration()
				cal.pixelWidth = pw
				cal.pixelHeight = ph
				cal.pixelDepth=pd
				tmp.setCalibration(cal)
				detect_objects_3D(tmp, size, i, zmin, zmax)
				tmp.close()
	#		analyze_cells(tmp,quantile,size,dict[i],zmin)
			else:
				file = open(folder+exp_name+'/'+'cells/'+str(i).zfill(4)+'_cells.csv','w')
				file.close()


	if struct_exp or vasc_exp:
		IJ.log("segment spatial image structures")
		imp = IJ.getImage()
		extract_structures(imp, folder+exp_name+'/', int(ost_ch), int(vasc_ch))
		#extract_osteoblast_structures(imp,folder+exp_name+'/',int(ost_ch))
'''
	if vasc_exp:
		IJ.log("segment vasculature")
		imp = IJ.getImage()
		extract_vaculature_structures(imp,folder+exp_name+'/',int(vasc_ch))
'''

#TODO: fix problems with stack calibrations throughout processing
#TODO: make single parts workable in isolation, i.e. check if reg not selected, then get reg. movie and do the tracking etc...
# main part
proc_mode, ch_pre, opts, cc, ost_ch, vasc_ch, size, zmin, zmax, del_black,  black_ch, preregflag, regflag, detect_cells, struct_exp, vasc_exp, save_res, del_tmp = get_opts()

zmin = False

#if quantile > 1:
#	quantile = 0.99999
#if quantile <0:
#	quantile = 0

if proc_mode == 'manual':
	black=[]
	dc = DirectoryChooser("Choose a folder to save results")
	folder = dc.getDirectory()

	if folder is None:
		print "urgh"
	else:
		print "let's go detecting... ", folder

	imp = IJ.getImage();
	exp_name = imp.getTitle()[:-4];

	if not os.path.exists(folder+exp_name+'/'):
		os.mkdir(folder+exp_name+'/')

	olength = imp.getNFrames();

	cal = imp.getCalibration()
	pw = cal.pixelWidth
	ph = cal.pixelHeight
	pd = cal.pixelDepth

	c_file = open(folder+exp_name+'/'+'calibration.csv','w')
	exporter = csv.writer(c_file)
	exporter.writerow([pw,ph,pd])
	c_file.close()

	proc_movie(imp, exp_name, olength, pw, ph, pd)
	to_clean = WM.getIDList()
	for i in to_clean:
		imp = WM.getImage(i)
		FileSaver(imp).saveAsTiffStack(folder+exp_name+'/trash.tif')
		imp.close()

else:
	dc = DirectoryChooser("choose folder for processing")
	folder = dc.getDirectory()

	files = os.listdir(folder)

	for m_file in files:
		if m_file.endswith(".tiff"):
			black=[]
			imp = BF.openImagePlus(folder+m_file)[0]
			imp.show()
			exp_name = imp.getTitle()[:-4]
			if not os.path.exists(folder+exp_name+'/'):
				os.mkdir(folder+exp_name+'/')
			olength = imp.getNFrames();

			cal = imp.getCalibration()
			pw = cal.pixelWidth
			ph = cal.pixelHeight
			pd = cal.pixelDepth

			c_file = open(folder+exp_name+'/'+'calibration.csv','w')
			exporter = csv.writer(c_file)
			exporter.writerow([pw,ph,pd])
			c_file.close()

			proc_movie(imp, exp_name, olength, pw, ph, pd)
			imp.close()

			to_clean = WM.getIDList()
			for i in to_clean:
				imp = WM.getImage(i)
				FileSaver(imp).saveAsTiffStack(folder+exp_name+'/trash.tif')
				imp.close()

#TODO export settings as csv



to_clean = WM.getIDList()
for i in to_clean:
	imp = WM.getImage(i)
	FileSaver(imp).saveAsTiffStack(folder+exp_name+'/trash.tif')
	imp.close()
IJ.log("_____________________________________done____________________________________________")
