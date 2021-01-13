#!/usr/bin/python
'''
                               ** dmcalc **
Estimates the Dispersion Measure (DM) from the data in psrfits file format.

Returns the DM value with its uncertainty and reduced chi-square from tempo2 
DM fit.

Dependencies 
-------------
PSRCHIVE with python interface: http://psrchive.sourceforge.net/
TEMPO2: https://bitbucket.org/psrsoft/tempo2
SKLEARN: https://scikit-learn.org/stable/install.html

Parameters
----------
file(s)    : Input file(s) in psrfits format

ephem      : Ephemeris (or parameter) file  of the  pulsar. This is  required 
             to update the model. It can be  given as a command line argument. 
             If it is available in "PWD/ephemerides" folder, one can use that.
             Giving the file with this option overrides the default one.

model      : Template profile for cross-correlating  with the observation  to
             obtain DM. It can be given as a command line argument, otherwise
             it will look  for a matching one in  "PWD/ephemerides" directory
             and if found, will use that instead. One can use this  option to
             override the default selection.
           
fscrunch   : int, optional, default: None. Factor for scrunching the frequency 
             channels before passing it to DM estimation.

b3fscrunch : int, optional, default: None. Factor for scrunching the BAND3 
             data of uGMRT before passing it to DM estimation.

b3fscrunch : int, optional, default: None. Factor for scrunching the BAND5 
             data of uGMRT before passing it to DM estimation.

offset     : float, optional, default: None. Fix for jump between BAND3 and 
             BAND5 of uGMRT bands. 

writeout   : bool, optional,  default: False. Writes  out the file  corrected 
             for DM in a default directory (PWD/PSRJ_{site}_final), using the
             following options to reduce the file.

plot       : bool, optional, default: True. Prints the data analysis plot in
             a PDF file. ToA rejection steps and DM corrected ToAs are shown
             in addition to  DM corrected frequency evolution of the profile.

ptoa       : bool, optional, default: False. Prints the outliers cleaned ToAs 
             to  a file in the TEMPO2  readable format, so that, if required, 
             it can be used for other purposes.
           
Fscrunch   : bool, optional, default: False. Collapse all frequency  channels
             to produce one profile.

Tscrunch   : bool, optional,  default: False. Collapse  all  sub-integrations
             to produce one profile.

tscrunch   : int, optional, default: None. Factor to  scrunch  sub-integrations
             for writing out the DM corrected file.
           
quiet      : bool, optional,  default: False. Supresses all print  statements
             except warnings and errors.

Returns
-------
Dispersion Measure with uncertainty.


Examples
--------
# (a) for DM estimation with files in default directories:
#
dmcalc.py inputfile.fits
#
# (c) to use different ephemeris and template files:
#
dmcalc.py -E ephemeris.par -M model.fits data_file.fits
#
# (d) to write the DM corrected fits file and ToAs:
#
./dmcalc2.py -w -ptoa inputfile.fits

'''


# import modules...
import os
import sys
import numpy as np
import psrchive
import argparse
import time
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec

start = time.time()

parser = argparse.ArgumentParser(description='Code for measuring in-band '+ 
                                 'DM for pulsar data in psrfits format.')
parser.add_argument('files', nargs='+', type=str, 
					help='The list of fits file(s) for processing')
parser.add_argument('-E', '--ephem', type=str, 
					help='Ephemeris file to update the model. Exits if not '+
					      'given or is not available in "PWD/ephemerides" '+
					      'directory')
parser.add_argument('-M', '--model', nargs='+', type=str,
					help='Model template for ToA generation. Exits if not '+ 
					     'given or is not available in "PWD/templates" '+
					     'directory')
parser.add_argument('-f','--fscrunch', type=int, default=1,
					help='Factor to scrunch the number of channels for '+ 
					     'doing DM estimation (Def: 1)')
parser.add_argument('-b3f','--b3fscrunch', type=int, default=1,
					help='Factor to scrunch the number of channels for '+ 
					     'band3 GMRT data (Def: 1)')
parser.add_argument('-b5f','--b5fscrunch', type=int, default=1,
					help='Factor to scrunch the number of channels for '+ 
					     'band5 GMRT data (Def: 1)')
parser.add_argument('-w','--writeout', action='store_true',
					help='Writes out the DM corrected file. Def: False')
parser.add_argument('-ptoa','--print_toas', action='store_true',
					help='Print the prefit ToAs to file in tempo2 format. '+
					     'Def: False')
parser.add_argument('-F','--Fscrunch', action='store_true',
					help='Fully scrunch the number of channels for the '+
						 'final output archive (Def: False)')
parser.add_argument('-T','--Tscrunch', action='store_true',
					help='Completely time scrunch all the integrations')
parser.add_argument('-t','--tscrunch', type=int, default=1,
					help='Factor to scrunch the number of integrations for '+ 
					     'the final output archive (Def: None)')
parser.add_argument('-o','--offset', type=float, default=0.670520675,
					help='Offset to shift band 5 ToAs (in secs)')
parser.add_argument('-q', '--quiet', action='store_true', 
							help='Only print warnings')


def main():
	
	# parses the input arguments
	args = parser.parse_args()

	# checks status of quiet and ptoa
	quiet=False
	if args.quiet:
		quiet=True
	tempo2=True
	ptoa=False
	if args.print_toas:
		ptoa=True
		
	if not quiet:
		print("Loading the archive files for DM estimation")

	# loads the psrfits file
	archives = []
	for filename in args.files:
		archives.append(psrchive.Archive_load(filename))
	narch = len(archives)
	if narch >= 1:
		if not quiet:
			print("Appending the archives ..."),
		# append data
		ar = freq_appendData(narch, archives, args.offset, 
							args.b3fscrunch, args.b5fscrunch)
		if not quiet:
			print(" done!")
	else:
		if not quiet:
			print("Only one archive was given, so nothing to frequency-append.")
	# ar is the final archive after performing frequency append
	ar = archives[0]
	del archives
	
	# extracts relevant information from the archive
	ar_psr = ar.get_source()
	ar_nbins = ar.get_nbin()
	ar_tel = ar.get_telescope()
	mjd_start=ar.get_Integration(0).get_start_time().in_days()
	mjd_end=ar.get_Integration(0).get_end_time().in_days()
	ar_mjd = mjd_start + (mjd_end-mjd_start)/2.
	length = ar.integration_length()
	ar.update_centre_frequency()
	ar_centfr = ar.get_centre_frequency()
	ar_nchan = ar.get_nchan()
	ar_bw = ar.get_bandwidth()
	ar_chnwdth = ar_bw / ar_nchan
	ffrac = args.fscrunch
	if not quiet:
		print("\nNow preparing for DM estimation\n")

	pwd=os.getcwd()

	# checks for ephemeris file and exit if not given or is not available
	# in the default directory "PWD/ephemerides".
	if args.ephem != None:
		ephemeris = args.ephem
	else:
		ephemeris = "ephemerides/"+ar_psr+".par"
		if not (os.path.exists(ephemeris)):
			sys.exit(1)
	if not quiet:
		print ("\nEphemeris file is:"+ephemeris+'\n')
	
	# if template is given as input argument load and process them
	model = []
	for filename in args.model:
		model.append(psrchive.Archive_load(filename))
	if args.model != None:
		if len(args.model) == 1:
			model = freq_appendModel(1,model,args.offset, args.b3fscrunch, args.b5fscrunch)
		if len(args.model) > 1:
			model = freq_appendModel(1,model,args.offset, args.b3fscrunch, args.b5fscrunch)
	# If the template is not given, looking for a matching template in the templates directory
	if args.model == None:
		if not quiet:
			print("Looking for matching template in templates directory..."),
		import subprocess
		tempdir="templates/*.sm"
		tempfile=ar_psr+'_tmp.txt'
		a=subprocess.call("psredit -c name,nbin,bw,nchan,freq -Q '%s' > '%s'"
							 % (tempdir,tempfile), shell=True)

		tempnchan=""
		t1=str(ar_nbins)
		if ar_tel=='gmrt':
			t2=str(int(ar_bw))
		else:
			t2=str((ar_bw))
		t3=('%.2f'%ar_centfr)
		f = open(tempfile,'r')
		for line in f:
			line = line.strip()
			columns=line.split()
			t4 = float(columns[5])
			t4 = ('%.2f'%t4)
			if ar_tel=='gmrt':
				if (columns[1]==ar_psr and columns[2]==t1 and str(int(columns[3]))==t2 and t4==t3):
					modeltempl=columns[0]
					tempnchan=columns[4]
					if not quiet:
						print (' done\n')
			else:
				if (columns[1]==ar_psr and columns[2]==t1 and str((columns[3]))==t2 and t4==t3):
					modeltempl=columns[0]
					tempnchan=columns[4]
					if not quiet:
						print (' done\n')
		if modeltempl=='' and tempnchan=='':
			
			print("\n** No matching template found for DM fitting. Exiting. **\n")
			sys.exit(1)
		f.close()
		os.remove(tempfile)
		if not quiet:
			print("Found matching template: "+modeltempl)
		model.append(psrchive.Archive_load(modeltempl))
	if not quiet:
		print("\nEstimating the DM from the observation")
	model.update_centre_frequency()

	# cloning the original file for passing to DMCalc() routine
	arch = ar.clone()

	# Calling the DM estimation routine	
	dmval, dmverr, fitchisq, pre_rms, post_rms, ToA_Err  = DMCalc(arch, ar_nchan, ar_centfr, 
														   ar_bw, ar_psr, ar_tel, ar_mjd, model, 
									 					   ephemeris, pwd, ffrac, quiet, tempo2,
														   ptoa, narch)
	
	# writing out the final DM corrected file, if requested
	if args.writeout:
		# removing the DM and DMEPOCH from the ephemeris file for uptation
		infile = open(ephemeris,"r")
		tmpeph = ar_psr+'.eph'
		output = open(tmpeph,"w+")
		for i, line in enumerate(infile):
			if not line.lstrip().startswith('DM'):
					if not line.lstrip().startswith('DMEPOCH'):
						output.write(line)
		infile.close()
		output.close()
		# updating the ephemeris file with measured DM
		dmline = "DM			  "+str(dmval)+"\t\t"+str(dmverr)
		dmepochline  = "DMEPOCH		 "+str(round(ar_mjd,2))
		if not args.quiet:
			print("Updating the ephemeris with new DM... "),
		f = open(tmpeph,'a')
		f.write("%s\n %s\n" % (dmline, dmepochline))
		if not args.quiet:
			print(" done!")
		f.close()

		# updating the ephemeris in the archive with the measured DM
		if not quiet:
			print("Correcting the DM of the observed file and writing it out... "),
		os.remove(tmpeph)
		# creating the directory for writing the file
		dirfinal=os.path.join(pwd,ar_psr+"_"+ar_tel+"_final")
		if not os.path.exists(dirfinal):
			os.makedirs(dirfinal)
		# filename with path of the DM corrected file
		outfile = dirfinal+"/"+ar_psr + "_" + str(ar_mjd) + "_" + ar_tel + ".ar"

		# Setting the DMC flag to 1. In other words, doing the DM correction.
		ar.set_dispersion_measure(dmval)
		ar.dedisperse()
		# Performing different scrunching in the archive for writing out
		if not args.Tscrunch:
			ar.tscrunch(args.tscrunch)
		else:
			ar.tscrunch()
		if not args.Fscrunch:
			ar.fscrunch(ffrac)
		else:
			ar.fscrunch()
		# Writing out the DM corrected, time/frequency scrunched file.
		ar.unload(outfile)
		if not args.quiet:
			print(" done!")
		del ar
		if not quiet:
			print("The file is corrected for DM and is written out to\n"+outfile)
	# Printing the results to the file and also in the terminal
	f= open(ar_psr+"_DM_timeseries.txt",'a')
	f.write('%s %.4f %.6f %.6f %.2f %.4f %.4f %.4f %.2f %.2f %s\n' %( filename, \
			ar_mjd, dmval, dmverr, fitchisq, pre_rms, post_rms, ToA_Err, ar_centfr, \
			ar_bw, ar_tel))
	f.close()

	import time
	end = time.time()
	total = end - start
	print ('-----------------------------------------------------------------------------')
	print ('MJD\t\tDM\t\tDMerr\t\tChisq\tC_Fr\tBW\tTel')
	print ('%.6f\t%.6f\t%.6f\t%.2f\t%.1f\t%.1f\t%s' % (ar_mjd, dmval, dmverr, 
			fitchisq, ar_centfr, ar_bw, ar_tel) )
	
	print ('-----------------------------------------------------------------------------')

	print("\nThe program took %.1f seconds to finish"%total)
#-------------------------------------------------------------------------------------------#

''' Main function that performs the DM estimation '''
def DMCalc(ar, ar_nchan, ar_centfr, ar_bw, ar_psr, ar_tel, ar_mjd, model, ephemeris, pwd, ffrac, quiet, tempo2, ptoa, narch): 
	# Checks if model file is available.
	if model == None:
		sys.exit(1)
	init_dm = ar.get_dispersion_measure()
	# setting up the ToA estimation routine using the psrchive ArrivalTime()
	if not quiet:
		print("Using the ArrivalTime (pat) with PGS in Tempo2 format")
	arrtim = psrchive.ArrivalTime()
	arrtim.set_shift_estimator('PGS')
	arrtim.set_format('Tempo2')
	arrtim.set_format_flags('IPTA')
	if not quiet:
		print("Loading the template file for processing... "),
	std = model.clone()
	std.pscrunch()
	std.tscrunch()
	std_nchan = std.get_nchan()
	
	std.dedisperse()
	std.fscrunch(ffrac)
	arrtim.set_standard(std)
	if not quiet:
		print(" done!")
	ar.fscrunch(ffrac)
	ar.pscrunch()
	ar.tscrunch()
	arrtim.set_observation(ar)
	if not quiet:
		print("Finding the ToAs... "),

	# Finding the ToAs and reading it into numpy arrays
	toas = arrtim.get_toas()
	toas_filtered = [x.split()[:5] for x in toas] 
	str_filename,str_freq,str_mjd,str_toaErr,str_site = zip(*toas_filtered)
	freq = np.asarray(str_freq, dtype=np.float64)
	amjd = np.asarray(str_mjd, dtype=np.float64)
	terr = np.asarray(str_toaErr, dtype=np.float64)
	if not quiet:
		print(" done!")
		print("Removing the bad ToAs using Huber Regression... "),
	# removing the ToAs with zero errors
	condition1 = terr < 3*np.median(terr)
	freqnew = np.extract(condition1,freq)
	amjdnew = np.extract(condition1,amjd)
	terrnew = np.extract(condition1,terr)
	# writing the ToAs to a temporary file for getting the non-phase resolved ToAs using general2
	tempfile = ar_psr+"_tmp.txt"
	f = open(tempfile,"w+")
	head="FORMAT 1\n"
	f.write('%s' % head)
	for i in range(0,np.size(freqnew)):
		f.write('%s %.12f %.20f %.8f %s\n' % 
				(str_filename[0], freqnew[i], amjdnew[i], terrnew[i], str_site[0]))
	f.close()
	tmpstr="tempo2 -output general2 -f"
	tmp = os.popen(tmpstr+" %s %s -s \"1111111 {freq} {pre} {err}\n\" | grep '1111111'" %
					 (ephemeris,tempfile)).read()
	os.remove(tempfile)

	# extracting the data from general2 output
	tmp1 = tmp.split('\n')
	freqtmp = np.zeros(np.size(amjdnew))
	toastmp = np.zeros(np.size(amjdnew))
	TErrtmp = np.zeros(np.size(amjdnew))
	for i in range(np.size(amjdnew)):
		_,freqtmp[i],toastmp[i],TErrtmp[i] = (tmp1[i].split())
	TErrtmp /= 1e+6
	# importing libraries for outlier removal
	from sklearn import linear_model
	from sklearn.linear_model import HuberRegressor
	from sklearn.preprocessing import PolynomialFeatures
	from sklearn.pipeline import make_pipeline
	# changing the shape of frequency array
	freqarr = freqtmp.reshape(-1,1)
	# making a nu^2 model and fitting using Huber Regression
	toastmp *= 1e+6
	toashift = (np.min(toastmp)*-1.5)
	toastmp += toashift
	Terrtmp = TErrtmp*1e+6
	model = make_pipeline(PolynomialFeatures(2), HuberRegressor())
	model.fit(freqarr,toastmp,
			  huberregressor__sample_weight=np.ravel(1./Terrtmp))
	y_pred = model.predict(freqarr)
	residuals = toastmp - y_pred
	median = np.median(residuals)
	MAD = np.median(np.abs(residuals-np.median(residuals)))/0.6744897501960817
	# filtering the good ToAs
	condition2 = (residuals > median - 3*MAD) & (residuals < median + 3*MAD)
	freqf = np.around(np.extract(condition2,freqarr),3)
	amjdf = np.extract(condition2,amjdnew)
	toasf = np.extract(condition2,toastmp)
	terrf = np.extract(condition2,TErrtmp)
	prefit_rms = np.sqrt(np.cov(toasf, aweights=terrf))
	
	terrf *= 1e+6
	if not quiet:
		print(" done!")
	# writing out the ToAs in proper format
	if ptoa:
		if not quiet:
			print ('Writing out ToAs into a file in tempo2 format'),
		dirtoas=os.path.join(pwd,ar_psr+"_"+ar_tel+"_ToAs")
		if not os.path.exists(dirtoas):
		   os.makedirs(dirtoas)
		outfile=dirtoas+"/"+ar_psr+"_"+str(ar_mjd)+"_"+ar_tel+"_ToAs.txt"
		f = open(outfile,"w+")
		head="FORMAT 1"
		f.write('%s\n' % head)
		for i in range(0,np.size(freqf)):
			f.write('%s %.8f %.18f %.6f %s\n' % (str_filename[0], freqf[i], amjdf[i], terrf[i], str_site[0]))
		f.close()
		if not quiet:
			print("done!")

	# Fitting the ToAs with tempo2 for DM
	if not quiet:
		print("\nWriting the ToAs to a temporary file for tempo2 fitting..."),
	outfiletmp=ar_psr+"tmp_ToAs.txt"
	f = open(outfiletmp,"w+")
	head="FORMAT 1"
	f.write('%s\n' % head)
	for i in range(0,np.size(freqf)):
		f.write('%s %.8f %.18f %.6f %s\n' % (str_filename[0], freqf[i], amjdf[i], terrf[i], str_site[0]))
	f.close()
	if not quiet:
		print(" done!\n")
	# performing the fit
	dmstr=os.popen("tempo2 -f %s %s -nofit -fit dm | grep 'DM (cm^-3 pc)'| awk \'{print $5,$6}\'" 
					% (ephemeris, outfiletmp)).read()
	(dm, dmerr) = dmstr.split()
	dmval = float(dm)
	dmverr = float(dmerr)
	# doing the fit again to read the chisquare
	chisqstr=os.popen("tempo2 -f %s %s -nofit -fit dm | grep 'Fit Chisq'| awk \'{print $9}\'" 
					% (ephemeris, outfiletmp)).read()
	fitchisq = float(chisqstr)
	os.remove(outfiletmp)

	# Preparing the data for plotting the residuals, prefit and postfit
	infile = open(ephemeris,"r")
	tmpeph1 = ar_psr+'_tmpeph.eph'
	output = open(tmpeph1,"w+")
	for i, line in enumerate(infile):
		if not line.lstrip().startswith('DM'):
				if not line.lstrip().startswith('DMEPOCH'):
					output.write(line)
	infile.close()
	output.close()
	# updating the ephemeris file with measured DM
	dmline = "DM             "+str(dmval)+"\t1\t"+str(dmverr)
	dmepochline  = "DMEPOCH	       "+str(round(ar_mjd,2))
	f = open(tmpeph1,'a')
	f.write('%s\n%s\n' % (dmline, dmepochline))
	f.close()
	newarch = ar.clone()
	newarch.tscrunch()
	newarch.set_dispersion_measure(dmval)
	arrtim.set_observation(newarch)
	arrtim.set_standard(std)
	toas1 = arrtim.get_toas()
	toas1_filtered = [x.split()[:5] for x in toas1] 
	str_filename1,str_freq1,str_mjd1,str_toaErr1,str_site1 = zip(*toas1_filtered)
	freq1 = np.asarray(str_freq1, dtype=np.float64)
	amjd1 = np.asarray(str_mjd1, dtype=np.float64)
	terr1 = np.asarray(str_toaErr1, dtype=np.float64)
	freqnew1 = np.extract(condition1,freq1)
	amjdnew1 = np.extract(condition1,amjd1)
	terrnew1 = np.extract(condition1,terr1)
	tempfile1 = ar_psr+"_tmp1.txt"
	f = open(tempfile1,"w+")
	head="FORMAT 1\n"
	f.write('%s' % head)
	for i in range(0,np.size(freqnew1)):
		f.write('%s %.12f %.20f %.8f %s\n' % (str_filename1[0], freqnew1[i], amjdnew1[i], terrnew1[i], str_site1[0]))
	f.close()

	tmp2 = os.popen("tempo2 -output general2 -f %s %s -s \"1111111 {freq} {pre} {err}\n\" | grep '1111111'" 
					% (tmpeph1,tempfile1)).read()
	os.remove(tempfile1)
	os.remove(tmpeph1)
	# extracting the data from general2 output
	tmp3 = tmp2.split('\n')
	freqtmp2 = np.zeros(np.size(amjdnew1))
	toastmp2 = np.zeros(np.size(amjdnew1))
	TErrtmp2 = np.zeros(np.size(amjdnew1))
	for i in range(np.size(amjdnew1)):
		_,freqtmp2[i],toastmp2[i],TErrtmp2[i] = (tmp3[i].split())
	freqf1 = np.around(np.extract(condition2,freqtmp2),3)
	amjdf1 = np.extract(condition2,amjdnew1)
	toasf1 = np.extract(condition2,toastmp2)
	terrf1 = np.extract(condition2,TErrtmp2)
	toasf1 *= 1e+6
	postfit_rms = np.sqrt(np.cov(toasf1, aweights=terrf1))
	ar_nbin = newarch.get_nbin()
	ar_nchn = newarch.get_nchan()
	if (narch == 1):
		freq_bot = (ar.get_centre_frequency() - ar_bw/2.0)
		freq_top = (ar.get_centre_frequency() + ar_bw/2.0)
	if (narch > 1):
		if (ar_bw == 200.):
			freq_bot = 400.0
			freq_top = 1460.0
		if (ar_bw == 400.):
			freq_bot = 300.0
			freq_top = 1460.0
	# Getting the profile data for plotting
	newarch.dedisperse()
	newarch.remove_baseline()
	profdata2D = newarch.get_data()[:,0,:,:].flatten().reshape(ar_nchn,ar_nbin)
	prof = newarch.clone()
	prof.fscrunch()
	profdata1D = prof.get_data().flatten()
	profdata1D /= np.max(profdata1D)
	residDM = init_dm - dmval
	dmcurve = 4.15 * 1000. * residDM * ( (1./(np.min(freqf)/1000.)**2) - (1./(freqf/1000.)**2) )
	dmoff = np.median(toasf) - np.median(dmcurve)
	dmcurve += dmoff
	# Now does the actual plotting	
	fig = plt.figure(3, figsize=(8, 6))
	fig.subplots_adjust(hspace=0.05)
	ax0 = plt.subplot2grid((3, 8), (0,0), rowspan=2, colspan=3)
	ax1 = plt.subplot2grid((3, 8), (2,0), rowspan=1, colspan=3)
	ax2 = plt.subplot2grid((3, 8), (0,4), colspan=4)
	ax3 = plt.subplot2grid((3, 8), (1,4), colspan=4)
	ax4 = plt.subplot2grid((3, 8), (2,4), colspan=4)
	ax0.imshow((np.sqrt(profdata2D**2))**0.5, origin='lower', extent=(0,ar_nbin-1,freq_bot,freq_top), aspect='auto', cmap='hot')
	ax0.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=12)
	ax0.tick_params(axis='x', which='both', bottom=True, top=True, 
			labelbottom=False)
	ax1.plot(np.arange(ar_nbin, dtype=float),profdata1D, color='black', linewidth=0.5)
	ax1.set_xlim(0,ar_nbin-1)
	ax1.set_xlabel('Pulse Phase (bins)', fontweight='bold', fontsize=12)
	ax1.set_ylabel('Intensity', fontweight='bold', fontsize=12)
	ax2.errorbar(freqtmp, toastmp, yerr=Terrtmp,fmt='.', color='gray', label='Prefit: Unfiltered', capsize=2)
	ax2.plot(freqtmp, y_pred,'--r', label='Polynomial Fit')
	ax2.set_xlim(freq_bot, freq_top)
	ax2.grid()
	ax2.legend(loc='upper right')
	ax2.axes.xaxis.set_ticklabels([])
	ax3.yaxis.set_label_position("right")
	ax3.errorbar(freqf, toasf-np.median(toasf), terrf,fmt='.k', label='Prefit: Filtered', capsize=2)
	ax3.set_xlim(freq_bot, freq_top)
	ax3.grid()
	ax3.legend(loc='upper right')
	ax3.axes.xaxis.set_ticklabels([])
	ax3.set_ylabel(r'ToA Residuals ($\mu$s)', fontweight='bold', fontsize=12)
	ax4.errorbar(freqf1, toasf1-np.median(toasf1), terrf1, fmt='.r', label='Postfit', capsize=2)
	ax4.set_xlim(freq_bot, freq_top)
	ax4.grid()
	ax4.legend(loc='upper right')
	ax4.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=12)
	fig.suptitle('Source: PSR %s;  MJD: %.4f;  Prefit Wrms: %.2f $\mu$s; Postfit Wrms: %.2f $\mu$s\nMedian ToA Err: %.2f $\mu$s; DM: %.6f $\pm$ %.6f pc cm$^{-3}$;  Reduced $\chi^2$: %.2f' % (ar.get_source(), ar_mjd, prefit_rms, postfit_rms, np.median(terrf1), dmval, dmverr, fitchisq), fontsize=11, fontweight='bold')
	dirplot=os.path.join(pwd,ar_psr+"_"+ar_tel+"_plots")
	if not os.path.exists(dirplot):
	   os.makedirs(dirplot)
	plotfile=dirplot+"/"+ar_psr+"_"+str(ar_mjd)+"_"+str(ar_centfr)+"_"+ar_tel+"_DMfitResid.pdf"
	plt.savefig(plotfile, format='pdf')
	plt.close()
	if not quiet:
		print ('done!')
	del ar
	return(dmval, dmverr, fitchisq, prefit_rms, postfit_rms, np.median(terrf1))


''' Frequency appending the data archives '''
def freq_appendData(narch, archives, offset, b3scrunch, b5scrunch):

	for i in range(narch):
		archives[i].tscrunch()
	# GMRT specific Jump. This is not ideal, as these jumps calculated by tempo2 
	# will be dependent on the pulsar period. Default values of this jump given 
	# is from the timing of PSR J1643-1224. 
	# PS: this jump is valid for only cycle 37 dataset (or the given MJD limits).
	if (archives[0].get_telescope() == 'GMRT'):
		for i in range(narch):
			ar_mjd = archives[i].get_Integration(0).get_start_time().in_days()
			ar_frq = archives[i].get_centre_frequency()
			ar_bw  = archives[i].get_bandwidth()
			period = (archives[i].get_Integration(0).get_folding_period())
			offset = 0.670520675
			jump = (offset/period) - int(offset/period)
			if (ar_frq >= 1260. and ar_frq < 1460.):
				if (ar_mjd >=58810. and ar_mjd < 58991.):
					archives[i].rotate_phase(-jump)
	freq_append = psrchive.FrequencyAppend()
	ttfreq = archives[0].get_centre_frequency()
	if (300. < ttfreq < 500.):
		archives[0].fscrunch(b3scrunch)
	if (1160. < ttfreq < 1460.):
		archives[0].fscrunch(b5scrunch)

	freq_append.init(archives[0])
	while len(archives) > 1:
		ttfreq = archives[1].get_centre_frequency()
		if (300. < ttfreq < 500.):
			archives[1].fscrunch(b3scrunch)
		if (1160. < ttfreq < 1460.):
			archives[1].fscrunch(b5scrunch)
		
		freq_append.append(archives[0],archives[1])
		del archives[1]
	return(archives[0])

''' Frequency Appending the Templates '''
def freq_appendModel(narch, archives, offset, b3scrunch, b5scrunch):

	for i in range(narch):
		archives[i].tscrunch()
	# GMRT specific Jump. This is not ideal, as these jumps calculated by tempo2 
	# will be dependent on the pulsar period. Default values of this jump given 
	# is from the timing of PSR J1643-1224. 
	# PS: this jump is valid for only cycle 37 dataset (or the given MJD limits).
	if (archives[0].get_telescope() == 'GMRT'):
		for i in range(narch):
			ar_mjd = archives[i].get_Integration(0).get_start_time().in_days()
			ar_frq = archives[i].get_centre_frequency()
			ar_bw  = archives[i].get_bandwidth()
			period = (archives[i].get_Integration(0).get_folding_period())
			offset = 0.670520675
			jump = (offset/period) - int(offset/period)
			if (ar_frq >= 1260. and ar_frq < 1460.):
				if (ar_mjd >=58810. and ar_mjd < 58991.):
					archives[i].rotate_phase(-jump)

	freq_append = psrchive.FrequencyAppend()
	ttfreq = archives[0].get_centre_frequency()
	if (300. < ttfreq < 500.):
		archives[0].fscrunch(b3scrunch)
	if (1160. < ttfreq < 1460.):
		archives[0].fscrunch(b5scrunch)
	freq_append.init(archives[0])
	while len(archives) > 1:
		ttfreq = archives[1].get_centre_frequency()
		if (300. < ttfreq < 500.):
			archives[1].fscrunch(b3scrunch)
		if (1160. < ttfreq < 1460.):
			archives[1].fscrunch(b5scrunch)
		freq_append.append(archives[0],archives[1])
		del archives[1]
	return(archives[0])

#----------------------------------------------------------------------------------#

main()
