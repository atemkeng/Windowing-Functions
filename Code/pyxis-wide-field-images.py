"""
Module import here
"""
try:
        import os
        import sys
        import MSResampler
        import Pyxis
        from Pyxis.ModSupport import *
        import mqt
        import pyrap.tables
	from pyrap.tables import table
        import numpy as np
        import ms
        import imager
        import scipy.special
        import pylab
	import pyfits
	import Tigger
	from scipy import signal;
except:
        print "Error: cannot load one of the modules. Check if  install."
        exit()

"""
     Conversion modules
"""
Radian2Arcmin = lambda x: x * 180 * 60 /np.pi;
Radian2deg = lambda x: x * 180./np.pi
Radian2min = lambda x:x*12*60/np.pi;
#v.DESTDIR_Template = '${OUTDIR>/}plots-${MS:BASE}${-stage<STAGE}'
def givephasecenter (msname = None):

	"""Gets  MS Dec and RA.
                Return couple: (Dec, RA) in (radian,radian)
        """;
	msname = msname or v.MS
        tab = ms.ms(msname, "FIELD").getcol("PHASE_DIR");
        dec = Radian2Arcmin (tab[0,0,1]);
        ra = Radian2min(tab[0,0,0]);
        return tab[0,0,1], tab[0,0,0];

def changePhaseCenter (msname=None, LSM=None, SUBLSM=None, DESTDIR=None):

	""" Change the phase centre of the sky model
	to feed the one of your MS
	""";
	DESTDIR = DESTDIR or "."
	msname = msname or v.MS
	dec0, ra0 = givephasecenter(msname);
	dec_d = math.floor(Radian2deg (dec0)); 
	dec_m = Radian2deg (dec0) - dec_d 
	ra_d = math.floor(Radian2deg (ra0))
	ra_m = Radian2deg (ra0) - ra_d 
	SUBLSM = SUBLSM or "%s-ra%dh%dmin0s-dec%ddeg%dmin0s.lsm.html"%(LSM.split(".")[0],ra_d,ra_m,dec_d,dec_m)
	info("******* %s phase center at ra=%dh%dmin0s, dec=%ddeg%darmin0arcsec *******"%(SUBLSM,ra_d,ra_m,dec_d,dec_m))
	command = "tigger-convert -f %s %s --recenter=%s,%dh%dm0s,%dd%dm0s"%(LSM,SUBLSM,"j2000",ra_d,ra_m,dec_d,dec_m)
	os.system(command)
	v.SUBLSM = II("${DESTDIR}/$SUBLSM")
	if os.path.exists(DESTDIR) is not None:
		makedir("$DESTDIR")
	if DESTDIR != "." :
		if os.path.exists(v.SUBLSM) is not None:
			x.mv("$SUBLSM $DESTDIR");
		else:
			x.sh("rm -fr ${DESTDIR}/$SUBLSM");	
			x.mv("$SUBLSM $DESTDIR");
	return v.SUBLSM;
    

def allSourceToSameFlux (LSM=None, SUBLSM=None, defaulflux=1., DESTDIR=None):
	"""
	From the sky model LSM, create a new sky model and put 
	all sources to default flux in  Jy, save to  DESTDIR
	"""
	SUBLSM = SUBLSM or "%s-allsrc%.3fJy.lsm.html"%(LSM.split(".")[0],defaulflux)
	DESTDIR = DESTDIR or "."
	model = Tigger.load(LSM)
	for src in model.sources:
		src.flux.I = defaulflux;
	model.save(SUBLSM);
	v.SUBLSM = II("${DESTDIR}/$SUBLSM")
	if os.path.exists(DESTDIR) is not None:
		makedir("$DESTDIR")
	if DESTDIR != "." :
		if os.path.exists(v.SUBLSM) is not None:
			x.mv("$SUBLSM $DESTDIR");
		else:
			x.sh("rm -fr ${DESTDIR}/$SUBLSM");	
			x.mv("$SUBLSM $DESTDIR");
	return v.SUBLSM;

def selectFieldBrightestSRC (LSM = None, SUBLSM = None, NUMSRC = None,DESTDIR=None):
	"""
	From the sky model LSM, create a new sky model that contain numbersrc brightest sources of LSM
	save to  DESTDIR/saveslsmname
	"""
	SUBLSM = SUBLSM or "%ibrightestsrc-%.3fjy-%s"%(NUMSRC,LSM)
	DESTDIR = DESTDIR or "."
	### complete this after

	v.SUBLSM = II("${DESTDIR}/$SUBLSM")
	if os.path.exists(DESTDIR) is not None:
		makedir("$DESTDIR")
	model.save(SUBLSM);
	x.mv("$SUBLSM $DESTDIR");
	return v.SUBLSM;
	

def selectModelByMaxFlux (LSM=None, SUBLSM = None, MAXFLUX=None,DESTDIR=None):
	""""
	this function selects a subset of sources less 
	than MAXFLUX from the sky model LSM an save the new model to DESTDIR,
   	""";
	SUBLSM = SUBLSM or "%s-%.3fmaxflux.lsm.html"%(LSM.split(".")[0],MAXFLUX)
	DESTDIR = DESTDIR or "."
	fluxlessthan = "I.le.%f"%MAXFLUX
	command = "tigger-convert -f %s %s --select=%s"%(LSM,SUBLSM,fluxlessthan)
	os.system(command)
	v.SUBLSM = II("${DESTDIR}/$SUBLSM")
	if os.path.exists(DESTDIR) is not None:
		makedir("$DESTDIR")
	if DESTDIR != "." :
		if os.path.exists(v.SUBLSM) is not None:
			x.mv("$SUBLSM $DESTDIR");
		else:
			x.sh("rm -fr ${DESTDIR}/$SUBSLM");	
			x.mv("$SUBLSM $DESTDIR");
	return v.SUBLSM;

def selectModelByMinFlux (LSM=None, SUBLSM = None, MINFLUX=None,DESTDIR=None):
	""""
	this function selects a subset of sources greater 
	than MINFLUX from the sky model LSM an save the new model to DESTDIR,
   	""";
	SUBLSM = SUBLSM or "%s-%.3fminflux.lsm.html"%(LSM.split(".")[0],MINFLUX)
	DESTDIR = DESTDIR or "."
	fluxgreater = "I.ge.%f"%MINFLUX
	command = "tigger-convert -f %s %s --select=%s"%(LSM,SUBLSM,fluxgreater)
	os.system(command)
	v.SUBLSM = II("${DESTDIR}/$SUBLSM")
	if os.path.exists(DESTDIR) is not None:
		makedir("$DESTDIR")
	if DESTDIR != "." :
		if os.path.exists(v.SUBLSM) is not None:
			x.mv("$SUBLSM $DESTDIR");
		else:
			x.sh("rm -fr ${DESTDIR}/$SUBSLM");	
			x.mv("$SUBLSM $DESTDIR");
	return v.SUBLSM;

def selectModelByRadius (LSM=None, SUBLSM = None, RADIUS=None, DESTDIR=None):
	""""
	this function selects a subset of sources in the area 2*RADIUS of the sky model 
	and save the new model to DESTDIR
   	""";
       
	SUBLSM = SUBLSM or "%s-%.2fdeg.lsm.html"%(LSM.split(".")[0],RADIUS)
	DESTDIR = DESTDIR or "."
	radius = "r.lt.%fd"%RADIUS
	command = "tigger-convert -f %s %s --select=%s"%(LSM,SUBLSM,radius)
	os.system(command)
	v.SUBLSM = II("${DESTDIR}/$SUBLSM")
	if os.path.exists(DESTDIR) is not None:
		makedir("$DESTDIR")
	if DESTDIR != "." :
		if os.path.exists(v.SUBLSM) is not None:
			x.mv("$SUBLSM $DESTDIR");
		else:
			x.sh("rm -fr ${DESTDIR}/$SUBSLM");	
			x.mv("$SUBLSM $DESTDIR");
	return v.SUBLSM;


def antennasDiameter(msname = None):
	"""
	this function gives the diameter in meter of the interferometer antenna
	""";
	msname = msname or v.MS;
        tab = table(msname+"/ANTENNA")
	tab = tab.getcol("DISH_DIAMETER");
	#info("***** Antenna size %f meter*******"%tab[0])
        return tab[0];
	
def simModelVisibility (msname = None, LSM = None, COLUMN = None, BEAMFACTOR = None, DISHDIA = None):
	"""
	The function simulate the visibilities of a 
	      sky model LSM, and put them to the specify MS column 
	"""
	COLUMN = COLUMN or "DATA"
	msname = msname or v.MS
	DISHDIA = DISHDIA or antennasDiameter(msname = v.MS)
	BEAMFACTOR = BEAMFACTOR or 1.;
	try:
		DISHDIA = DISHDIA or antennasDiameter(msname);
		options = {};
		# options['ms_sel.msname'] = msname;
		# options['ms_sel.output_column'] = COLUMN;turbo-sim:modelVisibility
		# options['tiggerlsm.filename'] = LSM;
		# options['analytic_beams.circular_aperture_beam.bf'] = BEAMFACTOR;
		# options['analytic_beams.circular_aperture_beam.dish_sizes'] = DISHDIA;
		options['ms_sel.msname'] = msname
		options['ms_sel.output_column'] = COLUMN
		options['tiggerlsm.filename']= LSM
		options['analytic_beams.wsrt_cos3_beam.dish_sizes'] = DISHDIA
		mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS", config="tdlconf.profiles", section="turbo-sim:noncircularbeam", options=options);
	except:
		print "Error: the MS and the LSM must be given" 
		exit()


def makems (conf=None,name=None, destdir = None, Ntime=8,
	                                        integration=60,dec=-45,ra=0,startfq=1400,nchan=1,chanwidth=10,starttime=None):
	"""Makes an MS using the specified parameters.
	Ntime: total synthesis duration
	integration: timeslot in seconds
	dec: declination, in degrees
	freq0: starting freq in Hz
	nchan: number of channels
	chanwidth: channel width, in Hz
	"""
	destdir = destdir or "TEMPLE"
	conf,name,destdir = interpolate_locals("conf name destdir");
	for anttab in conf,"${conf}_ANTENNA","Layouts/${conf}_ANTENNA":
		if exists(anttab):
			anttab = II(anttab);
			break;
	else:
		abort("configuration $conf not found");
	if not name:
		name = os.path.basename(anttab);
		if name.endswith("_ANTENNAS"):
			name = name.rsplit("_",1)[0];
	msname = "%s-%.1fs-%.1fMHz.MS"%(name,integration,chanwidth*1e-6);
	info("ms $msname, configuration $conf, antenna table $anttab");
	if exists(msname):
		x.sh("rm -fr $msname");
	conffile = II("makems.${msname:BASE}.cfg");
	file(conffile,"w").write(II("""
WriteAutoCorr=Tlo
StartFreq=%g
StepFreq=%g
NFrequencies=$nchan
WriteImagerColumns=T
StepTime=$integration
#TileSizeRest=1 0
NParts=1
MSDesPath=.
AntennaTableName=$anttab
Declination=$dec
NBands=1
RightAscension=$ra
StartTime=$starttime
MSName=$msname
NTimes=$Ntime 
 #TileSizeFreq=16
"""%(startfq,chanwidth)));
	info("""creating $msname: ${Ntime} timslots, ${integration}s integration, Dec=$dec, Ra=$ra
                                                          $nchan channels of $chanwidth Hz starting at $startfreq Hz""");
        # run makems
	x.makems(conffile);
	if exists(msname+"_p0") and not exists(msname):
		x.mv("${msname}_p0 $msname");
	if os.path.exists(destdir) is not None:
		makedir("$destdir")
	try:
		x.mv ("$msname $destdir")
	except:
                x.sh("rm -fr $destdir/$msname");
		x.mv ("$msname $destdir")
		info("overight and saved $msname")
	v.MS = II("${destdir}/$msname");
	return v.MS

def mscriticalsampling (msname=None, DISHDIA = None):
	"""
	Work out MS critical sampling time of the Longest baseline;
	This should be the low resolution MS maximun integration
	""";
	msname = msname or v.MS;
	DISHDIA = DISHDIA or antennasDiameter(msname)
	tab = ms.ms(msname);
	times = sorted(set(tab.getcol("TIME")));
	dt = times[1]-times[0];
	tottime = times[-1]-times[0]+dt;
	uvw = tab.getcol("UVW");
	maxbl = math.sqrt((uvw[:,:2]**2).sum(1).max());
	# max baseline sweeps out a circle of length pi*d over 24h
	arclen = math.pi*maxbl*(tottime/(24.*3600));
	# critical sampling is at half the dish size
	nsamp = arclen/(DISHDIA/2)
	# corresponding sampling interval
	critint = tottime/nsamp;
	info("(NB: critical sampling interval for this baseline is %.3fs, diam=%f, *****PLEASE USED THIS INTEGRATION FOR THE LOW RESOLUTION MS)"%(critint,DISHDIA))
	return critint;

def msinfo (msname=None):
	"""Gets various MS info.
	Return values: NUM_TIMESLOTS,NUM_CHANNELS,
	MAX_BASELINE,TOTAL_TIME,TIMESLOT_SIZE,FIRST_TIMESTAMP 
	""";
	msname = msname or v.MS;
	tab = ms.ms(msname);
	chanwidths = ms.ms(msname,"SPECTRAL_WINDOW").getcol("CHAN_WIDTH",0,1);
	nchan = chanwidths.size;
	times = sorted(set(tab.getcol("TIME")));
	dt = times[1]-times[0];
	tottime = times[-1]-times[0]+dt;
	ntime = len(times);
	uvw = tab.getcol("UVW");
	maxbl = math.sqrt((uvw[:,:2]**2).sum(1).max());
	return ntime,nchan,maxbl,tottime,dt,chanwidths[0,0];
	
v.FOV = 4.
#def windowing(funcion=None, X=None, Y=None, imagePlaneMainLobeWidth = None, fov = FOV):
#	imagePlaneMainLobeWidth = imagePlaneMainLobeWidth or fov;
#OUTFILE_Template = "${DESTDIR>/}${MS:BASE}.fits"
#LOG_Template = "log-${MS:BASE}.txt"
def sinc2 (x, y):
	"""
	Two dimensional sinc window
	""";
	deg=(v.FOV*np.pi**2)/(180.)
	x1,y1 = x*deg,y*deg; wx,wy = (np.sin(x1)/x1),(np.sin(y1)/y1);
	wx[x1==0] = 1; wy[y1==0] = 1; 
	return wx*wy;
def sincsinc2 (x, y):
	"""
	Two dimensional sinc window
	""";
	deg=(v.FOV*np.pi**2)/(180.)
	r = np.sqrt(x**2+y**2)
	x1 = r*deg; wx = (np.sin(x1)/x1);
	wx[x1==0] = 1; 
	return wx

def airy2 (x, y):
	"""
	two dimension bessel window
	""";
	a=1.*np.pi*v.FOV*np.pi/180
	r = np.sqrt(x**2+y**2); w = scipy.special.j1(r*a)/(r*a);
	w[r==0] = 0.5; return w;
def butterw2(x, y):
	"""
	two dimensional Butterwordth window
	""";
	a=2*np.pi*v.FOV*np.pi/180
	r = np.sqrt(x**2+y**2); w = scipy.special.j1(r*a)/(r*a);
	w[r==0] = 0.5; return w;

def sinchamming2(x, y):
	"""
	two dimensional sinc*hamming; this function is helpfull when you 
	want to smooth your passband; but the disavantage is that the 
	transition band width become wider that the one of sinc
	"""
	deg=(v.FOV*np.pi**2)/(180.)
	x1,y1 = x*deg,y*deg;
	sx,sy = (np.sin(x1)/x1),(np.sin(y1)/y1);
	sx[x1==0] = 1; sy[y1==0] = 1;
	hx = signal.hamming(len(sx))*np.max(sx)
	hy = signal.hamming(len(sy))*np.max(sy)
	hz = hx[:,np.newaxis]*hy[np.newaxis,:];
	wz = sx*sy
	return wz*hz;

def sincblackman2(x, y):
	"""
	two dimensional sinc*hamming; this function is helpfull when you 
	want to smooth your passband; but the disavantage is that the 
	transition band width become wider that the one of sinc
	"""
	deg=(v.FOV*np.pi**2)/(180.)
	x1,y1 = x*deg,y*deg;
	sx,sy = (np.sin(x1)/x1),(np.sin(y1)/y1);
	sx[x1==0] = 1; sy[y1==0] = 1;
	hx = signal.blackman(len(sx))*np.max(sx)
	hy = signal.blackman(len(sy))*np.max(sy)
	hz = hx[:,np.newaxis]*hy[np.newaxis,:];
	wz = sx*sy
	return wz*hz;

def sinchanning2(x, y):
	"""
	two dimensional sinc*hamming; this function is helpfull when you 
	want to smooth your passband; but the disavantage is that the 
	transition band width become wider that the one of sinc
	"""
	deg=(v.FOV*np.pi**2)/(180.)
	x1,y1 = x*deg,y*deg;
	sx,sy = (np.sin(x1)/x1),(np.sin(y1)/y1);
	sx[x1==0] = 1; sy[y1==0] = 1;
	hx = signal.hanning(len(sx))*np.max(sx)
	hy = signal.hanning(len(sy))*np.max(sy)
	hz = hx[:,np.newaxis]*hy[np.newaxis,:];
	wz = sx*sy
	return wz*hz;

def besselhamming2(x, y):
	"""
	two dimensional sinc*hamming; this function is helpfull when you 
	want to smooth your passband; but the disavantage is that the 
	transition band width become wider that the one of sinc
	"""
	a=2*np.pi*v.FOV*np.pi/180
	r = np.sqrt(x**2+y**2); 
	wb = scipy.special.j1(r*a)/(r*a);
	r = 2
	w[r==0] = 0.5; return w;
	
	deg=(v.FOV*np.pi**2)/(180.)
	x1,y1 = x*deg,y*deg;
	sx,sy = (np.sin(x1)/x1),(np.sin(y1)/y1);
	hx= signal.hamming(len(sx))*np.max(sx);
	hy = signal.hamming(len(sy))*np.max(sy);
	wx,wy = np.zeros_like(sx),np.zeros_like(sy)
	wx[:,-1] = sx[:,-1]*hx;
	wy[-1,:] = sy[-1,:]*hy;
	wx[x1==0] = 1; wy[y1==0] = 1; 
	return wx*wy;

# def bestpassband_sampling_rate (msname = None, time0 = None, ntime = None, freq0 = None, nfreq = None, window_function = sinc2):
# 	"""
# 	Work out the Fourier Transform best passband of Windowing function
# 	for an integration equal to the critical sampling time of the Longest baseline;
# 	the beam size of the interferometer at the first null, in degree
# 	""";
	
# 	msname = msname or v.MS;
# 	#DISHDIA = DISHDIA or antennasDiameter(msname)
# 	tab = ms.ms(msname);
# 	uvw = tab.getcol("UVW");
# 	# take the next integer that follow the real x, this mean that if the sampling rate,sr is 9.0<sr<=9.999  dtime = 10
# 	dtime  = 100;# int(math.ceil(mscriticalsampling (msname= msname)//1.5));
# 	dfreq = 50;
# 	p = 17;
# 	q = 26;
# 	ntime1 = ntime/dtime;
# 	nfreq1 = nfreq/dfreq;
# 	A0 = tab.getcol("ANTENNA1")
# 	A1 = tab.getcol("ANTENNA2")
# 	input_index = (A0==p)&(A1==q)
# 	uv = uvw[input_index,:2].copy();
# 	uv = uv[time0:time0+dtime,:2];
# 	if dtime%2:  # odd dtime: we have a uv-bin at the centre, so start at n/2, then take every n-th uv point
#             uv0 = uv[dtime/2::dtime,:];
# 	else: # even dtime -- take the middle of the two centre bins
#             uv0 = (uv[dtime/2-1::dtime,:] + uv[dtime/2::dtime,:])/2;
# 	info ("*********%s*********"%str(uv.shape))
# 	info ("*********%s*********"%str(uv0[0,:]))
# 	uvd1 = np.sqrt(((uv - uv0)**2).sum(1))
# 	info ("*********%s*********"%str(uvd1.shape))
# 	t2 = table(tab.getkeyword("SPECTRAL_WINDOW"),readonly=False)
# 	channels = t2.getcol("CHAN_FREQ",0);
# 	freqs = channels[0,freq0:freq0+dfreq];
# 	wavel = 3e8/freqs;
# 	wl = wavel.copy()
# 	# wl0: wavelength of centre channel of each frequency bin, shape NF1
# 	if dfreq%2:
# 		wl0 = wavel[dfreq/2::dfreq];
# 	else:
# 		wl0 = (wavel[dfreq/2-1::dfreq] + wavel[dfreq/2::dfreq])/2;
# 	info ("*********%s*********"%str(wl0.shape))
# 	uv0r = np.sqrt((uv0**2).sum(1));  
# 	# work out second component of uv-distance, shape will be NT1,NF1,DF
# 	uvd1 = uvd1 / wl0
# 	uvd2 = uv0r/wl - uv0r/wl0
# 	info ("*********%s*********"%str(wavel.shape))
# 	info ("*********%s*********"%str(uv0r.shape))
# 	info ("*********%s*********"%str(uvd2.shape))
# 	# evaluate windowing function
# 	# shapes are: uvd1 is NT1,DT,DF and uvd2 is NT1,NF1,DF, so insert new axes appropriately
# 	# wf shape is then NT1,DT,NF1,DF
# 	info ("******time***%s*********"%str(uvd1[:,np.newaxis]))
# 	info ("***freq******%s*********"%str(uvd2[np.newaxis,:]))
# 	wf = window_function(uvd1[:,np.newaxis],uvd2[np.newaxis,:]);
# 	info ("the shape of the window %s"%str(wf))
# 	import scipy.signal as signal
# 	w,h = signal.freqz(wf[:,0],np.sum(wf[:,0]))
# 	from scipy import signal;
# 	ham = butter(np.arange(len(wf[:,0])),p=1.,deg=1.)*np.max(wf[:,0])#signal.blackman
# 	#pylab.plot(wf[:,0],label="sinc")
# 	#pylab.plot(ham,label="haming")
# 	multi = ham*wf[:,0]
# 	sc = np.convolve(wf[:,0],wf[:,0])#*np.max(wf[:,0])
# 	#pylab.plot(-wf[:,0],label='sinc')
# 	pylab.plot(sc,label="sincconsinc")
# 	wl,hmul = signal.freqz(multi,np.sum(multi))
# 	h_dBm = 10*signal.log10(abs(sc))
# 	h_dB = 10*signal.log10(abs(h))
# 	#pylab.plot(h_dB,label="sinc")
# 	#pylab.plot(h_dBm,label="multi")
# 	pylab.legend()
# 	pylab.show()
# 	stop
# 	w1,h1 = signal.freqz(wf[0,:],np.sum(wf[0,:]))
# 	info ("***dtime******%i*********"%dtime)
# 	h_dB = 10*signal.log10(abs(h))
# 	h_dB1 = 10*signal.log10(abs(h1))
# 	#info ("***s******%s*********"%str((s.real).shape))
# 	pylab.plot(h_dB,label="time")
# 	#pylab.plot(h_dB1,label="freq")
# 	pylab.legend()
# 	#pylab.plot(wf[:,0])
# 	pylab.show()

def run_window_function (hiresms=None,lowresms=None,inputcolumn=None,outputcolumn=None,dtime=None,time0=0,\
				 ntime=None,dfreq=None,freq0=0,nfreq=None,window=None,overlaptime=None,\
				 overlapfreq=None):
	
	mshi =MSResampler.MSResampler(hiresms+"/",column=inputcolumn,time0=time0,ntime=ntime,freq0=freq0,nfreq=nfreq)
	arrays=mshi.overlap_window (window,dtime,dfreq,overlap_time=overlaptime,overlap_freq=overlapfreq,dump=None,dumpfig=None)
	MSResampler.save_visibility_arrays (lowresms,arrays,column=outputcolumn)

def imaging (msname=None,restore=None,dirty=None,weightsheme='natural',npix=None,cellsize=None,
				 niter=None,wprojplanes=0,gain=.1,threshold="1mJy", CLEAN_ALGORITHM="csclean",stokes="I", imagingcolumn="CORRECTED_DATA"):
       	imager.npix = npix;
	imager.cellsize = cellsize;
	#imager.stokes   = stokes;
	imager.wprojplanes = wprojplanes;
	imager.niter    =  niter;
       	imager.weight = weightsheme;
	#imager.gain     = gain;
	#imager.threshold = threshold;
	imager.CLEAN_ALGORITHM = CLEAN_ALGORITHM;
       	imager.make_image(msname=msname,column=imagingcolumn,restore = restore,dirty = dirty);

def simWithdefauldsMS(outputdir=None,hiresms=None,lowresms=None,skymodel=None,radius=0,defauld_brightes=False,\
				   src_grthan=0, src_lsthan=0,all_srcto=0, feed_modelmsphasecenter=None,\
				   ovlaptime=None,ovlapfreq=None,weightsheme="natural",window=None,\
				   fov=None,optimalw=None,crs=None,intput_column=None,output_column=None,\
				   dtime=None,dfreq=None,time0=0,ntime=None,freq0=0,nfreq=None,\
				   freqdirect=None,timedirect=None,optimal=None,dish_sizes=None,b_f=1.,\
				   npix=None,cellsize=None,niter=None,wprojection=None,\
				   gain=.1,threshold="1mJy",algorithm="csclean",stokes="I",restore=None,dirty=None):
	
	lowresms = lowresms or v.MS;
	hnum_timeslots,hnum_channels,maxbaseline,htotaltime,hdt,hdf = msinfo(msname=hiresms)
	lnum_timeslots,lnum_channels,maxbaseline,ltotaltime,ldt,ldf = msinfo(msname=lowresms)
	ntime = ntime  or  hnum_timeslots - time0; nfreq = nfreq or hnum_channels - freq0; weightsheme = weightsheme or "natural";
	if crs:
                ## get the critical sampling rate and evaluate the number of time bins to be average
		critisr =  mscriticalsampling(hiresms); dtime =  int((mscriticalsampling(hiresms))//hdt) #hintegration
	        ## get the low rest integration time, number of timeslots and the hires number of time slots to use in averaging
		info("********************here****************************%f"%critisr)
		Ntimel = (ntime/dtime); ntime = Ntimel*dtime
		if Ntimel != lnum_timeslots :
			info("***the lowres number of timeslots %i <>  %i,the number of timeslots require for a critical time of %.2fs ***"\
				     %(lnum_timeslots,Ntimel,critisr))
			exit();     
				     
	elif dtime == 0:  ## check lores save capacities
		dtime = int(ldt//hdt);
	if dfreq==0:
		dfreq = int(ldf//hdf);
	info("*******************dtime=%i, dfreq=%i***************"%(dtime,dfreq))
	if ntime%dtime !=0 or nfreq%dfreq !=0:
		info("The number of time bin (ntime=%i) and frequency bins (nfreq=%i) have to be multiple of dtime=%i and dfreq=%i respectively"\
			     %(ntime,nfreq,dtime,dfreq));
		exit();
	if os.path.exists(outputdir) is not None:
		makedir("$outputdir")
	if feed_modelmsphasecenter:
		skymodel = changePhaseCenter(msname=hiresms,LSM=skymodel,DESTDIR=outputdir)
	if radius != 0:
		subsetmodel = "%s-%.2fradius.lsm.html"%(skymodel.split('/')[1][0:-9],radius)
		skymodel_delete = skymodel
		skymodel = selectModelByRadius(LSM=skymodel,SUBLSM=subsetmodel,RADIUS=radius,DESTDIR=outputdir)
		x.sh("rm -fr $skymodel_delete");
	if defauld_brightes:
		info("not yet implement!!!!");
		exit;
	if src_grthan != 0:
		subsetmodel = "%s-%.2fminflux.lsm.html"%(skymodel.split('/')[1][0:-9],src_grthan)
		skymodel_delete = skymodel
		skymodel = selectModelByMinFlux (LSM=skymodel,SUBLSM=subsetmodel,MINFLUX=src_grthan,DESTDIR=outputdir)
		x.sh("rm -fr $skymodel_delete");
	if src_lsthan != 0:
		subsetmodel = "%s-%.2fmaxflux.lsm.html"%(skymodel.split('/')[1][0:-9],src_lsthan)
		skymodel_delete = skymodel
		skymodel = selectModelByMaxFlux (LSM=skymodel,SUBLSM=subsetmodel,MAXFLUX=src_lsthan,DESTDIR=outputdir)
		x.sh("rm -fr $skymodel_delete");
	if all_srcto !=0:
		skymodel_delete = skymodel
		subsetmodel = "%s-%.2fallsource.lsm.html"%(skymodel.split('/')[1][0:-9],src_lsthan)
		skymodel = allSourceToSameFlux (LSM=skymodel,SUBLSM=subsetmodel,defaulflux=all_srcto,DESTDIR=outputdir)
		x.sh("rm -fr $skymodel_delete");
	simModelVisibility (msname=hiresms,LSM=skymodel,COLUMN=intput_column,BEAMFACTOR=b_f,DISHDIA=dish_sizes)
	run_window_function (hiresms=hiresms,lowresms=lowresms,inputcolumn=intput_column,outputcolumn=output_column,dtime=dtime,time0=time0,\
				 ntime=ntime,dfreq=dfreq,freq0=freq0,nfreq=nfreq,window=window,overlaptime=ovlaptime,\
				 overlapfreq=ovlaptime)
	imaging(msname=lowresms,restore=restore,dirty=dirty,weightsheme=weightsheme,npix=npix,cellsize=cellsize,\
       			 niter=niter,wprojplanes=wprojection,gain=gain,threshold=threshold, CLEAN_ALGORITHM=algorithm,stokes=stokes);

def simWithmakingMS(skymodel=None,radius=0,defauld_brightes=False,src_grthan=0,\
			          src_lsthan=0,all_srcto=0, feed_modelmsphasecenter=None,\
			          ovlaptime=None,ovlapfreq=None,weightsheme=None,window=None,\
				  fov=None,imageplanerespone=None,optimalw=None,crs=None,\
			          intput_column=None,output_column=None,dtime=None,\
				  dfreq=None,time0=0,ntime=None,freq0=0,nfreq=None,use=None,\
				  freqdirect=None,timedirect=None,optimal=None,dish_sizes=None,\
				  b_f=None,confh=None,nameh=None,outputdir=None,Ntimeh=None,\
				  integrationh=None,dech=None,rah=None,startfreqh=None,nchanh=None,\
				  chanwidthh=None,starttimeh=None,confl=None,namel=None,\
				  startfreql=None, starttimel=None,chanwidthl=None,nchanl=None,\
			          Ntimel=None,integrationl=None,npix=None,cellsize=None,niter=None,wprojection=None,\
			          gain=.1,threshold="1mJy",algorithm="csclean",stokes="I",restore=None,dirty=None):

	
	## actual number of timeslot to use given by time0 and ntime, natural weighting of not precised
	ntime = ntime  or Ntimeh - time0; nfreq = nfreq or nchanh - freq0; weightsheme = weightsheme or "natural";
	intput_column = intput_column or "DATA"; output_column = output_column or "CORRECTED_DATA";
	ovlaptime = ovlaptime or 0; ovlapfreq = ovlapfreq or 0 ; fov = fov or 0.
	print "ovlatim=%f, overlapfr=%f,intput_column=%s,output_column=%s"%(ovlaptime,ovlapfreq,intput_column,output_column)
	## make hight resolution measurment set with the provide parameters;
	hiresms = makems (conf=confh,name=nameh,destdir=outputdir,Ntime=Ntimeh,integration=integrationh,
	 		dec=dech,ra=rah,startfq=startfreqh,nchan=nchanh,chanwidth=chanwidthh,starttime=starttimeh);
	## work out the parameters of the low resolution ms;
	startfreql = startfreql or startfreqh;
	##evaluate the lowres integration and the number of timeslots if using 
	##hires  critical sample rate of the longest baseline or not
	if crs:
		try:
			## get the critical sampling rate and evaluate the number of time bins to be average
			critisr =  mscriticalsampling(hiresms); dtime =  int((mscriticalsampling(hiresms))//integrationh)
	                ## get the low rest integration time, number of timeslots and the hires number of time slots to use in averaging
			integrationl = dtime*integrationh; Ntimel = (ntime/dtime)*ntime/ntime; ntime = Ntimel*dtime 
			info("***critical sampling rate=%d, number of bin averaged=%d***"%(integrationl,dtime))
		except:
			print "**Error: hires integration time is greater or approximate the critical sampling time=%f \
                                  or your synthesis time haven't reach the telescope critical sampling time**"%(critisr)
			exit()
	else:  
		integrationl =integrationl or dtime*integrationh;
		dtime = dtime or int(integrationl//integrationh);
		Ntimel = Ntimel or ntime//dtime;
		ntime = Ntimel*dtime;
		
	chanwidthl = chanwidthl or dfreq*chanwidthh; 
	dfreq = dfreq or int(chanwidthl//chanwidthh);
	nchanl = nchanl or nfreq//dfreq; 
	nfreq = nchanl*dfreq;
	info("****dt=%i,df=%i,ntime=%i,nfreq=%i*******"%(dtime,dfreq,ntime,nfreq))
        lowresms = makems (conf=confl,name=namel,destdir=outputdir,Ntime=Ntimel,integration=integrationl,\
	  		dec=dech,ra=rah,startfq=startfreql,nchan=nchanl,chanwidth=chanwidthl,starttime=starttimel);
	## if the dish size is not given the take the dish diameter of the ms with beam factor 1
	dish_sizes = dish_sizes or antennasDiameter(lowresms); b_f = b_f or 1.;
	## if the FOV is not given then take the Primary Beam size at the FWHM
	FOV = fov if fov > 0 else Radian2deg(2*1.2*(3e8/startfreqh)/dish_sizes);
	if os.path.exists(outputdir) is not None:
		makedir("$outputdir")
	if feed_modelmsphasecenter:
		skymodel = changePhaseCenter(msname=hiresms,LSM=skymodel,DESTDIR=outputdir)
	if radius != 0:
		subsetmodel = "%s-%.2fradius.lsm.html"%(skymodel.split('/')[1][0:-9],radius)
		skymodel_delete = skymodel
		skymodel = selectModelByRadius(LSM=skymodel,SUBLSM=subsetmodel,RADIUS=radius,DESTDIR=outputdir)
		x.sh("rm -fr $skymodel_delete");
	if defauld_brightes:
		info("not yet implement!!!!");
		exit;
	if src_grthan != 0:
		subsetmodel = "%s-%.2fminflux.lsm.html"%(skymodel.split('/')[1][0:-9],src_grthan)
		skymodel_delete = skymodel
		skymodel = selectModelByMinFlux (LSM=skymodel,SUBLSM=subsetmodel,MINFLUX=src_grthan,DESTDIR=outputdir)
		x.sh("rm -fr $skymodel_delete");
	if src_lsthan != 0:
		subsetmodel = "%s-%.2fmaxflux.lsm.html"%(skymodel.split('/')[1][0:-9],src_lsthan)
		skymodel_delete = skymodel
		skymodel = selectModelByMaxFlux (LSM=skymodel,SUBLSM=subsetmodel,MAXFLUX=src_lsthan,DESTDIR=outputdir)
		x.sh("rm -fr $skymodel_delete");
	if all_srcto !=0:
		skymodel_delete = skymodel
		subsetmodel = "%s-%.2fallsource.lsm.html"%(skymodel.split('/')[1][0:-9],src_lsthan)
		skymodel = allSourceToSameFlux (LSM=skymodel,SUBLSM=subsetmodel,defaulflux=all_srcto,DESTDIR=outputdir)
		x.sh("rm -fr $skymodel_delete");
	simModelVisibility (msname=hiresms,LSM=skymodel,COLUMN=intput_column,BEAMFACTOR=b_f,DISHDIA=dish_sizes)
	run_window_function (hiresms=hiresms,lowresms=lowresms,inputcolumn=intput_column,outputcolumn=output_column,dtime=dtime,time0=time0,\
				 ntime=ntime,dfreq=dfreq,freq0=freq0,nfreq=nfreq,window=window,overlaptime=ovlaptime,\
				 overlapfreq=ovlaptime)
	imaging(msname=lowresms,restore=restore,dirty=dirty,weightsheme=weightsheme,npix=npix,cellsize=cellsize,\
       			 niter=niter,wprojplanes=wprojection,gain=gain,threshold=threshold, CLEAN_ALGORITHM=algorithm,stokes=stokes);
	

def src_supresion_fradius(ovlaptime=None,ovlapfreq=None,weightsheme=None,window=None,\
				  fov=0.,optimalw=None,crs=None,stm=None,stepm=None,\
				  edm=None,intput_column=None,output_column=None,dtime=None,\
				  dfreq=None,time0=0,ntime=None,freq0=0,nfreq=None,use=None,\
				  freqdirect=None,timedirect=None,optimal=None,dish_sizes=0.,\
				  b_f=None,confh=None,nameh=None,outputdir=None,Ntimeh=None,\
				  integrationh=None,dech=None,rah=None,startfreqh=None,nchanh=None,\
				  chanwidthh=None,starttimeh=None,confl=None,namel=None,\
	 			  startfreql=None, starttimel=None,integrationl=None,chanwidthl=None,\
				  Ntimel=None,nchanl=None,noise=None,ramdomseed=None,nameoutputfile=None ):
	"""
	This function     
	"""
	info("direction time=%s, freq=%s"%(timedirect,freqdirect))
	ntime = ntime  or Ntimeh - time0; nfreq = nfreq or nchanh - freq0; weightsheme = weightsheme or "natural";
	intput_column = intput_column or "DATA"; output_column = output_column or "CORRECTED_DATA";
	ovlaptime = ovlaptime or 0; ovlapfreq = ovlapfreq or 0 ;  
	print "ovlatim=%f, overlapfr=%f,intput_column=%s,output_column=%s"%(ovlaptime,ovlapfreq,intput_column,output_column)
	## make hight resolution measurment set with the provide parameters;
	hiresms = makems (conf=confh,name=nameh,destdir=outputdir,Ntime=Ntimeh,integration=integrationh,
	 		dec=dech,ra=rah,startfq=startfreqh,nchan=nchanh,chanwidth=chanwidthh,starttime=starttimeh);
	## work out the parameters of the low resolution ms;
	startfreql = startfreql or startfreqh;
        ## prepare the output list of fluxes 
        src_locat_m = np.arange(stm,edm,stepm); 
        Listflux = np.zeros(len(src_locat_m), dtype=float)
	##evaluate the lowres integration and the number of timeslots if using 
	##hires  critical sample rate of the longest baseline or not
	if crs:
		try:
			## get the critical sampling rate and evaluate the number of time bins to be average
			critisr =  mscriticalsampling(hiresms);
			dtime =  int((mscriticalsampling(hiresms))//integrationh)
	                ## get the low rest integration time, number of timeslots and the hires number of time slots to use in averaging
			integrationl = dtime*integrationh; Ntimel = (ntime/dtime)*ntime/ntime; ntime = Ntimel*dtime 
			info("***critical sampling rate=%d, number of bin averaged=%d***"%(integrationl,dtime))
		except:
			print "**Error: hires integration time is greater or approximate the critical sampling time=%f \
                                  or your synthesis time haven't reach the telescope critical sampling time**"%(critisr)
			exit()
	else:  
		integrationl =integrationl or dtime*integrationh;
		dtime = dtime or int(integrationl//integrationh);
		Ntimel = Ntimel or ntime//dtime;
		ntime = Ntimel*dtime;
	## work out the Low resolution MS channel width, number of channel
	chanwidthl = chanwidthl or dfreq*chanwidthh; 
	dfreq = dfreq or int(chanwidthl//chanwidthh);
	nchanl = nchanl or nfreq//dfreq; 
	nfreq = nchanl*dfreq;	
	info("****dt=%i,df=%i,ntime=%i,nfreq=%i*******"%(dtime,dfreq,ntime,nfreq))
        lowresms = makems (conf=confl,name=namel,destdir=outputdir,Ntime=Ntimel,integration=integrationl,\
	  		dec=dech,ra=rah,startfq=startfreql,nchan=nchanl,chanwidth=chanwidthl,starttime=starttimel);
	## if the dish size is not given the take the dish diameter of the ms with beam factor 1
	dish_sizes = dish_sizes if dish_sizes > 0. else  antennasDiameter(lowresms); b_f = b_f if b_f > 0. else 1.;
	## if the FOV is not given then take the Primary Beam size at the FWHM
	v.FOV = fov if fov > 0. else Radian2deg(2*1.2*(3e8/startfreqh)/dish_sizes);
	info("Field of view not specified, the program use a field of view at the FWHM of the primary beam (%f degree)"%(v.FOV))
	dec = dech.split("."); ra = rah.split(":"); windowing = "%s"%window;
	sign = (dec[-3][0] == '-' and -1 or 1 )
	dec = sign*(float(dec[-1])/60. + float(dec[-2])) +float(dec[-3])*60.
	for index  in range(len(src_locat_m)):
		# prepare options for Meqserver
		options = {}
		options['gridded_sky.grid_m0'] = src_locat_m[index]
		options['ms_sel.msname'] = hiresms;
		options['ms_sel.output_column'] = intput_column;
		mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",
			config="tdlconf.profiles",section="Sim_source_radius",
			options=options);
		# options = {};
		# options['ms_sel.msname'] = hiresms; 
		# #options['ms_sel.output_column'] = intput_column;
		# options['analytic_beams.circular_aperture_beam.bf'] = b_f;
		# options['analytic_beams.circular_aperture_beam.dish_sizes'] = dish_sizes;
		# options['gridded_sky.grid_m0'] = src_locat_m[index];
		# mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",config="tdlconf.profiles",
		# 	section="turbo-sim:source_radius", options=options);
		center_min = dec+src_locat_m[index];
		center_deg = math.ceil(center_min/60);
                center_min = abs(center_min - center_deg*60);
		mshi =MSResampler.MSResampler(hiresms+"/",column=intput_column,time0=time0,ntime=ntime,freq0=freq0,nfreq=nfreq)
		if window == "boxcar":
			arrays=mshi.boxcar (dtime,dfreq)
			windowing = "boxcar %s"%window;
		else:
			arrays=mshi.overlap_window (window,dtime,dfreq,overlap_time=ovlaptime,overlap_freq=ovlapfreq,dump=None,dumpfig=None)
                MSResampler.save_visibility_arrays (lowresms,arrays,column=output_column)
		imager.npix= 256;imager.cellsize = "2arcsec"
		info ("****************dec=%ddeg%dm ***************"%(center_deg,center_min))
                imager.make_image(msname=lowresms,column =output_column,phasecenter = "j2000,%dh%dm,%dd%dm" %(float(ra[0]),
		       float(ra[1]),center_deg,center_min),restore = False,dirty = True, restore_lsm = False, weight = weightsheme);
		Listflux[index] = np.max(pyfits.open(imager.DIRTY_IMAGE)[0].data[0][0])
	## pickling functionality is not loaded by default; you have to import it

	if noise != None:
		
		
		evaluate_noise = noise_estimate (noise,ramdomseed,ovlaptime,ovlapfreq,weightsheme,window,\
				  fov,optimalw,crs,intput_column,\
			          output_column, dtime,\
				  dfreq,time0,ntime,freq0,nfreq,use,\
				  freqdirect,timedirect,optimal,dish_sizes,\
				  b_f,confh,nameh,outputdir,Ntimeh,\
				  integrationh,dech,rah,startfreqh,nchanh,\
				  chanwidthh,starttimeh,confl,namel,\
	 			  startfreql, starttimel,integrationl,chanwidthl,\
				  Ntimel,nchanl);
		

	
	import pickle;
	if window == "boxcar":
		filename = "%s/%s-int%.2fs-bandwidth%.2fMHz.data"%(outputdir,windowing.split()[1], integrationl,chanwidthl*1e-6)	
	else:
		filename = "%s/%s-%ix%i-FoV%.2fdeg-int%.2fs-bandwidth%.2fMHz.data"%(outputdir,windowing.split()[1], ovlaptime,ovlapfreq, v.FOV, integrationl,chanwidthl*1e-6)
	#nameoutputfile = nameoutputfile or "data-save-here"
	#filename ="%s/%s.data"%(outputdir,nameoutputfile)
	#info("********************here noise %s****************************"%filename)

	f = open(filename, 'wb'); 
	src_locat_m = src_locat_m/60.;
	type_window = "%s-%ix%i"%(windowing.split()[1],int((dtime+2*ovlaptime)/dtime),int((dfreq+2*ovlapfreq)/dfreq))
	dict = {'window':type_window,'flux':Listflux,'noise':evaluate_noise,'radius':src_locat_m,'fov':v.FOV,'Bandwidth':chanwidthl*1e-6,'int':integrationl}
	pickle.dump(dict, f)
	f.close()
	x.sh("rm -fr $outputdir/*.MS")
	x.sh("rm -fr $outputdir/*.MS")
        x.sh("rm -rf  plotsvlac-lores-*")
	x.sh("rm makems.vlac-*")
	x.sh("rm vlac-lores-*")
        #pylab.plot(src_locat_m, Listflux,"b-o", label="%s-%i-%i, n=%f"%(windowing.split()[1],ovlaptime,ovlapfreq,noise))
	#pylab.savefig("%s/%s-%s-%iovlaptimex%iovlapfreq-fluxvsradius"%(outputdir,outputdir,windowing.split()[1],ovlaptime,ovlapfreq));



def noise_estimate (noise=None,ramdomseed=None,ovlaptime=None,ovlapfreq=None,weightsheme=None,window=None,\
				  fov=0.,optimalw=None,crs=None,intput_column=None,\
			          output_column=None, dtime=None,\
				  dfreq=None,time0=0,ntime=None,freq0=0,nfreq=None,use=None,\
				  freqdirect=None,timedirect=None,optimal=None,dish_sizes=0.,\
				  b_f=None,confh=None,nameh=None,outputdir=None,Ntimeh=None,\
				  integrationh=None,dech=None,rah=None,startfreqh=None,nchanh=None,\
				  chanwidthh=None,starttimeh=None,confl=None,namel=None,\
	 			  startfreql=None, starttimel=None,integrationl=None,chanwidthl=None,\
				  Ntimel=None,nchanl=None):
			
	ntime = ntime  or Ntimeh - time0; nfreq = nfreq or nchanh - freq0; 
	## make hight resolution measurment set with the provide parameters;
	hiresms = makems (conf=confh,name=nameh,destdir=outputdir,Ntime=Ntimeh,integration=integrationh,
	 		dec=dech,ra=rah,startfq=startfreqh,nchan=nchanh,chanwidth=chanwidthh,starttime=starttimeh);
	## work out the parameters of the low resolution ms;
	startfreql = startfreql or startfreqh;
	##evaluate the lowres integration and the number of timeslots if using 
	##hires  critical sample rate of the longest baseline or not
	if crs:
		try:
			## get the critical sampling rate and evaluate the number of time bins to be average
			critisr =  50.#mscriticalsampling(hiresms); 
			dtime =  int(50./integrationh)#int((mscriticalsampling(hiresms))//integrationh)
	                ## get the low rest integration time, number of timeslots and the hires number of time slots to use in averaging
			integrationl = dtime*integrationh; Ntimel = (ntime/dtime)*ntime/ntime; ntime = Ntimel*dtime 
			info("***critical sampling rate=%d, number of bin averaged=%d***"%(integrationl,dtime))
		except:
			print "**Error: hires integration time is greater or approximate the critical sampling time=%f \
                                  or your synthesis time haven't reach the telescope critical sampling time**"%(critisr)
			exit()
	else:  
		integrationl =integrationl or dtime*integrationh;
		dtime = dtime or int(integrationl//integrationh);
		Ntimel = Ntimel or ntime//dtime;
		ntime = Ntimel*dtime;
	## work out the Low resolution MS channel width, number of channel
	chanwidthl = chanwidthl or dfreq*chanwidthh; 
	dfreq = dfreq or int(chanwidthl//chanwidthh);
	nchanl = nchanl or nfreq//dfreq; 
	nfreq = nchanl*dfreq; windowing = "%s"%window;
	lowresms = makems (conf=confl,name=namel,destdir=outputdir,Ntime=Ntimel,integration=integrationl,\
	  		dec=dech,ra=rah,startfq=startfreql,nchan=nchanl,chanwidth=chanwidthl,starttime=starttimel);
	## if the dish size is not given the take the dish diameter of the ms with beam factor 1
	dish_sizes = dish_sizes if dish_sizes > 0. else  antennasDiameter(lowresms); b_f = b_f if b_f > 0. else 1.;
	## if the FOV is not given then take the Primary Beam size at the FWHM
	v.FOV = fov if fov > 0. else Radian2deg(2*1.2*(3e8/startfreqh)/dish_sizes);
	info("Field of view not specified, the program use a field of view at the FWHM of the primary beam (%f degree)"%(v.FOV))
	# noise simulation, prepare options for Meqserver
	info("MS========================%s"%hiresms)
	options = {}
	options['ms_sel.msname'] = hiresms;
	options['gridded_sky.source_flux'] = 0.;
	options['noise_stddev'] = noise;
	options['random_seed'] = ramdomseed;
	options['ms_sel.output_column'] = intput_column;
	mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",config="tdlconf.profiles",section="Sim_source_radius", options=options);
	#options = {};
	#options['ms_sel.msname'] = hiresms; 
	#options['ms_sel.output_column'] = intput_column;
	#options['analytic_beams.circular_aperture_beam.bf'] = b_f;
	#options['analytic_beams.circular_aperture_beam.dish_sizes'] = dish_sizes;
	#options['gridded_sky.grid_m0'] = 0;
	#options['gridded_sky.source_flux'] = 0.
	#options['noise_stddev'] = noise;
	#options['random_seed'] = ramdomseed;
	#mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",config="tdlconf.profiles",
	#	 	section="turbo-sim:source_radius", options=options);
	info("intpucol=%s,outputcol=%s,time0=%s,ntime=%s,freq0=%s,nfreq=%s"%(intput_column,output_column,str(time0),str(ntime),str(freq0),str(nfreq)))
	mshi =MSResampler.MSResampler(hiresms+"/",column=intput_column,time0=time0,ntime=ntime,freq0=freq0,nfreq=nfreq)
	if window == "boxcar":
		arrays=mshi.boxcar (dtime,dfreq)
		windowing = "boxcar %s"%window;
	else:
		arrays=mshi.overlap_window (window,dtime,dfreq,overlap_time=ovlaptime,overlap_freq=ovlapfreq,dump=None,dumpfig=None)
	MSResampler.save_visibility_arrays (lowresms,arrays,column=output_column)
	imager.npix= 256;imager.cellsize = "2arcsec"
	imager.make_image(msname=lowresms,column =output_column,restore = False,dirty = True, restore_lsm = False, weight=weightsheme);
	window_noise = np.std(pyfits.open(imager.DIRTY_IMAGE)[0].data[0][0])
	return window_noise
	
def compression_noise (noise=None,ramdomseed=None,ovlaptime=None,ovlapfreq=None,weightsheme=None,window=None,\
				  fov=0.,optimalw=None,crs=None,intput_column=None,\
			          output_column=None, dtime=None,\
				  dfreq=None,time0=0,ntime=None,freq0=0,nfreq=None,use=None,\
				  freqdirect=None,timedirect=None,optimal=None,dish_sizes=0.,\
				  b_f=None,confh=None,nameh=None,outputdir=None,Ntimeh=None,\
				  integrationh=None,dech=None,rah=None,startfreqh=None,nchanh=None,\
				  chanwidthh=None,starttimeh=None,confl=None,namel=None,\
	 			  startfreql=None, starttimel=None,integrationl=None,chanwidthl=None,\
				  Ntimel=None,nchanl=None):
			
	ntime = ntime  or Ntimeh - time0; nfreq = nfreq or nchanh - freq0; 
	## make hight resolution measurment set with the provide parameters;
	hiresms = makems (conf=confh,name=nameh,destdir=outputdir,Ntime=Ntimeh,integration=integrationh,
	 		dec=dech,ra=rah,startfq=startfreqh,nchan=nchanh,chanwidth=chanwidthh,starttime=starttimeh);
	## work out the parameters of the low resolution ms;
	startfreql = startfreql or startfreqh;
	##evaluate the lowres integration and the number of timeslots if using 
	##hires  critical sample rate of the longest baseline or not
	if crs:
		try:
			## get the critical sampling rate and evaluate the number of time bins to be average
			critisr =  50.#mscriticalsampling(hiresms); 
			dtime =  int(50./integrationh)#int((mscriticalsampling(hiresms))//integrationh)
	                ## get the low rest integration time, number of timeslots and the hires number of time slots to use in averaging
			integrationl = dtime*integrationh; Ntimel = (ntime/dtime)*ntime/ntime; ntime = Ntimel*dtime 
			info("***critical sampling rate=%d, number of bin averaged=%d***"%(integrationl,dtime))
		except:
			print "**Error: hires integration time is greater or approximate the critical sampling time=%f \
                                  or your synthesis time haven't reach the telescope critical sampling time**"%(critisr)
			exit()
	else:  
		integrationl =integrationl or dtime*integrationh;
		dtime = dtime or int(integrationl//integrationh);
		Ntimel = Ntimel or ntime//dtime;
		ntime = Ntimel*dtime;
	## work out the Low resolution MS channel width, number of channel
	chanwidthl = chanwidthl or dfreq*chanwidthh; 
	dfreq = dfreq or int(chanwidthl//chanwidthh);
	nchanl = nchanl or nfreq//dfreq; 
	nfreq = nchanl*dfreq; windowing = "%s"%window;
	lowresms = makems (conf=confl,name=namel,destdir=outputdir,Ntime=Ntimel,integration=integrationl,\
	  		dec=dech,ra=rah,startfq=startfreql,nchan=nchanl,chanwidth=chanwidthl,starttime=starttimel);
	## if the dish size is not given the take the dish diameter of the ms with beam factor 1
	dish_sizes = dish_sizes if dish_sizes > 0. else  antennasDiameter(lowresms); b_f = b_f if b_f > 0. else 1.;
	## if the FOV is not given then take the Primary Beam size at the FWHM
	v.FOV = fov if fov > 0. else Radian2deg(2*1.2*(3e8/startfreqh)/dish_sizes);
	info("Field of view not specified, the program use a field of view at the FWHM of the primary beam (%f degree)"%(v.FOV))
	# noise simulation, prepare options for Meqserver
	info("MS========================%s"%hiresms)
	options = {}
	options['ms_sel.msname'] = hiresms;
	options['gridded_sky.source_flux'] = 0.;
	options['noise_stddev'] = noise;
	options['random_seed'] = ramdomseed;
	options['ms_sel.output_column'] = intput_column;
	mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",config="tdlconf.profiles",section="Sim_source_radius", options=options);
	#options = {};
	#options['ms_sel.msname'] = hiresms; 
	#options['ms_sel.output_column'] = intput_column;
	#options['analytic_beams.circular_aperture_beam.bf'] = b_f;
	#options['analytic_beams.circular_aperture_beam.dish_sizes'] = dish_sizes;
	#options['gridded_sky.grid_m0'] = 0;
	#options['gridded_sky.source_flux'] = 0.
	#options['noise_stddev'] = noise;
	#options['random_seed'] = ramdomseed;
	#mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",config="tdlconf.profiles",
	#	 	section="turbo-sim:source_radius", options=options);
	info("intpucol=%s,outputcol=%s,time0=%s,ntime=%s,freq0=%s,nfreq=%s"%(intput_column,output_column,str(time0),str(ntime),str(freq0),str(nfreq)))
	mshi =MSResampler.MSResampler(hiresms+"/",column=intput_column,time0=time0,ntime=ntime,freq0=freq0,nfreq=nfreq)
	if window == "boxcar":
		arrays=mshi.boxcar (dtime,dfreq)
		windowing = "boxcar %s"%window;
	else:
		arrays=mshi.overlap_window (window,dtime,dfreq,overlap_time=ovlaptime,overlap_freq=ovlapfreq,dump=None,dumpfig=None)
	MSResampler.save_visibility_arrays (lowresms,arrays,column=output_column)
	imager.npix= 256;imager.cellsize = "2arcsec"
	imager.make_image(msname=lowresms,column =output_column,restore = False,dirty = True, restore_lsm = False, weight=weightsheme);
	window_noise = np.std(pyfits.open(imager.DIRTY_IMAGE)[0].data[0][0])
	memoryHires = np.array([ntime,nfreq,])
	memoryLores = np.array([Ntimel,nchanl])
	import pickle;
	type_window = "%s-%ix%i"%(windowing.split()[1],int((dtime+2*ovlaptime)/dtime),int((dfreq+2*ovlapfreq)/dfreq))
	filename = "%scompression-noise.data"%type_window
	f = open(filename, "a"); 
	
	dict = {'window':type_window,'HiresMS':memoryHires,'LowresMS':memoryLores,'Noise':window_noise}
	pickle.dump(dict, f)
	f.close()
        
	x.sh("rm -fr $outputdir/*.MS")
	x.sh("rm -fr $outputdir/*.MS")
        x.sh("rm -rf  plotsvlac-lores-*")
	x.sh("rm makems.vlac-*")
	x.sh("rm vlac-lores-*")
	
	
define("HIRES_VLAC",[100,1.5,-45,1400,1,10,"2011/11/16/17:00"],"synthesis(h), integr(s),dec(d.m.s),freq(MHz),Nchannels(int),Chanwidth(MHz),starttime(date)");
define("LORES_VLAC",[1,150,-45,1400,1,10,'2011/11/16/17:00'],"synthesis(h), integr(s),dec(d.m.s),freq(MHz),Nchannels(int),Chanwidth(MHz),starttime(date)");
define("SRC",[0.,10.,60.],"source suppression as a function of radius from 0arcmin to 2000arcmin, spacing by 10arcmin")
define("FIGURE_WIDTH",4,"with")
define("FIGURE_HEIGHT",4.,"height")
define("FIGURE_DPI",1,"dpd")
define("HIRES",[100,1.5,-45,1400,1,10,"2011/11/16/17:00"],"synthesis(h), integr(s),dec(d.m.s),freq(MHz),Nchannels(int),Chanwidth(MHz),starttime(date)");
define("FREQDIRECTION",True,'work accross the frequency direction')
define("TIMEDIRECTION",True,'work accross the time direction')
define("TF",[0,0,None,None],'(time0,freq0,ntime,nfreq) is the defauld, in this case take all frequency and time bins')
define("OPTIMAL",True,'use the optiomal window for this telescope')
define("CRISPL",True,'use the optiomal window for this telescopetelescope critical sampling rate')	
define("USE",True,'run this window')	
define("DTDF",[100,10],"number of time and frequency bins to compress")
def run_serie_function (conf=None,name=None, destdir = None, hours=None,
	                                        integration=None,dec=None,freq0=None,nchan=None,chanwidth=None,starttime=None):
	"""
	run a set of makems and simulation
	"""

	makems (conf=conf or MAKEMS_CONF,name=name or MAKEMS_NAME, destdir=destdir or MAKEMS_DIR,
		hours= hours or HIRES_VLAC[0], integration= integration or HIRES_VLAC[1],
		dec=dec or HIRES_VLAC[2],freq0=freq0 or HIRES_VLAC[3],nchan=nchan or HIRES_VLAC[4],
		chanwidth=chanwidth or HIRES_VLAC[5],starttime=starttime or HIRES_VLAC[6]);

	#simModelVisibility (LSM = "nvss-4degree-across.lsm.html")
	crit = mscriticalsampling (msname=v.MS,DISHDIA=antennasDiameter(msname=v.MS))
	info ("***********diameter disk VLAC%s critical simpling %f*************"%(antennasDiameter(msname=v.MS),crit))
	
	

	
