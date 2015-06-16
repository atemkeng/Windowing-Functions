import ConfigParser
import sys;
import os;
from ConfigParser import SafeConfigParser

def ConfigSectionMap(Config,section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

def read_hires (section,configfile):
#	try:
		return  section['antennatablename'],section['hires_telescop_name'],section['ouput_directory'],\
			int(section['num_time_steps']),float(section['step_time']),section['phase_centre_dec'],section['phase_centre_ra'],\
			float(section['start_frequency_hz']),int(section['num_frequencies']),float(section['frequency_inc_hz']),\
			section['start_time_utc'],section['name_file_output_data'];
#	except:
#        	print "Error: check  the section hires of the config file %s "%(configfile)
#        	exit()	

def read_lores(section,configfile):
	try:
        	return section['antennatablename'],section['lores_telescop_name'],float(section['start_frequency_hz']),\
		section['start_time_utc'],float(section['frequency_inc_hz']),float(section['step_time']);
        except:
        	print "Error: check  the section lores of config file %s "%(configfile)
        	exit()

def read_ptsrcsky (section,configfile):
#	try:
        	return  section['list_start_step_stop.grid_m0'],section['imaging_column'],\
                	section['simulate_column'],float(section['dish_sizes']),float(section['beam_factor']),\
			section['usecriticalsamplingrate'],int(section['time0']),int(section['freq0']),\
			int(section['ntimebins']),int(section['nfreqbins']),int(section['numbertimebintocompress']),\
			int(section['numberfreqbintocompress']),section['weightingsheme'],int(section['overlapfreq']),\
			int(section['overlaptime']),float(section['fov']);
 #       except:
  #      	print "Error: check the section %s in the config file %s "%(section,configfile)
   #     	exit()

def compres_rate_noise(section,configfile):
	return section,float(section['noise']),section['random_seed'];

def read_windows (config,section,configfile):
        #try:
        	return  section,section['use'],section['useoptimalwindow'];

        #except:
        #	print "Error: check the section %s in the config file %s "%(section,configfile)

def get_parameter_telescope(hiresms,loresms,parawind,window,cflname):
	conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,chanwidth1,starttime1,nameoutputfile = read_hires (hiresms,cflname)
	conf2,name2,startfreq2,starttime2,step_freq, integrationl = read_lores (loresms,arg[1])	

	g_m0,outputcol,intputcol,dish_s,b_f,crsampl,t0,f0,ntime,nfreq,dt,df,weightshm,ovf,ovt,fov = read_ptsrcsky(ptsrcsky,cflname) 
	g_m0 = g_m0.split(',')
	sect,noise,ramdom = compres_rate_noise(rate_noise,cflname);
	
	sect1,use1,opt1 = read_windows (Config,parawind,cflname)
	command1 = "";
        if crsampl.upper()=="FALSE":
	    print "here crp", 
            if step_freq <= 0.:
                if integrationl <= 0.:
                    if dt <= 0 or df <= 0:
                        print """Please give step_time or numbertimebintocompress 
                                  and frequency_inc_hz or  numberfreqbintocompress""";
                        sys.exit()
                    else:
                        command = "pyxis src_supresion_fradius[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,stm=%.1f,stepm=%.1f,edm=%.1f,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,dtime=%i,dfreq=%i,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,intput_column=%s,output_column=%s,noise=%f,ramdomseed=%s,nameoutputfile=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
		chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,\
		float(g_m0[0]),float(g_m0[1]),float(g_m0[2]),crsampl,t0,f0,ntime,nfreq,dt,df,window,fov,opt1,weightshm,ovt,ovf,intputcol,outputcol,noise,ramdom,nameoutputfile)
			#print "**************************************here****************
			#if rate_noise['use'].upper()=="TRUE":
			#	command1 = "pyxis compression_noise[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,dtime=%i,dfreq=%i,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,noise=%f,ramdomseed=%s,intput_column=%s,output_column=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
                #chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,crsampl,t0,f0,ntime,nfreq,dt,df,window,fov,opt1,weightshm,ovt,ovf,noise,ramdom,intputcol,outputcol)
                else:
                    if df <= 0:
                        print """ Please give frequency_inc_hz or  numberfreqbintocompress"""
                        sys.exit()
                    else:
                        command = "pyxis src_supresion_fradius[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,stm=%.1f,stepm=%.1f,edm=%.1f,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,integrationl=%f,dfreq=%i,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,intput_column=%s,output_column=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
		chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,\
		float(g_m0[0]),float(g_m0[1]),float(g_m0[2]),crsampl,t0,f0,ntime,\
		nfreq,integrationl,df,window,fov,opt1,weightshm,ovt,ovf,intputcol,outputcol)
			if  rate_noise['use'].upper()=="TRUE":
				command1 = "pyxis compression_noise[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,integrationl=%f,dfreq=%i,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,noise=%f,ramdomseed=%s,intput_column=%s,output_column=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
                chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,crsampl,t0,f0,ntime,\
                nfreq,integrationl,df,window,fov,opt1,weightshm,ovt,ovf,noise,ramdom,intputcol,outputcol)
            else:
                if integrationl <= 0:
                    if dt <= 0:
                        print """ Please give the step_time or the  numbertimebintocompress """
                        sys.exit()
                    else:
                        command = "pyxis src_supresion_fradius[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,stm=%.1f,stepm=%.1f,edm=%.1f,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,dtime=%i,chanwidthl=%f,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,intput_column=%s,output_column=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
		chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,\
		float(g_m0[0]),float(g_m0[1]),float(g_m0[2]),crsampl,t0,f0,ntime,nfreq,dt,step_freq,window,fov,opt1,weightshm,ovt,ovf,intputcol,outputcol)
			if rate_noise['use'].upper()=="TRUE":
				command1 = "pyxis compression_noise[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,dtime=%i,chanwidthl=%f,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,noise=%f,ramdomseed=%s,intput_column=%s,output_column=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
                chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,crsampl,t0,f0,ntime,nfreq,dt,step_freq,window,fov,opt1,weightshm,ovt,ovf,noise,ramdom,intputcol,outputcol)
                else:
                    command = "pyxis src_supresion_fradius[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,stm=%.1f,stepm=%.1f,edm=%.1f,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,integrationl=%f,chanwidthl=%f,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,intput_column=%s,output_column=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
		chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,float(g_m0[0]),float(g_m0[1]),float(g_m0[2]),crsampl,t0,f0,ntime,nfreq,integrationl,step_freq,window,fov,opt1,weightshm,ovt,ovf,intputcol,outputcol)
		    if rate_noise['use'].upper()=="TRUE": 
			    command1 = "pyxis compression_noise[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,integrationl=%f,chanwidthl=%f,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,noise=%f,ramdomseed=%s,intput_column=%s,output_column=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
                chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,crsampl,t0,f0,ntime,nfreq,integrationl,step_freq,window,fov,opt1,weightshm,ovt,ovf,noise,ramdom,intputcol,outputcol)

        else:
	    print "here------------------------------------"
            if step_freq <= 0.:
                if df <= 0:
                     print """ Using Telescope simpling rate; Please give frequency_inc_hz or  numberfreqbintocompress """
                     sys.exit()
                else:
                	command = "pyxis src_supresion_fradius[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,stm=%.1f,stepm=%.1f,edm=%.1f,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,dfreq=%i,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,intput_column=%s,output_column=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
		chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,\
		float(g_m0[0]),float(g_m0[1]),float(g_m0[2]),crsampl,t0,f0,ntime,nfreq,df,window,fov,opt1,weightshm,ovt,ovf,intputcol,outputcol)
			if rate_noise['use'].upper()=="TRUE":
                            command1 = "pyxis compression_noise[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,dfreq=%i,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,noise=%f,ramdomseed=%s,intput_column=%s,output_column=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
		chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,crsampl,t0,f0,ntime,nfreq,df,window,fov,opt1,weightshm,ovt,ovf,noise,ramdom,intputcol,outputcol)

            else:
                command = "pyxis src_supresion_fradius[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,stm=%.1f,stepm=%.1f,edm=%.1f,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,chanwidthl=%f,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,intput_column=%s,output_column=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
		chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,\
		float(g_m0[0]),float(g_m0[1]),float(g_m0[2]),crsampl,t0,f0,ntime,nfreq,step_freq,window,fov,opt1,weightshm,ovt,ovf,intputcol,outputcol)
		print "here"
		
		if rate_noise['use'].upper()=="TRUE":
			 command1 = "pyxis compression_noise[confh=%s,nameh=%s,outputdir=%s,Ntimeh=%i,integrationh=%.1f,dech=%s,rah=%s,startfreqh=%.1f,nchanh=%i,chanwidthh=%.1f,starttimeh=%s,confl=%s,namel=%s,startfreql=%.1f,starttimel=%s,crs=%s,time0=%i,freq0=%i,ntime=%i,nfreq=%i,chanwidthl=%f,window=%s,fov=%f,optimalw=%s,weightsheme=%s,ovlaptime=%i,ovlapfreq=%i,noise=%f,ramdomseed=%s,intput_column=%s,output_column=%s]"%(conf1,name1,destdir1,hours1,integ1,dec1,ra1,startfreq1,nchan1,\
                chanwidth1,starttime1,conf2,name2,startfreq2,starttime2,crsampl,t0,f0,ntime,nfreq,step_freq,window,fov,opt1,weightshm,ovt,ovf,noise,ramdom,intputcol,outputcol)

        os.system(command)
	if command1 !="":
		os.system(command1)
        
if __name__ == "__main__":
    	arg = sys.argv
	Config = ConfigParser.ConfigParser()
	for i in range(1,len(arg)):
		Config.read("%s"%arg[i])
		sections = Config.sections()	
		hiresms = ConfigSectionMap(Config, sections[0])
		loresms = ConfigSectionMap(Config, sections[1])
		ptsrcsky = ConfigSectionMap(Config, sections[2])
		rate_noise = ConfigSectionMap(Config, sections[3])
		wsinc = ConfigSectionMap(Config, sections[4])
		wbessel = ConfigSectionMap(Config, sections[5])
		wbutter = ConfigSectionMap(Config, sections[6])
		sinchamming = ConfigSectionMap(Config, sections[7])
		sincbackman = ConfigSectionMap(Config, sections[8])
		sinchanning = ConfigSectionMap(Config, sections[9])
		besselhamming = ConfigSectionMap(Config, sections[10])
		besselblackman = ConfigSectionMap(Config, sections[11])
		besselhanning = ConfigSectionMap(Config, sections[12])
		boxcar = ConfigSectionMap(Config, sections[13])
		
        	if  wsinc['use'] == "True":
            		get_parameter_telescope(hiresms,loresms,wsinc,'sinc2',arg[i])
        	if  wbessel['use'] == "True":
            		get_parameter_telescope(hiresms,loresms,wbessel,'airy2',arg[i]);
		if  sinchamming['use'] == "True":
                        get_parameter_telescope(hiresms,loresms,sinchamming,'sinchamming2',arg[i]);
		if  sincbackman['use'] == "True":
                        get_parameter_telescope(hiresms,loresms,sincbackman,'sincblackman2',arg[i]);
		if  sinchanning['use'] == "True":
                        get_parameter_telescope(hiresms,loresms,sinchanning,'sinchanning2',arg[i]);
		if  besselhamming['use'] == "True":
                        get_parameter_telescope(hiresms,loresms,besselhamming,'besselhamming2',arg[i]);
		if  besselblackman['use'] == "True":
                        get_parameter_telescope(hiresms,loresms,besselblackman,'besselblackman2',arg[i]);
		if  besselhanning['use'] == "True":
                        get_parameter_telescope(hiresms,loresms,besselhanning,'besselhanning2',arg[i]);
		if  boxcar['use'] == "True":
                        get_parameter_telescope(hiresms,loresms,boxcar,'boxcar',arg[i]);


		


            
