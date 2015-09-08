This readMe described the structure of our configurations files.

We have simulated a hi-resolution dataset at 14GHz sampled at each 1s over a duration of 30min30s. see below the configuration:


start Frequency or centre frequency = 14 GHz

step frequency or channels width of this hires dataset = 1000000.Hz=1MHz

total number of frequency =  220

integration time = 1s

number of timeslots = 1830

From this configuration,  30 timeslots (pershap 30s) will be allow for overlap Wfs across time and 20 channels (pershap 20MHz) will be alow for overlap Wfs accross channels. This is equivalent to time0 = 15 and freq0 = 10 in the configuration file.

Now, we then applied some suitable cases of BDWF using FoV of 1deg, 2deg, 4deg. See below:


1-) 

	- (10s, 10MHz) resampled across 10s and 10MHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 10 and numberfreqbintocompress = 10
	
	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency
	
	- Window-4x3:
	
		 we overlap over 30s across time and 20MHz across frequency. In the configuration file, this is equivalent to 
	  	 overlaptime = 15 and overlapfreq = 10 (15timeslots in the left hand side of the window and 15timeslots in the right hand 			 side. 10 channels in the left hand side and 10 channels in the right hand side)

		 calculations:
		
			total window width in timeslots=2*overlaptime+numbertimebintocompress=2*15+10=40timeslots
			total window width in second =40*1s=4* resampling time
	
			total window width in number of channels=2*overlapfreq+numberfreqbintocompress=2*10+10=30channels
			total window width in MHz=30*1000000Hz=30MHz=3*resampling frequency

****Notice that I used the same configuration for the  uv amplitude Vs the source radius, here I taugh the amplitude of the lores centre uv bin and a 4deg sinc-4x3

