This readMe describedt he structure of our configurations files.

We have simulated a hi-resolution dataset at 1.4GHz sampled at each 1s over a duration of 35min. see below the configuration:


start Frequency or centre frequency = 1.4 GHz

step frequency or channels width of this hires dataset = 1000000.Hz=1MHz

total number of frequency =  220

integration time = 1s

number of timeslots = 2100

From this configuration,  300 timeslots (pershap 300s) will be allow for overlap Wfs across time and 20 channels (pershap 20MHz) will be alow for overlap Wfs accross channels. This is equivalent to time0 = 150 and freq0 = 10 in the configuration file.

Now, we then applied some suitable cases of BDWF using FoV of 1deg, 2deg, 4deg. See below:


1-) 

	- (100s, 10MHz) resampled across 100s and 10MHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 100 and numberfreqbintocompress = 10
	
	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency
	
	- Window-4x3:
	
		 we overlap over 300s across time and 20MHz across frequency. In the configuration file, this is equivalent to 
	  	 overlaptime = 150 and overlapfreq = 10 (150timeslots in the left hand side of the window and 150timeslots in the right hand 			 side. 10 channels in the left hand side and 10 channels in the right hand side)

		 calculations:
		
			total window width in timeslots=2*overlaptime+numbertimebintocompress=2*150+100=400timeslots
			total window width in second =400*1s=4* resampling time
	
			total window width in number of channels=2*overlapfreq+numberfreqbintocompress=2*10+10=30channels
			total window width in MHz=30*1000000Hz=30MHz=3*resampling frequency


