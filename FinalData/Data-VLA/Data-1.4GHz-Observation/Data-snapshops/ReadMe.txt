This readMe described the structure of our configurations files.

We have simulated a hi-resolution dataset at 1.4GHz sampled at each 1s over a duration of 6min40s. see below the configuration:


start Frequency or centre frequency = 1.4 MHz

step frequency or channels width of this hires dataset = 83400Hz

total number of frequency =  360

integration time = 1s

number of timeslots = 400

From this configuration, 300 timeslots (pershap 300s) will be allow for overlap Wfs acrross time and 240 channels (pershap 20MHz) will be alow for overlap Wfs accross channels. This is equivalent to time0 = 150 and freq0 = 120 in the configuration file.

Now, we then applied some suitable cases of BDWF using FoV of 1deg, 2deg, 4deg. See below:


1-) 

	- (100s, 10MHz) resampled across 100s and 10MHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 100 and numberfreqbintocompress = 120
	
	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency
	
	- Window-4x3:
	
		 we overlap over 300s across time and 20MHz across frequency. In the configuration file, this is equivalent to 
	  	 overlaptime = 150 and overlapfreq = 120 (150timeslots in the left hand side of the window and 150timeslots in the right hand 			 side. 120 channels in the left hand side and 120 channels in the right hand side)

		 calculations:
		
			total window width in timeslots=2*overlaptime+numbertimebintocompress=2*150+100=400timeslots
			total window width in second =400*1s=4* resampling time
	
			total window width in number of channels=2*overlapfreq+numberfreqbintocompress=2*120+120=360channels
			total window width in MHz=360*83400Hz=30MHz=3*resampling frequency
	

2-) 

	- (100s, 5MHz) resampled across 100s and 5MHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 100 and numberfreqbintocompress = 60

	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency

	- Window-4x3:

		we overlap over 300s across time and 10MHz across frequency. In the configuration file, this is equivalent to 
	  	overlaptime = 150 and overlapfreq = 60  (150timeslots in the left hand side of the window and 150timeslots 			in the right hand side. 60 channels in the left hand side and 60 channels in the right hand side)
		
		calculations:
		
			total window width in timeslots=2*overlaptime+numbertimebintocompress=2*150+100=400timeslots
			total window width in second =400*1s=4* resampling time
	
			total window width in number of channels=2*overlapfreq+numberfreqbintocompress=2*60+60=180channels
			total window width in MHz=180*83400Hz=15MHz=3*resampling frequency


3-) 

	- (100s, 2.5MHz) resampled across 100s and 2.5MHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 100 and numberfreqbintocompress = 30

	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency
	
	 Window-4x3:

		we overlap over 300s across time and 5MHz across frequency. In the configuration file, this is equivalent to 
	  	overlaptime = 150 and overlapfreq = 30  (150timeslots in the left hand side of the window and 150timeslots 			in the right hand side. 30 channels in the left hand side and 30 channels in the right hand side)
		
		calculations:
		
			total window width in timeslots=2*overlaptime+numbertimebintocompress=2*150+100=400timeslots
			total window width in second =400*1s=4* resampling time
	
			total window width in number of channels=2*overlapfreq+numberfreqbintocompress=2*30+30=90channels
			total window width in MHz=90*83400Hz=7.5MHz=3*resampling frequency


4-) 

	- (50s, 10MHz) resampled across 50s and 10MHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 50 and numberfreqbintocompress = 120
	
	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency
	
	- Window-4x3:
	
		 we overlap over 150s across time and 20MHz across frequency. In the configuration file, this is equivalent to 
	  	 overlaptime = 75 and overlapfreq = 120 (75timeslots in the left hand side of the window and 75timeslots in the right hand 			 side. 120 channels in the left hand side and 120 channels in the right hand side)

		 calculations:
		
			total window width in timeslots=2*overlaptime+numbertimebintocompress=2*75+50=200timeslots
			total window width in second =200*1s=4* resampling time
	
			total window width in number of channels=2*overlapfreq+numberfreqbintocompress=2*120+120=360channels
			total window width in MHz=360*83400Hz=30MHz=3*resampling frequency
	

5-) 

	- (50s, 5MHz) resampled across 50s and 5MHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 50 and numberfreqbintocompress = 60

	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency

	- Window-4x3:

		we overlap over 150s across time and 10MHz across frequency. In the configuration file, this is equivalent to 
	  	overlaptime = 75 and overlapfreq = 60  (75timeslots in the left hand side of the window and 75timeslots 			in the right hand side. 60 channels in the left hand side and 60 channels in the right hand side)
		
		calculations:
		
			total window width in timeslots=2*overlaptime+numbertimebintocompress=2*75+50=200timeslots
			total window width in second =200*1s=4* resampling time
	
			total window width in number of channels=2*overlapfreq+numberfreqbintocompress=2*60+60=180channels
			total window width in MHz=180*83400Hz=15MHz=3*resampling frequency


6-) 

	- (50s, 2.5MHz) resampled across 50s and 2.5MHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 50 and numberfreqbintocompress = 30

	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency
	
	 Window-4x3:

		we overlap over 150s across time and 5MHz across frequency. In the configuration file, this is equivalent to 
	  	overlaptime = 75 and overlapfreq = 30  (75timeslots in the left hand side of the window and 75timeslots 			in the right hand side. 30 channels in the left hand side and 30 channels in the right hand side)
		
		calculations:
		
			total window width in timeslots=2*overlaptime+numbertimebintocompress=2*75+100=200timeslots
			total window width in second =200*1s=4* resampling time
	
			total window width in number of channels=2*overlapfreq+numberfreqbintocompress=2*30+30=90channels
			total window width in MHz=90*83400Hz=7.5MHz=3*resampling frequency


7-) 

	- (25s, 10MHz) resampled across 25s and 10MHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 25 and numberfreqbintocompress = 120
	
	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency
	
	- Window-4x3:
	
		 we overlap over 75s across time and 20MHz across frequency. In the configuration file, this is equivalent to 
	  	 overlaptime = 38 and overlapfreq = 120 (38timeslots in the left hand side of the window and 38timeslots in the right hand 			 side. 120 channels in the left hand side and 120 channels in the right hand side)

		 calculations:
		
			total window width in timeslots=2*overlaptime+numbertimebintocompress=2*38+25=101timeslots
			total window width in second =101*1s~4* resampling time
	
			total window width in number of channels=2*overlapfreq+numberfreqbintocompress=2*120+120=360channels
			total window width in MHz=360*83400Hz=30MHz=3*resampling frequency
	

8-) 

	- (25s, 5MHz) resampled across 25s and 5MHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 25 and numberfreqbintocompress = 60

	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency

	- Window-4x3:

		we overlap over 75s across time and 10MHz across frequency. In the configuration file, this is equivalent to 
	  	overlaptime = 38 and overlapfreq = 60  (38timeslots in the left hand side of the window and 38timeslots 			in the right hand side. 60 channels in the left hand side and 60 channels in the right hand side)
		
		calculations:
		
			total window width in timeslots=2*overlaptime+numbertimebintocompress=2*38+25=101timeslots
			total window width in second =101*1s~4* resampling time
	
			total window width in number of channels=2*overlapfreq+numberfreqbintocompress=2*60+60=180channels
			total window width in MHz=180*83400Hz=15MHz=3*resampling frequency


9-) 

	- (25s, 2.5MHz) resampled across 25s and 2.5MHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 25 and numberfreqbintocompress = 30

	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency
	
	 Window-4x3:

		we overlap over 75s across time and 5MHz across frequency. In the configuration file, this is equivalent to 
	  	overlaptime = 38 and overlapfreq = 30  (38timeslots in the left hand side of the window and 38timeslots 			in the right hand side. 30 channels in the left hand side and 30 channels in the right hand side)
		
		calculations:
		
			total window width in timeslots=2*overlaptime+numbertimebintocompress=2*38+25=101timeslots
			total window width in second =101*1s~4* resampling time
	
			total window width in number of channels=2*overlapfreq+numberfreqbintocompress=2*30+30=90channels
			total window width in MHz=90*83400Hz=7.5MHz=3*resampling frequency



