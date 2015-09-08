
This readMe described the structure of our configurations files.

We have simulated a hi-resolution dataset at 16GHz sampled at each 0.01s for a snapshot of 1s. see below the configuration:


start Frequency or centre frequency = 16 MHz

step frequency or channels width of this hires dataset = 500.Hz

total number of frequency =  25

integration time = 0.01s

number of timeslots = 100

Overlaps windows are not applied, This is equivalent to time0 = 0 and freq0 = 0 in the configuration file.

Now, we then applied a BDWF with FoV of 0.7deg=42arcmin See below:


1-) 

	- (0.1s, 12.5KHz) resampled across 0.1s and 12.5KHz. In the configuration file this is equivalent to 
	    numbertimebintocompress = 10 and numberfreqbintocompress = 25
	
	- Window-1x1: 0s overlap in time and 0MHz overlap in frequency
