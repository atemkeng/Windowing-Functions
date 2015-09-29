
We used 9 EVN GROUND-BASED stations :

EF, HH, JB, NT, ON, TR, YS, WB, SH. The longest baseline here is ~10161Km and the shortest baseline is
~267Km. 

We sampled the uv domain during 0.001s time step during a period  2.5s at 1.6GHz, this resulting to 2500 timeslots with a total bandwidth of 69KHz channelized into 138 channels each of width 500Hz. 


1-) 0.5s and 23KHz resampling time and frequency: 

	for sinc-1x1 or bessel-1x1 : dtime=500 and dfreq=46, overlap_time =0s  overlap_freq=0HZ
	
	for sinc-4x3 : dtime=500 and dfreq=46, overlap_time =750*2=3xdtime  overlap_freq=46*2 = 2dfreq



2) 0.25s and 23KHz resampling time and frequency: 

	for sinc-1x1 or bessel-1x1 : dtime=250 and dfreq=46, overlap_time =0s  overlap_freq=0HZ
	
	for sinc-4x3 : dtime=250 and dfreq=46, overlap_time =375*2=3xdtime  overlap_freq=46*2 = 2dfreq



3) 0.1s and 23KHz resampling time and frequency: 

	for sinc-1x1 or bessel-1x1 : dtime=100 and dfreq=46, overlap_time =0s  overlap_freq=0HZ
	
	for sinc-4x3 : dtime=100 and dfreq=46, overlap_time =150*2=3xdtime  overlap_freq=46*2 = 2dfreq
