

Hires:

integration time = 0.001s, number of timeslots = 100

	I fixed 10 timeslots to average and fixed  overlap_left=45, overlap_rigth=45

channels widths = 500KHz, number of channels = 2500

	I fixed 500 channels to average  and overlap_left=1000, overlap_rigth=1000



LOres:

data:  vlbi-FoV0.26deg-int0.01s-bandwidth0.12MHz.data

 	resample time 0.01s (10 timeslots averaged), resample channels widths = 0.125MHZ (250 channels averaged)

	bessel-4x3:   overlap_time = 2*15, overlap_freq = 2*250
	
	bessel-4x5:   overlap_time = 2*15, overlap_freq = 2*500

	bessel-5x1:   overlap_time = 2*20, overlap_freq = 0
	
	bessel-10x1:  overlap_time = 2*45, overlap_freq = 0	


 
data:  vlbi-FoV0.26deg-int0.01s-bandwidth0.25MHz.data

	
	resample time 0.01s (10 timeslots averaged), resample channels widths = 0.25MHZ (500 channels averaged)
	
	bessel-4x3:   overlap_time = 2*15, overlap_freq = 2*500
	
	bessel-4x5:   overlap_time = 2*15, overlap_freq = 2*1000
	
	bessel-5x1:   overlap_time = 2*20, overlap_freq = 0
	
	bessel-10x1:  overlap_time = 2*45, overlap_freq = 0	


