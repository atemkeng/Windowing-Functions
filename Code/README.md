Source density as function of distance from the phase centre


The interface is not well construct. So please provide all input in your configuration (" One exemple is: wfconfigfileVLAC-Bandwidth-4MhZ_Sinc-0x0_FoV-1deg_int-100s.cfg")

Specify by "True" if you want the result of a given window in the configuration (NB: all windows are implemented except the butterword, so ignore this).

If you whant the FoV to be the Primary beam FWHM, put FoV=0 in the configuration file. 


If you want to integrate with the telescope critical time, put  "usecriticalsamplingrate=True".


Ignore "frequency_inc_hz=0. and step_time=0." if you specify already the resample time/frequency number of bins.


If you want to run windowing in one/deux directly speccify here: " exemple two directions freqdirection=True, timedirection=True "

The parameter: list_start_step_stop.grid_m0 =0,10,365. run from 0 to 365. spacing by 10arcmin

If you want to save a PDF file containing a table of comparison of noise ration specify here: [Compresion_rate_and_noise_analysis]. Make sure you have Tex install on your system


run the code as:

python runwf-fluxradius.py wfconfigfileVLAC-Bandwidth-4MhZ_Sinc-0x0_FoV-1deg_int-100s.cfg



If you want to run a bounch of telescopes / windows with differents integrations and/or frequency .... make many configuration files and run like

python runwf-fluxradius.py file1.cfg file2.cfg file3.cfg .......


Specify your output directory, the process will defined a robus name of your output pickle data  and save into the directory. The name depend strictly on your parameters in the configuration. An update of this work will allow you to give the naming your self.


A pipeline to run a full synthesis will be push very soon and the one of maximun integration and/or channel width.





