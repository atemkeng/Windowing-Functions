SOURCE DENSITY AS A FUNCTION OF RADIUS
=======================================

1) The interface is not well construct. So please provide all input in your configuration (" One exemple is: wfconfigfileVLAC-Bandwidth-4MhZ_Sinc-0x0_FoV-1deg_int-100s.cfg")

2) Specify by "True" if you want the result of a given window in the configuration (NB: all windows are implemented except the butterword, so ignore this).

3) If you whant the FoV to be the Primary beam FWHM, put FoV=0 in the configuration file.

4) If you want to integrate with the telescope critical time, put "usecriticalsamplingrate=True".

5) Ignore "frequency_inc_hz=0. and step_time=0." if you specify already the resample time/frequency number of bins.

6) If you want to run windowing in one/deux directly speccify here: " exemple two directions freqdirection=True, timedirection=True "

7) The parameter: list_start_step_stop.grid_m0 =0,10,365. run from 0 to 365. spacing by 10arcmin

8) If you want to save a PDF file containing a table of comparison of noise ration specify here: [Compresion_rate_and_noise_analysis]. Make sure you have Tex install on your system

RUN CODE AS:

python runwf-fluxradius.py "ls wfconfigfileVLAC-Bandwidth-4MhZ_Sinc-0x0_FoV-1deg_int-100s.cfg" : for one configuration file

python runwf-fluxradius.py "ls wfconfigfileVLAC-Bandwidth-4MhZ_Sinc-0x0_FoV-1deg_int-100s.cfg" "ls wfconfigfileVLAC-Bandwidth-4MhZ_Sinc-0x0_FoV-1deg_int-50s.cfg" : for two configuration files and so on



9) If you want to run a bounch of telescopes / windows with differents integrations and/or frequency .... make many configuration files and run like

python runwf-fluxradius.py "ls wf*.cfg" : for all configurations files in the current directory, stating with wf

python runwf-fluxradius.py "ls Rep1/wf*.cfg" "ls Rep2/wf*.cfg"  : for all configuration files in Rep1 and Rep2 stating with wf

and so on .....


10) Specify your output directory, the process will defined a robus name of your output pickle data and save into the directory. The name depend strictly on your parameters in the configuration. 

Also specify the file name of the combine output data: for all your configurations files. 

cc: Oleg you know the structure of this pickle data. It is simpling "picklefile" by default in the attached configurations files.

11) A pipeline to run a full synthesis will be push very soon and the one of maximun integration and/or channel width.
