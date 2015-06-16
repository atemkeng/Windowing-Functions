import numpy as np
from pyrap.tables import table
import os
import sys
import pylab

class MSUVDISTANCE (object):
    """Class for reading and resampling data from an MS"""

    def __init__ (self,msname,column="DATA",time0=0,ntime=None,freq0=0,nfreq=None):
        """Inits with given MS and column name.
        If time0/ntime and freq0/nfreq is given, handles only a subset of the data,
        otherwise reads all channels and timeslots.
        """;
        self.msname = msname;
        tab = table(msname,ack=False,readonly=False)
        self.A0 = A0 = tab.getcol('ANTENNA1')
        self.A1 = A1 = tab.getcol('ANTENNA2')
        data = tab.getcol(column);
        if nfreq is None:
            nfreq = data.shape[1]-freq0;
            #self.data = data = data[:,freq0:freq0+nfreq,:];
        self.data = data.copy();
        self.nfreq = nfreq;
        self.freq0=freq0
        print "Visibility column shape:",data.shape
        self.na = na = np.max(A1)+1
        self.nbl = (na*(na-1))/2+na
        self.ncorr = data.shape[2]
        self.nbins = data.shape[0]
        self.UVW = tab.getcol("UVW")
            # get number of timeslots. This assumes a regular MS staructure (i.e. first baseline is
            # present)
        ntimes = (data[(A0==0)&(A1==1)]).shape[0]
            # actual number of timeslot to use given by time0 and ntime
        self.ntime = ntime or ntimes - time0;
        self.time0 = time0;
        self.pmax = 17;
        self.qmax = 26;
        # get frequency and wavelength (per channel) from SPECTRAL_WINDOW subtable
        t2 = table(tab.getkeyword("SPECTRAL_WINDOW"),readonly=False)
        self.channels = t2.getcol("CHAN_FREQ",0);
        self.freqs = t2.getcol("CHAN_FREQ",0)[0,freq0:freq0+nfreq];
      # print "Frequencies are",self.freqs;
        self.wavel = 3e8/self.freqs;
      # print "Wavelengths are",self.wavel;
        t2.close()
        tab.close()

    def boxcar (self,dfreq=None, dtimelpq=None):
        """Downsamples data using normal boxcar averaging over windows of size dtime x dfreq.
        Returns list of (antenna1,antenna2,averaged_data) tuples, one per each baseline, where
        averaged_data has the shape (ntime/dtime,nfreq/dfreq,ncorr)."""
        uv = self.UVW[(self.A0==self.pmax)&(self.A1==self.qmax)];
        uvd_bmax = np.sqrt(((uv[self.time0:dtimelpq+self.time0,:2])**2).sum(axis=1))
        # this is the output shape
        nfreq1 = self.nfreq/dfreq;
        # make a list of per-baseline output arrays 
        result = [];
        # # loop over each baseline
        for p in range(self.na):
            for q in range(p+1,self.na):
                input_index = (self.A0==p)&(self.A1==q);
                uv = self.UVW[input_index].copy();
                uv = uv[self.time0:self.time0+self.ntime,:2];
                uvd = np.sqrt((uv**2).sum(axis=1));
                ## get the output shape of the time direction for this baseline
                timeslots = int(uvd.sum(axis=0) // uvd_bmax.sum(axis=0));
                ## get the corresponding data to work with
                datapq = self.data[input_index].copy();
                ## get the number of time bins to average for this baseline
                dtimepq = self.ntime//timeslots
                timeslots = self.ntime//dtimepq
                ## get the total number of time bins to work with
                ntimepq = dtimepq*timeslots
                datapq = datapq[self.time0:self.time0+ntimepq,self.freq0:self.freq0+self.nfreq,:];
                # reshape so that e.g. ntime becomes ntime1,dtime, so that we can then reduce over the second axis
                print "nfreq1 = %i, dfreq = %i, ntimepq = %f, dtimepq = %f, datapq.shape=%s"%(nfreq1,dfreq,ntimepq,dtimepq,str(datapq.shape))
                datapq = datapq.reshape((timeslots,dtimepq,nfreq1,dfreq,self.ncorr))
                # take mean over the new axes, this becomes our result
                result.append( (p,q,datapq.mean(3).mean(1)) );
        return result;
 
    def overlap_window (self,window_function,dtimelpq,dfreq,overlap_time=0,overlap_freq=0,dump=None,dumpfig=None):
      """Downsamples data using the specified window_function. Window size is dtime x dfreq nominally,
      but overlap_time makes the effective window extend by that number of timeslots (hi-res) in each
      direction. Overlap_freq does the same for frequency (currently not implemented).

      window_function(x,y) is a function taking two array of uv-distances (in wavelengths), and returning the
      corresponding weights.

      Returns list of (antenna1,antenna2,averaged_data) tuples, one per each baseline, where
      averaged_data has the shape (ntime/dtime,nfreq/dfreq,ncorr)."""
      uv = self.UVW[(self.A0==self.pmax)&(self.A1==self.qmax)];
      uvd_bmax = np.sqrt(((uv[self.time0:dtimelpq+self.time0,:2])**2).sum(axis=1))
      # effective window size in frequency
      dfreq1 = dfreq + overlap_freq*2;
      # this is the output shape for the frequency direction
      nfreq1 = self.nfreq/dfreq;
      # make a list of per-baseline output arrays 
      result = [];
      # loop over each baseline
      freqs = self.channels[0,:];
      wavel = 3e8/freqs;
      for p in range(self.na):
        for q in range(p+1,self.na):
          input_index = (self.A0==p)&(self.A1==q)
          data1 = self.data[input_index].copy();
          uv1 = self.UVW[input_index,:2].copy();
          uv2 = uv1[self.time0:self.time0+self.ntime,:2];
          uvd = np.sqrt((uv2**2).sum(axis=1));
          ## get the output shape of the time direction for this baseline
          timeslots = int(uvd.sum(axis=0) // uvd_bmax.sum(axis=0));
          ## get the number of time bins to average for this baseline
          dtimepq = self.ntime//timeslots
          ## get the total number of time bins to work with
          ntimepq = dtimepq*timeslots
          timeslots = self.ntime//dtimepq
    #      data1 = data1[self.time0:self.time0+ntimepq,self.freq0:self.freq0+self.nfreq,:];
    #      uv1 = uv1[self.time0:self.time0+ntimepq,:2];
          #print "nfreq1 = %i, dfreq = %i, ntimepq = %f, dtimepq = %f, datapq.shape=%s"%(nfreq1,dfreq,ntimepq,dtimepq,str(data1.shape))
          resdata = np.zeros((timeslots,nfreq1,self.ncorr),self.data.dtype);
          # loop over all output timeslots
          if p==17 and q==26:
                print data.shape, dtimepq,ntimepq,timeslots
          dtime1 = dtimepq + overlap_time*2;
          for out_time in range(timeslots):
            # make timeslice corresponding to extended window
            in_time0 = self.time0+out_time*dtimepq-overlap_time;
            timeslice = slice(in_time0,in_time0+dtime1);
            # extract data subset over the window
            data = data1[timeslice,...].copy();
	    #loop over all output frequencies
            #print "*********dfreq=%f; nfreq1=%f, data=%s*****"%(dfreq,nfreq1,str(data.shape))
       	    for out_freq in range(nfreq1):
		in_freq0 = self.freq0 + out_freq*dfreq - overlap_freq;
		freqslice = slice(in_freq0, in_freq0 + dfreq1);
		#print "in_freq=%f, freq0=%f,out_freq=%f,dfreq=%f,overlap_freq=%f; dfreq1=%f"%(in_freq0,self.freq0,out_freq,dfreq,overlap_freq,dfreq1)
		data2 = data[:,freqslice,...].copy();
                import Pyxis;
                Pyxis.info("data2=%s;in_freq=%i, freq0=%i,out_freq=%i,dfreq=%i,overlap_freq=%i; dfreq1=%i"\
                    %(str(data2.shape),in_freq0,self.freq0,out_freq,dfreq,overlap_freq,dfreq1))
                #print "data.shape=%s; data2.shape=%s"%(data.shape,data2.shape)
		#print "wavel shape******=",wavel.shape
		wavel1 = wavel[freqslice];
                # get UV coordinates for the window, shape is DT1,2
            	uv = uv1[timeslice];
                # get corresponding uv-center point (uv0), shape is 2  
            	if dtime1%2:  # odd dtime: we have a uv-bin at the centre, so start at n/2, then take every n-th uv point
              		uv0 = uv[dtime1/2];
            	else: # even dtime -- take the middle of the two centre bins
              		uv0 = (uv[dtime1/2-1,:] + uv[dtime1/2,:])/2;
                # take uv-distance component along time, shape is DT1
            	uvd1 = np.sqrt(((uv-uv0[np.newaxis,:])**2).sum(1));
                # convert it to wavelengths at centre channel of each window, this will now be of shape DT1,NF1
		if dfreq1%2:
			wl0 = wavel1[dfreq1/2];
		else:
			wl0 = (wavel1[dfreq1/2 -1] + wavel1[dfreq/2])/2;
            	uvd1 = uvd1[:,np.newaxis] / wl0;
                # now for the second component: ||uv0||/wl0 at each window is the centre channel uv-length,
                # the second component's distance is ||uv0||/wl - ||uv0||/wl0
                # first compute ||uv0||, this is scalar
            	uv0r = np.sqrt((uv0**2).sum());  
                # work out second component of uv-distance, shape will be NF1,DF
            	uvd2 = uv0r/wavel1 - uv0r/wl0;
                # evaluate windowing function
            	wf = window_function(uvd1[:,:],uvd2[np.newaxis,:]);
                # sum of wf over whole window and each frequency bin, shape NF1
            	wfsum = wf.sum(1).sum(0);
                #print"",out_time,out_freq,data2,wf[...,np.newaxis]#,((data2*wf[...,np.newaxis]).sum(1).sum(0))/wfsum
		resdata[out_time,out_freq,...] = ((data2*wf[...,np.newaxis]).sum(1).sum(0))/wfsum;
          result.append((p,q,resdata));
      return result;


def save_visibility_arrays (msname,arrays,column="DATA"):
  """Saves a set of visibility arrays to the MS.
  arrays is a list of (p,q,data) tuples, where p & q are antenna indices,
  and data is the corresponding array of shape ntime,nfreq,ncorr.
  The shape must match the shape of the MS, or errors will result. 
  """
  os.system("addbitflagcol %s"%msname)
  tab = table(msname,readonly=False,ack=False);
  # read in data column and antenna indices, fill data with 0
  data = tab.getcol(column)
  a1 = tab.getcol("ANTENNA1")
  a2 = tab.getcol("ANTENNA2")
  data.fill(0.);
  FLAG = tab.getcol("FLAG_ROW")
  #FLAG.fill(0);
  # fill data from arrays
  slotmax = np.max(np.array([len(dataxx) for p,q,dataxx in arrays]))
  dataout = np.zeros_like(data)
  intdex = (a1==0)&(a2==1)
  dat = data[intdex]
  datapq = np.zeros_like(dat)
  tab.close()
  arr1 = np.zeros_like(data[(a1==0)&(a2==1)])
  for p,q,arr in arrays:
      arr1[0:arr.shape[0],:,:] = arr.copy()
      data[(a1==p)&(a2==q),...] = arr1.copy();
      os.system("flag-ms.py -S %i,%i -T %i:%i -f %i:%i -c %s"%(p,q,arr.shape[0],arr1.shape[0],arr.shape[0],arr1.shape[0],msname));
      arr1.fill(0.)
  # for p,q,arr in arrays:
  #     import Pyxis;
  #     Pyxis.info("**********(%i,%i)*********"%(p,q))
  #     datapq.fill(0.)
  #     stepslots  = slotmax//len(arr)-1
  #     if len(arr)==len(datapq):
  #            data[(a1==p)&(a2==q),...]  = arr;
  #     else:
  #         timeslot = 0;
          
  #         for out_time in range(len(arr)):
  #             datapq[timeslot,:,:] = arr[out_time,:,:]
  #             if out_time%2 == 0:
  #                 os.system("flag-ms.py -S %i,%i -T %i:%i -f %i:%i -c %s"%(p,q,timeslot+1,timeslot+1+stepslots,timeslot+1,timeslot+1+stepslots,msname));
  #                 #FLAG[timeslot+1:timeslot+1+stepslots] = True;
  #                 timeslot+=stepslots
  #             else :
  #                 os.system("flag-ms.py -S %i,%i -T %i:%i -f %i:%i -c %s"%(p,q,timeslot+1,timeslot+2+stepslots,timeslot+1,timeslot+2+stepslots,msname));
  #                 #FLAG[timeslot+1:timeslot+2+stepslots] = True;
  #                 timeslot+=stepslots+1
  #         data[(a1==p)&(a2==q),...] = datapq.copy()
  # # write out
  tab = table(msname,readonly=False,ack=False);
  tab.putcol(column,data)
  tab.close();

# def flaggindata(self,msname,arrays,time1,column="CORRECTED_DATA",irgg=0):
#         """Saves a set of visibility arrays to the MS.
#         arrays is a list of (p,q,data) tuples, where p & q are antenna indices,
#         and data is the corresponding array of shape ntime,nfreq,ncorr.
#         The shape must match the shape of the MS, or errors will result. 
#         """
#         tab = table(msname,readonly=False,ack=False);
#   # read in data column and antenna indices, fill data with 0
#         data = tab.getcol(column)
#         a1 = tab.getcol("ANTENNA1")
#         a2 = tab.getcol("ANTENNA2")
#         data.fill(0.);
#         tab.close()
#   # fill data from arrays
#         arr1 = np.zeros_like(data[(a1==0)&(a2==1)])
#         for p,q,arr in arrays:
#             arr1[0:arr.shape[0],:,:] = arr[:,np.newaxis,:].copy();
#             data[(a1==p)&(a2==q),...] = arr1.copy();
#             #[0:arr.shape[0],:,:] = arr[:,np.newaxis,:].copy();
#             #print "shape of arr",p,q,#arr[:,np.newaxis,:]==data[(a1==p)&(a2==q)][0:arr.shape[0],:,:] #arr.shape[0]:time1
#             if irgg==0:
#                 os.system("flag-ms.py -S %i,%i -T %i:%i -f %i:%i -c MeerKAT-snapshot-21cm.MS"%(p,q,arr.shape[0],time1,arr.shape[0],time1));
#             #os.system(flag)
#   # write out
#         tab = table(msname,readonly=False,ack=False);
#         print "data",data
#         print "ar",(arr.reshape(len(arr),1,4)).dtype
#         tab.putcol(column,data)
#         print "save succeded"
#         tab.close();
   
