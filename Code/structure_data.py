import ConfigParser
import sys;
import os;
import pickle
import numpy

def list_data(listdataextension=None, destir=None, datafile=None, shellcom=None, savenoise=None):
        List    = []
        windows = []
        Fovs = []
        Timeint = []
        Bandwidth = []
	for i in range(1,len(listdataextension)):
		data = pickle.load( open( listdataextension[i], "rb"))
                bw = data['Bandwidth']
                flux = data['flux']
                fov = data['fov']
                wf = data['window']
                radius = data['radius']
                noise = data['noise']
                timeint = data['int']
                try :
                    windows.index(wf)
                except :
                    windows.append(wf)
                    
                try :
                    Fovs.index('%i'%int(fov))
                except :
                    Fovs.append('%i'%int(fov))
                try :
                    Timeint.index('int%is'%int(timeint))
                except :
                    Timeint.append('int%is'%int(timeint))
                try :
                    Bandwidth.index(int(bw))
                except :
                    Bandwidth.append(int(bw))
                
                List.append((int(bw),timeint,fov, wf ,flux,noise,radius ))

        bandwidth = dict(zip(range(0,len(Bandwidth)),Bandwidth))

        FoVdeg= dict(zip(range(0,len(Fovs)), Fovs))
        data100 = dict(zip(range(0,len(Timeint)),Timeint))
        # print "bandwidth",bandwidth
        # print "integaration",data100
        # print "FoV",FoVdeg
        # stop
        #for k in range(len(List)):
        trees ={}
        for keybw in bandwidth:
                bandwindths = {'bandwidth':'%iMhz'%bandwidth[keybw]}
                for keytime in data100:

                 varioustime = {} 
                 Avg = {}
                 for keyfov in FoVdeg:
                    dataFoV ={}
                        
                    for elt  in range(len(List)):
                        #print"len",len(List),elt
                        bw ,timeint, fov, window, flux, noise,radius = List[elt]
                        # print'', bandwidth[keybw],bw
                        if bandwidth[keybw] == bw and data100[keytime] == "int%is"%int(timeint) and FoVdeg[keyfov]=='%i'%int(fov):
                            #print'', bandwidth[keybw],bw,data100[keytime],"int%is"%int(timeint),FoVdeg[keyfov],'%i'%int(fov)
                            leaft = {"flux":flux, "noise":noise}
                            if window.split('-')[0] == "boxcar":
                                Avg = leaft
                            else:
                                    a=window.split('-')
                                    a1=a[0].split('2')
                                    
                                    
                                    if a1[0]=='sincsinc':
                                       a1[0] = 'sinc2'
                                    if a1[0]=='airy':
                                       a1[0]='bessel'
                                    #print a1[0]
                                    window=a1[0]+'-'+a[1]
                                    dataFoV.update({window: leaft})
                            
                            

                    varioustime.update({"FoV%sdeg"%FoVdeg[keyfov]:dataFoV}) 
                
                 varioustime.update({"avg":Avg})
                 varioustime.update({"radius":radius})
                 #print "keytime",data100[keytime]
                 bandwindths.update({data100[keytime]:varioustime})
                 trees.update({"bandwidth%iMhz"%bandwidth[keybw]:bandwindths})
                os.system(shellcom)
                f = open("%s/%s.data"%(destir, datafile), 'wb');    
                pickle.dump(trees, f)
		
                f.close()
	
