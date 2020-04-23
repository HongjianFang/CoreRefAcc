#/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
# File Name : main.py
#
# Purpose : ACC to extract PcP
#
# Creation Date : 05-07-2019
#
# Last Modified : Wed 22 Apr 2020 04:45:58 PM EDT
#
# Created By : Hongjian Fang: hfang@mit.edu 
#
#_._._._._._._._._._._._._._._._._._._._._.*/

#<<<<<<< HEAD:main_new.py
def autocorr_td(tr,stalign_a,stalign_b,conlen=15):
        N = tr.stats.npts
        delta = tr.stats.delta
        yf = fft(tr.data)
        ayf = np.abs(yf)
        myf = np.convolve(ayf, np.ones((conlen,))/conlen, mode='same')
        yf = yf/(myf+1.0e-10)
        xx = np.real(ifft(yf))
        
        stalign_a = int(stalign_a/delta)
        stalign_b = int(stalign_b/delta)
        winlen = stalign_b - stalign_a
        xxw = xx[stalign_a:stalign_b]
        acorr = np.correlate(xx,xxw,mode='full')
        #acorr = acorr[winlen:]
#       # maxloc = np.argmax(abs(acorr))
        #maxloc = stalign_a
        #acorr = np.roll(acorr,-maxloc)
        
#        tr.data[:len(acorr)] = acorr
#        return tr
        return acorr

def autocorr_fd(tr,conlen=10):
        N = tr.stats.npts
        delta = tr.stats.delta
        data = np.zeros(2*N,)
        data[:N] = tr.data
        yf = fft(data)
        ayf = np.abs(yf)
        myf = np.convolve(ayf, np.ones((conlen,))/conlen, mode='same')
        yf = yf/(myf+1.0e-10)
        myf = yf*np.conj(yf)
        xx = np.real(ifft(myf))
        tr.data = xx[:N]
        return tr

def process_event(ievent1):
    if ievent1>=nevent:
        return 0 
    ds = h5py.File(tempdata+'/'+str(ievent1)+'.h5','w')
#    ievent1 = ievent
    evname = evelist[ievent1].split('/')[-1]

    evidx1 = evid.index(evname)
    evlat1 = eqs[evidx1]['latitude']
    evlon1 = eqs[evidx1]['longitude']
    evdep1 = eqs[evidx1]['depth']
    if evdep1 > 80:
        return 0
    evtime1 = eqs[evidx1]['datetime']
    evmag1 = eqs[evidx1]['magnitude']  
    stalist1 = glob.glob('/'.join([eqdir,evname,datatp,'*.'+oricomp]))
    stalist = [stal.split('/')[-1] for stal in stalist1]
    irissta = pd.read_csv('/'.join([eqdir,evname,'info','station_event']),header=None,names=\
              ('net','sta','loc','channel','lat','lon','ele','tmp1','datasrc','startt',\
              'evlat','evlon','evdep','mg','tmp2','tmp3','tmp4','temp5'),keep_default_na=False)

    nsta = len(stalist1)
    stapos = np.zeros((nsta,2))
    print ('the',ievent,'th event of all:',nevent)
    print ('%s %s %4.1f  %7.2f %7.2f %7.2f %d\n' % ('evinfo :',evname,evmag1,evlat1,evlon1,evdep1,nsta))

    for ista in range(nsta):
                stanet = stalist[ista].split('.')[0]
                staname = stalist[ista].split('.')[1]
                location = stalist[ista].split('.')[2]
                stasub = irissta[(irissta['net']==stanet) & (irissta['sta']==staname)]
                if len(stasub)<1:
                    continue
                stlat = stasub.iloc[0]['lat']
                stlon = stasub.iloc[0]['lon']
                
                dis1 = distaz.DistAz(evlat1,evlon1,stlat,stlon)

                if dis1.delta<mindist or dis1.delta>maxdist:
                    continue

                parr = ftelep(dis1.delta,evdep1)[0]

                trace = '/'.join([eqdir,evname,datatp,stalist[ista]])
                if os.path.getsize(trace) < 10000:
                    continue
                strm = obspy.read(trace)

                #if comp == 'BHR' or comp == 'BHT':
                if comp in ['BHR','BHT']:
                    tracen = trace.split('.')[:-1]
                    tracen = '.'.join(tracen+['BHN'])
                    if not os.path.isfile(tracen) or os.path.getsize(tracen) < 10000:
                      continue
#                    print('begin...',ista)
                    strm += obspy.read(tracen)
                    strm.merge()
                    #strm.trim(evtime1+parr-trimb,evtime1+parr+trima,pad=True,fill_value=0)
                    strm.trim(evtime1+parr-trimb-50,evtime1+parr+trima+50,pad=True,fill_value=0)
                    strm.resample(rsample)
                    #print(strm)
                    strm.rotate('NE->RT',back_azimuth=dis1.baz)
                    tracet = strm.select(channel=comp)
#                    print('finish...')
                    strm = tracet
#                    print('finish...sampling')
#                    strm.remove(tracet[0])
                else:
                    strm.merge()
                    strm.trim(evtime1+parr-trimb-50,evtime1+parr+trima+50,pad=True,fill_value=0)
                    strm.resample(rsample)

                tr = strm[0]
                #if tr.stats.npts<100:
                #    continue
                #tr.resample(rsample)
                
                tr.stats.distance = dis1.delta
                #tr.trim(evtime1+parr-trimb-50,evtime1+parr+trima+50,pad=True,fill_value=0)                
                
                tr.detrend()
                tr.taper(tpratio)
                tr.data = bandpass(tr.data,frqmin,frqmax,tr.stats.sampling_rate,2,True)
                #tr.taper(tpratio)
                tr.trim(evtime1+parr-trimb,evtime1+parr+trima)
                tr.stats.starttime = 0
                npts = tr.stats.npts
                tr.normalize() 
                
                oridata = tr.data
                noisea = int(noise_a*rsample)
                noiseb = int(noise_b*rsample)
                siga = int(sig_a*rsample)
                sigb = int(sig_b*rsample)
                snr = max(abs(tr.data[noisea:noiseb]))/(1e-6+max(abs(tr.data[siga:sigb])))
                if tdomain == 1:
                    #<<<<<<< HEAD:main_new.py
                    #tr = autocorr_td(tr,windowb,windowa,fwin)
                    corr = autocorr_td(tr,windowb,windowa,fwin)
#=======
#                    tr = autocorr.autocorr_td(tr,windowb,windowa,fwin)
#>>>>>>> 8915670f06c8ef3611fd0d088d977e4b490d1a9b:main.py
                else:
                    tr = autocorr.autocorr_fd(tr,fwin)

                dist = dis1.delta
                #if np.isnan(tr.data).any():
                if np.isnan(corr).any():
                    continue

                evtag = '/'.join([stanet,staname,location+'_'+evname.split('.')[0]])
#                print(npts,len(tr.data))
                #ds.create_dataset('/'+evtag,(npts,),data=tr.data[:npts])
                npts = len(corr)
                #ds.create_dataset('/'+evtag,(npts,),data=tr.data[:npts])
                #ds.create_dataset('/'+evtag,(npts,),data=corr)
                savedata = np.zeros((npts,2))
                savedata[:,0] = corr
                savedata[:len(oridata),1] = oridata
                ds.create_dataset('/'+evtag,(npts,2),data=savedata)
                ds[evtag].attrs['snr'] = snr
                ds[evtag].attrs['evid'] = evname+'.'+tr.id
                ds[evtag].attrs['evlat'] = evlat1
                ds[evtag].attrs['evlon'] = evlon1
                ds[evtag].attrs['evdep'] = evdep1
                ds[evtag].attrs['stlat'] = stlat
                ds[evtag].attrs['stlon'] = stlon
                ds[evtag].attrs['dist'] = dis1.delta
                
    ds.close()
    return 0



from scipy import interpolate
import numpy as np
import pickle
import numpy as np
import obspy
import distaz
import os
from scipy.fftpack import fft,ifft
import glob
import pandas as pd
import autocorr
from obspy.signal.filter import bandpass
import h5py
import yaml
import time
import sys

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if len(sys.argv) == 1:
    parafile = 'parameters.in'
else:
    parafile = str(sys.argv[1])

with open(parafile,'r') as fin:
    par = yaml.load(fin)
    #<<<<<<< HEAD:main_new.py
    #par = yaml.load(fin,Loader=yaml.FullLoader)
#=======
    # py3
    #par = yaml.load(fin,Loader=yaml.FullLoader)
    # py2
#>>>>>>> 8915670f06c8ef3611fd0d088d977e4b490d1a9b:main.py
#    par = yaml.load(fin)

trimb   = par['trimb']
trima   = par['trima']
frqmin  = par['frqmin']
frqmax  = par['frqmax']
rsample = par['rsample']
tpratio = par['tpratio']
tdomain = par['timedomain']
datatp  = par['datatype']
mindist = par['mindist']
maxdist = par['maxdist']
fwin    = par['fwin']
eqdir   = par['eqdir']
outfile = par['outfilename']
comp = par['comp']

ndep = par['ndep']
ndis = par['ndis']
mindis = par['mindis']
maxdis = par['maxdis']
mindep = par['mindep']
maxdep = par['maxdep']
phasesl = par['phaselist']
firstarrtt = par['firstarrtt']
windowb = par['align_b']
windowa = par['align_a']
sig_b = par['sig_b']
sig_a = par['sig_a']
noise_b = par['noise_b']
noise_a = par['noise_a']
ncore = par['ncore']
tempdata = par['tempdata']
#phasesl = phasel.split(',')

#ndep = 96
#ndis = 111
dep = np.linspace(mindep,maxdep,ndep)
dis = np.linspace(mindis,maxdis,ndis)
filetele = './tables/'+firstarrtt
#filepcp = './tables/'+cmbreftt
if os.path.isfile(filetele):
    telep = np.load(filetele)
#    pcpdp = np.load(filepcp)
else:
    print('constructing tt table')
    from obspy.taup import TauPyModel
    mod = TauPyModel(model='ak135')
    telep = np.zeros((ndep,ndis))
#    pcpdp = np.zeros((ndep,ndis))
    for ii in range(ndep):
       for jj in range(ndis):        
            arr = mod.get_travel_times(source_depth_in_km=dep[ii],distance_in_degree=dis[jj],phase_list=phasesl)
            #pcpdp[ii,jj] = arr[-1].time-arr[0].time
            telep[ii,jj] = arr[0].time
    np.save(filetele,telep)
ftelep = interpolate.interp2d(dis, dep, telep , kind='linear')
#fpcpdp = interpolate.interp2d(dis, dep, pcpdp, kind='linear')

#pPdP = np.load('./tables/teleP50to600.npy')
#ndep = 551
#ndis = 651
#dep = np.linspace(50,600,ndep)
#dis = np.linspace(30,95,ndis)
#ftelep = interpolate.interp2d(dis, dep, pPdP, kind='linear')
#
#ndep = 121
#ndis = 101
#dep = np.linspace(0,600,ndep)
#dis = np.linspace(25,75,ndis)
#pPdP = np.load('./tables/PcPdP.npy')
#fpcpdp = interpolate.interp2d(dis, dep, pPdP, kind='linear')

eqs = open(''.join([eqdir,'/EVENTS-INFO/event_list_pickle']),'rb')
eqs = pickle.load(eqs)

evid = []
for ii in eqs:
    evid.append(ii['event_id'])
    
print ('begin stacking for each station')

evelist = glob.glob(eqdir+'/*.a')
nevent = len(evelist) 
#irissta = pd.read_csv('./tables/IRISSTA0319.txt',names=('net','sta','lat','lon'),header=0,delim_whitespace=True,keep_default_na=False)

samplevent = np.arange(0,nevent)
if comp == 'BHT' or comp == 'BHR':
    oricomp = 'BHE'
else:
    oricomp = comp

start = time.ctime()
if rank == 0:
    print('starting time:',start)
    if not os.path.isdir(tempdata):
        os.makedirs(tempdata)
ntime = int(np.ceil(float(nevent)/ncore))
for ievent in range(ntime):
    process_event(rank*ntime+ievent)
#for ievent in samplevent:
if rank==0:
    print ('finishing',time.ctime())
