#/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
# File Name : main.py
#
# Purpose : ACC to extract PcP
#
# Creation Date : 05-07-2019
#
# Last Modified : Sun Jul  7 08:57:29 2019
#
# Created By : Hongjian Fang: hfang@mit.edu 
#
#_._._._._._._._._._._._._._._._._._._._._.*/

def autocorr_fd(tr,conlen=10):
        N = tr.stats.npts
        delta = tr.stats.delta
        data = np.zeros(2*N,)
        data[:N] = tr.data
        yf = fft(data)
        ayf = np.abs(yf)
        myf = np.convolve(ayf, np.ones((conlen,))/conlen, mode='same')
        yf = yf/myf
        myf = yf*np.conj(yf)
        xx = np.real(ifft(myf))
        tr.data = xx[:N]
        return tr


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
from obspy.signal.filter import bandpass
import h5py
import yaml
import time

with open('parameters.in','r') as fin:
    par = yaml.load(fin)

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
bins    = par['bins']
fwin    = par['fwin']
eqdir   = par['eqdir']
outfile = par['outfilename']

pPdP = np.load('./tables/teleP50to600.npy')
ndep = 551
ndis = 651
dep = np.linspace(50,600,ndep)
dis = np.linspace(30,95,ndis)
ftelep = interpolate.interp2d(dis, dep, pPdP, kind='linear')

ndep = 121
ndis = 101
dep = np.linspace(0,600,ndep)
dis = np.linspace(25,75,ndis)
pPdP = np.load('./tables/PcPdP.npy')
fpcpdp = interpolate.interp2d(dis, dep, pPdP, kind='linear')

eqs = open(''.join([eqdir,'/EVENTS-INFO/event_list_pickle']))

eqs = pickle.load(eqs)
evid = []
for ii in eqs:
    evid.append(ii['event_id'])
    
print ('begin stacking for each station')
ds = h5py.File(outfile,'w')

evelist = glob.glob(eqdir+'/*.a')
nevent = len(evelist) 
irissta = pd.read_table('./tables/IRISSTA0319.txt',names=('net','sta','lat','lon'),header=0,delim_whitespace=True,keep_default_na=False)

samplevent = np.arange(nevent)

start = time.time()
print('starting time:',start)
for ievent in samplevent:
    ievent1 = ievent
    evname = evelist[ievent1].split('/')[-1]

    evidx1 = evid.index(evname)
    evlat1 = eqs[evidx1]['latitude']
    evlon1 = eqs[evidx1]['longitude']
    evdep1 = eqs[evidx1]['depth']
    evtime1 = eqs[evidx1]['datetime']
    evmag1 = eqs[evidx1]['magnitude']  
    stalist1 = glob.glob('/'.join([eqdir,evname,datatp,'*.BHZ']))
    stalist = [stal.split('/')[-1] for stal in stalist1]

    nsta = len(stalist1)
    stapos = np.zeros((nsta,2))
    print ('the',ievent,'th event pairs')
    print ('evinfo for eq 1:',evname,evmag1,evlat1,evlon1,evdep1,nsta)

    for ista in range(nsta):
                stanet = stalist[ista].split('.')[0]
                staname = stalist[ista].split('.')[1]
                location = stalist[ista].split('.')[2]
                stasub = irissta[(irissta['net']==stanet) & (irissta['sta']==staname)]
                stlat = stasub.iloc[0]['lat']
                stlon = stasub.iloc[0]['lon']
                
                dis1 = distaz.DistAz(evlat1,evlon1,stlat,stlon)

                if dis1.delta<mindist or dis1.delta>maxdist:
                    continue

                trace = '/'.join([eqdir,evname,datatp,stalist[ista]])
                if os.path.getsize(trace) < 1000:
                    continue
                strm = obspy.read(trace)
                tr = strm[0]
                if tr.stats.npts<100:
                    continue
                tr.resample(rsample)
                
                parr = ftelep(dis1.delta,evdep1)[0]
                tr.stats.distance = dis1.delta
                tr.trim(evtime1+parr-trimb-50,evtime1+parr+trima+50,pad=True,fill_value=0)                
                
                tr.detrend()
                tr.taper(tpratio)
                tr.data = bandpass(tr.data,frqmin,frqmax,tr.stats.sampling_rate,2,True)
                tr.taper(tpratio)
                tr.trim(evtime1+parr-trimb,evtime1+parr+trima)
                tr.stats.starttime = 0
                npts = tr.stats.npts
                tr.normalize() 
                
                tr = autocorr_fd(tr,fwin)

                dist = dis1.delta
                if np.isnan(tr.data).any():
                    continue

                evtag = '/'.join([stanet,staname,location+'_'+evname.split('.')[0]])
                ds.create_dataset('/'+evtag,(npts,),data=tr.data)
                ds[evtag].attrs['evid'] = evname+'.'+tr.id
                ds[evtag].attrs['evlat'] = evlat1
                ds[evtag].attrs['evlon'] = evlon1
                ds[evtag].attrs['evdep'] = evdep1
                ds[evtag].attrs['stlat'] = stlat
                ds[evtag].attrs['stlon'] = stlon
                ds[evtag].attrs['dist'] = dis1.delta
                
print 'finishing',time.time()
ds.close()

