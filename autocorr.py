#/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
# File Name : autocorr.py
#
# Purpose :
#
# Creation Date : 18-07-2019
#
# Last Modified : Thu Jul 18 09:11:18 2019
#
# Created By : Hongjian Fang: hfang@mit.edu 
#
#_._._._._._._._._._._._._._._._._._._._._.*/
import numpy as np
from scipy.fftpack import fft,ifft
def autocorr_td(tr,stalign_a,stalign_b,conlen=15):
        N = tr.stats.npts
        delta = tr.stats.delta
        yf = fft(tr.data)
        ayf = np.abs(yf)
        myf = np.convolve(ayf, np.ones((conlen,))/conlen, mode='same')
        yf = yf/myf
        xx = np.real(ifft(yf))
        
        stalign_a = int(stalign_a/delta)
        stalign_b = int(stalign_b/delta)
        winlen = stalign_b - stalign_a
        xxw = xx[stalign_a:stalign_b]
        acorr = np.correlate(xx,xxw,mode='full')
        acorr = acorr[winlen:]
        maxloc = np.argmax(abs(acorr))
        acorr = np.roll(acorr,-maxloc)
        
        tr.data[:len(acorr)] = acorr
        return tr

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



