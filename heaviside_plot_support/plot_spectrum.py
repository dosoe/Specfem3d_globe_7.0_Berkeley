import numpy as np

def plot_spectrum(data, dt):
    ndata = 2**np.ceil(np.log2(len(data)))
    ndata = int(ndata)  # numpy gives float so we convert it to int
    
    # Zero-padding is done automatically by the fft function in numpy when 
    # the size of the input is less than the fft size, hence the commented line in MATLAB script is not needed
    # data = np.pad(data, (0, int(ndata-len(data))))
    
    dataf = np.fft.fft(data, ndata)
    dataf = dataf[:ndata//2]

    freqs = np.fft.fftfreq(ndata, dt)
    freqs = freqs[:ndata//2]

    amp = np.abs(dataf)
    amp = amp[:ndata//2]

    return dataf, amp, freqs
