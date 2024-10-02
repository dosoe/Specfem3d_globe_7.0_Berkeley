import numpy as np

def costap_filter(data, w1, w2, w3, w4, dt):
    ndata = 2**np.ceil(np.log2(len(data)) + 1)
    ndata = int(ndata)

    filt = np.zeros(ndata, dtype=np.complex) # change the array to complex type

    y = np.fft.fft(data, ndata)
    x = np.fft.fftfreq(ndata, dt)

    ind1 = np.where((x > w1) & (x < w2))
    filt[ind1] = y[ind1] * 0.5 * (1 + np.cos(np.pi * (w2 - x[ind1]) / (w2 - w1)))

    ind2 = np.where((x > w3) & (x < w4))
    filt[ind2] = y[ind2] * 0.5 * (1 - np.cos(np.pi * (w4 - x[ind2]) / (w4 - w3)))

    ind3 = np.where((x > w2) & (x < w3))
    filt[ind3] = y[ind3]

    filt[ndata//2:] = np.conj(filt[:ndata//2][::-1])  # corrected line

    filt = np.real(np.fft.ifft(filt))
    filt = filt[:len(data)]

    return filt
