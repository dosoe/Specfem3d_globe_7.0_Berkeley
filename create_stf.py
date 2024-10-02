'''
Author: Utpal Kumar, BSL, UCB
Date: 2023-06-28
Email: utpalkumar@berkeley.edu
'''
import numpy as np
import matplotlib.pyplot as plt
from heaviside_plot_support.costap_filter import costap_filter
from heaviside_plot_support.plot_spectrum import plot_spectrum

def read_freq_band_par_file(parfilename):
    with open(parfilename, 'r') as file:
        lines = file.readlines()

    values = []

    # Extract the numeric values from the desired lines
    for line in lines:
        line = line.strip()  # Remove leading/trailing whitespace
        if line.startswith('SOURCE_T'):
            parts = line.split('=')
            if len(parts) >= 2:
                value = parts[1].split('.')[0].strip()
                try:
                    values.append(int(value))
                except ValueError:
                    continue

    # Print the extracted values
    # print(values)
    return values


def main(freqList):
    # f4 = 1/10 
    # f3 = 1/12.5
    # f2 = 1/30 
    # f1 = 1/50

    f4 = 1/freqList[3] 
    f3 = 1/freqList[2]
    f2 = 1/freqList[1] 
    f1 = 1/freqList[0]

    dt = 0.4
    time_shift = 500

    x = np.arange(0, 3601, dt)
    ind = np.argmin(np.abs(x - time_shift))
    heavi = np.zeros(len(x))
    heavi[ind] = 1/dt

    # Assuming costap_filter is a function defined elsewhere in your code
    stf = costap_filter(heavi, f1, f2, f3, f4, dt)

    # Assuming plot_spectrum is a function defined elsewhere in your code
    a, b, c = plot_spectrum(stf, dt)

    # fig, ax = plt.subplots()
    # ax.plot(x, heavi)
    # plt.savefig('heaviside.png', dpi=300, bbox_inches='tight')
    # plt.close()

    # fig, ax = plt.subplots()
    # ax.plot(x, stf)
    # plt.savefig('stf.png', dpi=300, bbox_inches='tight')
    # plt.close()

    fig, ax = plt.subplots()
    c = c + 1e-10 # to avoid divide by zero error
    ax.plot(1/c, b)
    ax.set_xlim([freqList[3]-10, freqList[0]+10])
    ax.set_xlabel('Period (s)')

    # Add vertical lines at each value in freqList
    for freq in freqList:
        ax.axvline(x=freq, color='red', linestyle='--')
        ax.text(freq, ax.get_ylim()[1], str(freq), ha='center', va='bottom')
    
    plt.savefig('stf_spectrum.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    freqList = read_freq_band_par_file(parfilename = 'DATA/Par_file')
    main(freqList)