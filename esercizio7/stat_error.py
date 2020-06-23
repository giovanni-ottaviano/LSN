#!/usr/bin/env python3

import numpy as np

#mean computed in range [min, max]
def ave_in_range(v, imin, imax):
    assert imax < len(v)
    assert imin < imax

    return np.mean(v[imin:imax])    #guarda se imax-1

#compute mean for each block an the square, the return them
#L = steps per block
def mean_blocks(v, L):
    assert len(v) % L == 0

    nblocks = len(v)//L
    ave = np.zeros(nblocks)

    for i in range(nblocks):
        ave[i] = np.mean(v[i*L:(i+1)*L])   #compute mean for each block

    ave2 = ave**2   #compute squares

    mean = np.mean(ave)
    mean2 = np.mean(ave2)

    return mean, mean2   #sono interessato solo al valore finale

#compute statistical error
def error(val, val2, n):
    if n == 1:
        return 0
    else:
        return np.sqrt((val2 - val**2)/(n-1.))


#*************************************MAIN**************************************
def main():
    L_begin, L_end = 10, 5 * 10**3
    L_incr = 1
    error_Epot, error_Pres = [], []
    L_list = []

    Epot = np.loadtxt('output.insta.epot.0', dtype=float)
    Pres = np.loadtxt('output.insta.pres.0', dtype=float)
    M = len(Epot)   #len(Epot) = len(Pres) ... number of MC steps

    for L in range(L_begin, L_end+1, L_incr):
        if M % L != 0:
            continue

        m_Epot, m2_Epot = mean_blocks(Epot, L)
        m_Pres, m2_Pres = mean_blocks(Pres, L)

        error_Epot.append(error(m_Epot, m2_Epot, M/L))
        error_Pres.append(error(m_Pres, m2_Pres, M/L))

        L_list.append(L)

    #print estimated functions
    np.savetxt('error_epot.dat', error_Epot, fmt='%.7f')
    np.savetxt('error_pres.dat', error_Pres, fmt='%.7f')
    np.savetxt('L.dat', L_list, fmt='%u')

    return 0



#RUN SCRIPT
if __name__ == "__main__":
    main()
