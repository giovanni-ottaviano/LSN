#!/usr/bin/env python3

import numpy as np

#mean value of the product at different times (t and t + t_prime)
def mean_prod(obs_array, t):
    sum = 0.
    last = len(obs_array) - t

#compute product at different times and return mean value
    for t_prime in range(last):
        sum += obs_array[t_prime]*obs_array[t_prime + t]

    return sum/last

#restricted mean - necessary for computing correctly AC
def restricted_mean(obs_array, index):
    if index >= 0:
        rst_array = obs_array[index:]
    else:
        last = len(obs_array) + index   #index is negative so len() - abs(index)
        rst_array = obs_array[:last]

    return np.mean(rst_array, dtype=np.float64)


#Time-displaced autocorrelation function
def AC(obs_array, time):
    mean_tprime = restricted_mean(obs_array, -time)
    mean_t_plus_tprime = restricted_mean(obs_array, time)
    mean_td = mean_prod(obs_array, time)
    sigma_2 = np.var(obs_array, dtype=np.float64, ddof=0)   #perhaps change delta degree of freedom = 1, one dof lost 'cos mean

    return (mean_td - mean_tprime * mean_t_plus_tprime)/sigma_2


#************************************MAIN***************************************
def main():
    t_end = 2500 #max value is len - 1 
    Epot = np.loadtxt('output.insta.epot.solid', dtype=float)
    Pres = np.loadtxt('output.insta.pres.solid', dtype=float)

    print("loading completed...")

    #calculate the autocorrelation function for the 2 observables
    AC_epot = np.array( [AC(Epot, t) for t in range(t_end)] )
    print("end first AC...")
    AC_pres = np.array( [AC(Pres, t) for t in range(t_end)] )

    print("computing completed...")

    #print estimated functions
    np.savetxt('ac_epot_solid.dat', AC_epot, fmt='%.7f')
    np.savetxt('ac_pres_solid.dat', AC_pres, fmt='%.7f')

    return 0



#RUN SCRIPT
if __name__ == "__main__":
    main()
