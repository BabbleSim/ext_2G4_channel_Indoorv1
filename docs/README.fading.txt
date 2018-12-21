************
*  FADING  *
************

It implements a model of the multipath / fast fading as described in:
2G4Indoorv1_Description.pdf
See the code for more information about the implementation.

Parameters: (for a longer description see the documentation)
    [-ds=<delay_spread>]   :RMS delay spreed of the channel (in ns) (default 15ns)
    [-speed=<doppler_speed]:Speed of the transmitter receiver and/or enviroment 
                            in m/s (default 1.1m/s = 4km/h = slow walk)
    [-RiceK=<riceK>]       :Ricean K factor: ratio of the line of sight power to
                            the scattered power.
                            Default 2, expected range from 0 to 3 = No line of
                            sight to approx 1.5m distance between tx and rx in a
                            small room

