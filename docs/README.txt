Model of the 2.4GHz Indoor channel
To be used as a channel library with BabbleSim, compatible with the 2.4GHz Phy v1

Note that this model requires the fftw3 library installed in your system.
In Debian/Ubuntu install libfftw3-dev (sudo apt-get install libfftw3-dev)

For information on its command line arguments call it with --help , i.e.
  ./bs_2G4_phy_v1 -s=banana  -D=2 -channel=Indoorv1 -argschannel --help

This channel models the fading as described in
2G4Indoorv1_Description.pdf
For more information about the details about the fading implementation refer to docs/README.fading.txt
and to the source code comments: 2G4_channel_Indoorv1.c

For more information about the path_loss part please refer to docs/README.path_loss.txt
