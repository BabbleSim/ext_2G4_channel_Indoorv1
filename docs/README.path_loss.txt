*************
* PATH LOSS *
*************

For the path loss two options are possible:
* A constant attenuation value, equal for all paths (default behavior):
  Can be set with -at=<attenuation> as a dB value
  
* A distance for each path.
  From this distance the path loss attenuation will be calculated. (The distance exponent can be specified with -exp=<distance_exp>)
  (see 2G4Indoorv1_Description.pdf for information on how this is calculated)

  The distances for each path is expecified with:
  -df=<distance_matrix_file>
  The distance for each of the paths may be either a constant (specified directly in the file)
  or a reference to another file <distance_file> which can contain either a constant or a set of distances over time
  For more information on these files' formats see below

On top of this, a delta can be added as -atextra=<at_extra>
This value will be added to all paths path loss (be it specified with -at or from files).

  * The <distance_matrix_file>: *

    For each valid path this file shall contain one line, of the format
    x y : {value|"<distance_file>"}
    
    Where
    * x is the transmitter number, x = 0..N-1
    * y is the receiver number     y != x, y > x
    * N = number of devices
    * value = a number of meters as a floating point value applying to that path (and the reciprocal)
    * "<distance_file>" : is the path to a <distance_file> as described below which applies to that path (and its reciprocal) 
      Note that this file name shall be provided in between "" 

  * The <distance_file>: *

    This file contains 2 columns, the first column represents time, and the second column the distance at that given time
    The times shall be in ascending order.
    If the simulation time is found in the file, the channel will use the corresponding distance value.
    If the simulation time falls in between 2 times in the file, the channel will interpolate accordingly the distance.
    If the simulation time is before the first time value in the file, the channel will use the first distance provided in the file.   
    If the simulation time is after the last time value in the file, the channel will use the last distance provided in the file.

The file names can be either absolute or relative to bin/  

For both files, '#' is treated as a comment mark (anything after a # is discarded)
Empty lines are ignored


These files can be generated automatically with utilities/distance_files_from_coordinates/create_distance_files.m
See the examples provided in that same folder for more info
