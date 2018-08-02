# Heat Transport in One Dimension

## About

While heat transport at large scales follows Fouriers law, at the level of nanoscale
systems, heat transport becomes anomolous. Defining rigorously even what is meant
by temperature of a system at these scales away from equilibrium is non-trivial.

The goal of my thesis is to investigate the anomaly of heat transport in one dimension.
In doing so it became necessary to investigate various systems numerically in order 
to understand the qualitative impact of different approximations.

There main numerical experiments of my are:

* A non-Markovian system driven by non-stationary Gaussian noise

* A one parameter family non-Markovian system driven by stationary Gaussian noise,
which limit to a Langevin equation

* A Langevin Equation

Only the first two of the above systems are included in this repo.

## Running this code

If in Linux, simply run the Makefile in either the directory

/dev/non_stationary_v_1.5

or

/dev/stationary_v_1.9.

Running the code results dumps several csv files. For large systems with long run
times these files can be up to gigs in size. The run time for these systems
can be upwards of three days on 8 cores.

To analyze the data in the csv. Place the csv's files in a folder
/dev/analysis_of_observables/*system_size*
where *system_size* is the numerical representation of the system size of the run.
Then in Matlab run process_raw_data.m and then analyze_processed_data.m.
You may have to change some parameters in process_raw_data.m in order depending on the
system size you have chosen. 


**This code requires a c++17 compliant compiler.** I use the built in Bessel functions
which are now built into to the standard library as of c++17. I use the GNU 7.2 compiler. 
If you cannot use that compiler for whatever reason, simply replace the Bessel functions 
with ones from a publicly available library, like the Boost libraries.