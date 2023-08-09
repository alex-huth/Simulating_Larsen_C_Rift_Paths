Compilation instructions:

Modify the makefile so that the NETCDFHOME variable is set to your local netcdf home directory

compile the CD-MPM-SSA code:
`make`

compile the user functions:
`elmerf90 USF_iGimp.F90 -o USF_iGimp`