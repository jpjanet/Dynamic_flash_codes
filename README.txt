These codes simulate dynamic flash vessels using the equations given in [1]. 
There are two sets of equations. The first set, "System 1", comprises of the
full set of equations. By contrast, "System 2" is a simplification based on the 
assumption of negligible vapour holdup.

For more detailed information, please consult the attached thesis, particularly
Appendix C, which provides more information on these codes. A short instructive  
summary is given here.

The codes are provided in a package that contains three sub-directories: common,
system1 and system2. The only file that needs to be edited and invoked is 
dynamic_flash_mainfile.m, which lives in the top level directory. The code 
supports two methods of operation. Firstly, the code can be provided with a 
single tolerance and test case, in which case the code executes a single 
simulation using both sets of equations - first System 1, and then 
System 2 are solved for the given parameters. 

Alternatively, if both tol is a vector and clock is set on, the code will 
detect this and switch into an alternate behaviour. If these two conditions 
are fulfilled, the code will instead sweep through the simulation for each 
tolerance, both solvers and each test case. 

As Test Case 3 takes longer to complete, the third test case can be excluded 
by setting the parameter include3 to 0. Effectively, this mode compares the 
solution time across all combinations of test case, solver and system for 
the vector of tolerances specified. 

More details about the inputs and outputs are given in the associated report.

These codes were written and tested in Matlab® 2014b, and compatibility with 
other versions has not been tested.

References
 [1]: Biegler, L. T. (2010). 
 "Nonlinear programming: concepts, algorithms, and applications to
  chemical processes, volume 10. SIAM, Philadelphia, PA.

