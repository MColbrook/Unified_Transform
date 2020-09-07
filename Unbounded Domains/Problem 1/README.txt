This code computes the errors (as defined in the paper) for problem (a) using:

UT_method.m       the unified transform
Spec_method.m     the propose spectral method
Sep_Var_method.m  separation of variables

Some additional functions are required for quick evaluation of functions and collocation points. These are provided in the subfolder "Function_evals" which should be added to the path before using the above routines. The routines can easily be adapted to plot the approximate solution or for different problems (such as different boundary data).

Copyright (c) 2019 Matthew Colbrook