This is the original readme from the repo this code is based on,
which is by Nicholas Tessore

The code is based on a fast routine in C to compute the Wigner d-matrix for a
range of angular momentum values l and fixed m. This is done via a recurrence
in l, which is described at http://www.jb.man.ac.uk/~tessore/notes/180831.html.
To use the routine, simply copy the `wigner_d.c` source file into your project,
and let your code know about the new function by putting

    extern void wigner_d(int l0, int l1, int m, int k, double x, double* d);

near the top of your source.

An additional program to print the `wigner_d` function values is also included.

    usage: showdl l0 l1 m k theta
    
    Compute the Wigner d-matrix elements `d^l_{m,k}(theta)` for angular
    momentum l = l0, .., l1 and fixed values m, k, theta.
