Overview
--------

This is a program that computes approximations to `float -> float` functions.
You give it a straight-line program in a certain form with some unspecified
constants, an interval, and a function giving, for each `x` in that interval,
the acceptable values for the function; it tries to fill in the blanks in your
straight-line program in order to get an acceptable approximation to the
function over the whole interval of interest.

Building
--------

This program depends on [`QSopt_ex` by Applegate, Cook, Dash, and
Espinoza.](http://www.math.uwaterloo.ca/~bico/qsopt/ex/).  This program also
depends on a complete and working installation of [GMP](https://gmplib.org/).
Download and build `QSopt_ex` first in directory `QSex090408`, then run `make`.

Using
-----

Hack `main()` to specify the straight-line program of interest, bounds on the
unspecified coefficients, and the interval.  Hack `get_bounds()` to specify the
function and error bounds of interest.  Then build and run the program.

The output
----------

The output of this program is a mishmash of various kinds of status information.

A chunk that looks like this:

    float c0 = 0x1.0ce86cp-4
    float c1 = -0x1.81a69cp-4
    float c2 = 0x1.07661ep-4
    float c3 = 0x1.f6d1fp-8
    float c4 = 0x1.c6cc72p-5
    float c5 = 0x1.10f968p-3
    float c6 = 0x1.55550ap-2

represents the current candidate list of coefficients.

A chunk that looks like this:

    vector<float> testpoints = {0x1p-1,0x1.bf4876p-2,0x1.47d9bep-1,0x1.f0b3c2p-2,0x1.4e15b4p-1,0x1.2915f8p-1,0x1.36870ep-3,0x1.67b29ap-2,0x1.7bb93cp-1,0x1.762f4cp-1};

represents the current list of test points.

Lines that look like this:

    0x1.7cp+10
    0x1.a724p+16
    0x1.3670a8p+22
    0x1.a3342ep+24
    0x1.1e7b58p+25

represent the current status of the diving search.  Each line gives the number of `float`s for some coefficient that have not yet been proven infeasible.

At the end of execution, the program prints out `0`, `-1`, or `-2`.  `0` means
that the most recent candidate list of coefficients provably works.  `-1` means
that there are provably no coefficients that work.  `-2` means that the
heuristic gave up trying to solve the problem because a recursive search took too long.
