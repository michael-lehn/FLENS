## FLENS: Flexible Library for Efficient Numerical Solutions


FLENS is a C++ library for the efficient implementation of numerical algorithms.

Compiles with

    clang++
    g++ (>= 4.7)
    icc (>= 12.1.2)

Read more on

    http://www.mathematik.uni-ulm.de/~lehn/FLENS/


If you want to test FLENS-LAPACK with Netlib's LAPACK test suites
simply run

    make check

Note that this requires a Fortran compiler.

## Project Status

FLENS is no longer actively developed.

In its time, FLENS had what I would describe as a small to moderate success. It was used in academic research, including PhD projects, found its way into other scientific software, and occasionally generated discussions in the C++ and scientific-computing communities. I was very happy to see that people outside our group found the library useful.

What I never managed to achieve, however, was to build a community of developers around FLENS. While there were contributors, the long-term development and maintenance of the library remained largely my responsibility.

Over time, a significant part of maintaining FLENS therefore became user support. I regularly received questions from researchers using it for their PhD work, often in areas such as numerical methods for PDEs. Many of these questions were not really about FLENS itself. They concerned which numerical method should be used, how an algorithm should be implemented, or sometimes more fundamental issues with C++.

Helping with such questions could be interesting and rewarding, but incorporating the resulting fixes, extensions, examples, and explanations into FLENS week after week was not sustainable as a one-person project. At some point I therefore stopped active development. Without other developers taking ownership of parts of the project, I no longer saw a good reason to keep expanding and maintaining the library indefinitely.

The repository remains online because I still think some of the ideas behind FLENS are interesting, and because the code and documentation may be useful as a historical reference for approaches to expressive, high-performance numerical linear algebra in C++.
