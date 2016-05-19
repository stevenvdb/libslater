# libslater
Library for the analytical evaluation of integrals involving Slater densities

libslater is a library for the analytical evalution of integrals involving Slater densities.
The code can be automatically generated using SymPy by running the script `gencode.py' in the folder `generate'. Source code files will be placed in the `src' folder. Once this is done, you should be able to install the library using something like
    python setup.py config_fc  --f90flags="-fno-range-check" build_src build_ext
    python setup.py install --home=~
The `gencode.py' script contains one important parameter: mmax. This indicates the highest value for m that will be considered, where m is the power of r in the Slater density:
    rho = 1/(4 pi (m+2)! a**(m+3)) r^m  exp(-r/a)
Note that even for mmax=3 the code generation takes a while: SymPy is very slow for computing some complicated expressions.
For convenience, the source code for mmax=3 is provided in src_mmax3; by copying this folder to src, you can skip the code generation step.
