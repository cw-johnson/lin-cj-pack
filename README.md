# lin-cj-pack
Linear Algebra Package by CJ


# Target Architectures
lin-cj-pack uses OS-specific CBLAS libraries. It currently supports...

MacOS (darwin) using Apple's "Accelerate" Framework

Linux using the generic CBLAS interface
    You can install architecture specific backends, i.e OpenBLAS, Intel MKL, AMD AOCL-BLIS for increased performance


# Prerequisites
For MacOS, ensure you have Apple Developer Tools installed. You can install them by running
'''
xcode-select --install
'''

For Linux ensure you have the CBLAS library installed. Using the apt package manager, install by running
'''
sudo apt-get install libblas-dev liblapack-dev
'''

Alternatively, run
'''

'''
# Compilation
Navigate to the repository directory, and type "make" to compile all demonstration programs


