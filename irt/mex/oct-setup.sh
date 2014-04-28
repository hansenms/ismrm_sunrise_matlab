# oct-setup.sh
# extra flags needed for compiling mex files for octave

octlibdir = /opt/local/lib/octave/3.6.4/
octincdir = /opt/local/include/octave-3.6.4/octave

# http://stackoverflow.com/questions/7806418/using-setenv-in-makefile
export XTRA_CFLAGS=-std=c99 -UCountAlloc -DMmex -DUse_simd -DUse_thread -O3
