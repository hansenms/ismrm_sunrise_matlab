MPATH = /usr/local/matlab/bin/
MEX = $(MPATH)mex
MEXEXT = $(shell $(MPATH)mexext)
dir = ../mex

all:	\
	$(dir)/BsplExpand.$(MEXEXT) \
	$(dir)/BsplReduce.$(MEXEXT) \
	$(dir)/BsplVal2CoMirr.$(MEXEXT) \
	$(dir)/BsplCo2ValMirr.$(MEXEXT) \
	$(dir)/BsplCo2GdXMirr.$(MEXEXT) \
	$(dir)/BsplCo2GdYMirr.$(MEXEXT) \
	$(dir)/BsplCo2GdZMirr.$(MEXEXT) \
	$(dir)/BsplCo2ValTranMirr.$(MEXEXT) \
	$(dir)/BsplCo2GdXTranMirr.$(MEXEXT) \
	$(dir)/BsplCo2GdYTranMirr.$(MEXEXT) \
	$(dir)/BsplCo2GdZTranMirr.$(MEXEXT) \
	$(dir)/BsplVal2CoZero.$(MEXEXT) \
	$(dir)/BsplCo2ValZero.$(MEXEXT) \
	$(dir)/BsplCo2GdXZero.$(MEXEXT) \
	$(dir)/BsplCo2GdYZero.$(MEXEXT) \
	$(dir)/BsplCo2GdZZero.$(MEXEXT) \
	$(dir)/BsplCo2ValTranZero.$(MEXEXT) \
	$(dir)/BsplCo2GdXTranZero.$(MEXEXT) \
	$(dir)/BsplCo2GdYTranZero.$(MEXEXT) \
	$(dir)/BsplCo2GdZTranZero.$(MEXEXT) \
	$(dir)/BsplCo2ValZeroFilt.$(MEXEXT) \
	$(dir)/BsplCo2ValTranZeroFilt.$(MEXEXT)

PYRAM = BIG/pyramids/messagedisplay.c BIG/pyramids/pyramidtools.c BIG/pyramids/pyramidfilters.c BIG/pyramidtools3D.c

$(dir)/BsplExpand.$(MEXEXT): BsplExpand.c $(PYRAM)
	$(MEX) -o $@ BsplExpand.c $(PYRAM)

$(dir)/BsplReduce.$(MEXEXT): BsplReduce.c $(PYRAM)
	$(MEX) -o $@ BsplReduce.c $(PYRAM)

INTERP = batchinterpol2D.c batchinterpol3D.c interpol2D.c interpol3D.c

$(dir)/BsplVal2CoMirr.$(MEXEXT): BsplVal2CoMirr.c $(INTERP)
	$(MEX) -o $@ BsplVal2CoMirr.c $(INTERP)

$(dir)/BsplCo2ValMirr.$(MEXEXT): BsplCo2ValMirr.c $(INTERP)
	$(MEX) -o $@ BsplCo2ValMirr.c $(INTERP)

$(dir)/BsplCo2GdXMirr.$(MEXEXT): BsplCo2GdXMirr.c $(INTERP)
	$(MEX) -o $@ BsplCo2GdXMirr.c $(INTERP)

$(dir)/BsplCo2GdYMirr.$(MEXEXT): BsplCo2GdYMirr.c $(INTERP)
	$(MEX) -o $@ BsplCo2GdYMirr.c $(INTERP)

$(dir)/BsplCo2GdZMirr.$(MEXEXT): BsplCo2GdZMirr.c $(INTERP)
	$(MEX) -o $@ BsplCo2GdZMirr.c $(INTERP)

TINTERP =  batchinterpol2Dtran.c batchinterpol3Dtran.c interpol2Dtran.c interpol3Dtran.c

$(dir)/BsplCo2ValTranMirr.$(MEXEXT): BsplCo2ValTranMirr.c $(TINTERP)
	$(MEX) -o $@ BsplCo2ValTranMirr.c $(TINTERP)

$(dir)/BsplCo2GdXTranMirr.$(MEXEXT): BsplCo2GdXTranMirr.c $(TINTERP)
	$(MEX) -o $@ BsplCo2GdXTranMirr.c $(TINTERP)

$(dir)/BsplCo2GdYTranMirr.$(MEXEXT): BsplCo2GdYTranMirr.c $(TINTERP)
	$(MEX) -o $@ BsplCo2GdYTranMirr.c $(TINTERP)

$(dir)/BsplCo2GdZTranMirr.$(MEXEXT): BsplCo2GdZTranMirr.c $(TINTERP)
	$(MEX) -o $@ BsplCo2GdZTranMirr.c $(TINTERP)


$(dir)/BsplVal2CoZero.$(MEXEXT): BsplVal2CoZero.c
	$(MEX) -o $@ BsplVal2CoZero.c

$(dir)/BsplCo2ValZero.$(MEXEXT): BsplCo2ValZero.c $(INTERP)
	$(MEX) -o $@ BsplCo2ValZero.c $(INTERP)

$(dir)/BsplCo2GdXZero.$(MEXEXT): BsplCo2GdXZero.c $(INTERP)
	$(MEX) -o $@ BsplCo2GdXZero.c $(INTERP)

$(dir)/BsplCo2GdYZero.$(MEXEXT): BsplCo2GdXZero.c $(INTERP)
	$(MEX) -o $@ BsplCo2GdYZero.c $(INTERP)

$(dir)/BsplCo2GdZZero.$(MEXEXT): BsplCo2GdZZero.c $(INTERP)
	$(MEX) -o $@ BsplCo2GdZZero.c $(INTERP)

$(dir)/BsplCo2ValTranZero.$(MEXEXT): BsplCo2ValTranZero.c $(TINTERP)
	$(MEX) -o $@ BsplCo2ValTranZero.c $(TINTERP)

$(dir)/BsplCo2GdXTranZero.$(MEXEXT): BsplCo2GdXTranZero.c $(TINTERP)
	$(MEX) -o $@ BsplCo2GdXTranZero.c $(TINTERP)

$(dir)/BsplCo2GdYTranZero.$(MEXEXT): BsplCo2GdYTranZero.c $(TINTERP)
	$(MEX) -o $@ BsplCo2GdYTranZero.c $(TINTERP)

$(dir)/BsplCo2GdZTranZero.$(MEXEXT): BsplCo2GdZTranZero.c $(TINTERP)
	$(MEX) -o $@ BsplCo2GdZTranZero.c $(TINTERP)

FINTERP = BIG/BsplnTrfKer.c BIG/reg/convolve.c BIG/reg/getPut.c BIG/phil.c

$(dir)/BsplCo2ValZeroFilt.$(MEXEXT): BsplCo2ValZeroFilt.c filtinterpol2D.c filtinterpol3D.c $(FINTERP)
	$(MEX) -o $@ BsplCo2ValZeroFilt.c filtinterpol2D.c filtinterpol3D.c $(FINTERP)

$(dir)/BsplCo2ValTranZeroFilt.$(MEXEXT): BsplCo2ValTranZeroFilt.c filtinterpol2Dtran.c filtinterpol3Dtran.c $(FINTERP)
	$(MEX) -o $@ BsplCo2ValTranZeroFilt.c filtinterpol2Dtran.c filtinterpol3Dtran.c $(FINTERP)

reset:
	rm -i ../mex/*.$(MEXEXT)
