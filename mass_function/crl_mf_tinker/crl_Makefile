FC	= g95
FFLAGS	= -O3 -freal-loops
OBJS	= cosmo.o linearpk.o sigma.o spline.o hunt.o qromb.o chebyshev.o growth.o dverk.o mf_tinker.o
.SUFFIXES: .f90
.f90.o: 
	$(FC) $(FFLAGS) -c $<
default:  sample compute_mf_tinker compute_mf_tinker2
compute_mf_tinker: $(OBJS) compute_mf_tinker.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)
compute_mf_tinker2: $(OBJS) compute_mf_tinker2.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)
sample: $(OBJS) sample.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)
clean:
	-rm -f *.o *.mod
tidy: clean
	-rm -f sample
	-rm -f compute_mf_tinker
	-rm -f compute_mf_tinker2
	-rm -f Rh_dndlnRh.txt
	-rm -f Mh_dndlnMh.txt
	-rm -f Mhom0_dndlnMh.txt
