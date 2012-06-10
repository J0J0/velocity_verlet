# vim: noet

.PHONY : all clean cleanall underscorify 

PROGRAMS = velocity_verlet_algorithm

FC       = gfortran
GDEBUG   = -g
OMP      = -fopenmp
OPTIMIZE = -O0
FOTHER   =
FFLAGS   = $(GDEBUG) -fdefault-real-8 $(OPTIMIZE) -fimplicit-none $(FOTHER)  # $(OMP)

all : $(PROGRAMS) underscorify

%: %.o
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.F90
	$(FC) $(FFLAGS) -c $<


underscorify:
	@for f in $(PROGRAMS) ; do if test -e $$f; then echo mv $$f _$$f; mv $$f _$$f; fi; done

plot:
	./gen_frames.sh
	./gen_plot_script.sh
	gnuplot plot_particle_movement.gpl

clean:
	rm -f *.o *.mod

cleanoutput:
	rm -f output/particle_*.txt
	rm -f output/frames.txt output/frames_info.txt
	rm -f output/energy.txt

cleanall: clean cleanoutput
	rm -f $(PROGRAMS)
	rm -f $(addprefix _,$(PROGRAMS))
	rm -f plot_particle_movement.gpl
	rm -f particle_movement.gif

