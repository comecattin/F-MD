FC = gfortran
FFLAGS ?= -O2
OBJ = units.o input_parser.o initialization.o forces.o contraints.o integration.o output_module.o  energies.o main.o

md_simulation: $(OBJ)
	$(FC) $(FFLAGS) -o md_simulation $(OBJ)

units.o: units.f90
	$(FC) $(FFLAGS) -c units.f90

main.o: main.f90 input_parser.o initialization.o forces.o integration.o output_module.o energies.o
	$(FC) $(FFLAGS) -c main.f90

initialization.o: initialization.f90
	$(FC) $(FFLAGS) -c initialization.f90

forces.o: forces.f90 units.o
	$(FC) $(FFLAGS) -c forces.f90

integration.o: integration.f90 forces.o units.o
	$(FC) $(FFLAGS) -c integration.f90

energies.o: energies.f90 units.o
	$(FC) $(FFLAGS) -c energies.f90

output_module.o: output_module.f90
	$(FC) $(FFLAGS) -c output_module.f90

contraints.o: contraints.f90
	$(FC) $(FFLAGS) -c contraints.f90

input_parser.o: input_parser.f90
	$(FC) $(FFLAGS) -c input_parser.f90

clean:
	rm -f *.o *.mod md_simulation trajectories.dat


.PHONY: clean
