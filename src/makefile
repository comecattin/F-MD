FC = gfortran
OBJ = input_parser.o initialization.o forces.o contraints.o integration.o output_module.o  energies.o main.o

md_simulation: $(OBJ)
	$(FC) -o md_simulation $(OBJ)

main.o: main.f90
	$(FC) -c main.f90

initialization.o: initialization.f90
	$(FC) -c initialization.f90

forces.o: forces.f90
	$(FC) -c forces.f90

integration.o: integration.f90
	$(FC) -c integration.f90

energies.o: energies.f90
	$(FC) -c energies.f90

output_module.o: output_module.f90
	$(FC) -c output_module.f90

contraints.o: contraints.f90
	$(FC) -c contraints.f90

input_parser.o: input_parser.f90
	$(FC) -c input_parser.f90

clean:
	rm -f *.o *.mod md_simulation trajectories.dat

