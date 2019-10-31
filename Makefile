.PHONY: main
main: compute_energy_ewald.o global_variables.o main.o \
	input_output.o initialize_update.o
	
	gfortran -g -Wall -o main compute_energy_ewald.o global_variables.o \
	main.o input_output.o initialize_update.o -lfftw3

# .PHONY: main.o
main.o: initialize_update.mod global_variables.mod \
	input_output.mod compute_energy_ewald.mod main.f90
	
	gfortran -g -c main.f90
	
# .PHONY: initialize_update.o initialize_update.mod
initialize_update.o initialize_update.mod: compute_energy_ewald.mod \
  global_variables.mod input_output.mod initialize_update.f90
	
	gfortran -g -c initialize_update.f90
	
# .PHONY: compute_energy_ewald.o compute_energy_ewald.mod
compute_energy_ewald.o compute_energy_ewald.mod: global_variables.mod \
  input_output.mod compute_energy_ewald.f90
	
	gfortran -g -c compute_energy_ewald.f90

# .PHONY: input_output.o input_output.mod
input_output.o input_output.mod: global_variables.mod input_output.f90
	
	gfortran -g -c input_output.f90
	
# .PHONY: global_variables.o global_variables.mod
global_variables.o global_variables.mod: global_variables.f90
	
	gfortran -g -c global_variables.f90
	
clean:
	
	rm compute_energy_ewald.o global_variables.o main.o \
		 input_output.o initialize_update.o \
		 compute_energy_ewald.mod global_variables.mod \
		 input_output.mod initialize_update.mod \ 
		 main