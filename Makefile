FC = gfortran
FCFLAGS = -pedantic-errors -std=f95 -Wall -Wextra -Werror -fcheck=all -fbacktrace -g -Og
FLFLAGS =

driver: driver.o matrix.o
	$(FC) -o $@ $^ $(FLFLAGS)

matrix.o: matrix.f90
	$(FC) -c $< $(FCFLAGS)

driver.o: driver.f90 matrix.o
	$(FC) -c $< $(FCFLAGS)

.PHONY: clean
clean:
	$(RM) driver *.o *.mod *~
