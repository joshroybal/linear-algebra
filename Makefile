FC = gfortran
FCFLAGS = -pedantic-errors -std=f95 -Wall -Wextra -Werror -fcheck=all -fbacktrace -O2
FLFLAGS = -static -s
OBJ = driver.o linearalgebra.o
BIN = driver

$(BIN): $(OBJ)
	$(FC) -o $@ $^ $(FLFLAGS)

%.o: %.f90
	$(FC) -c $< $(FCFLAGS)

linearalgebra.mod := linearalgebra.o
driver.o: $(linearalgebra.mod)

.PHONY: clean
clean:
	$(RM) driver *.o *.mod *~
