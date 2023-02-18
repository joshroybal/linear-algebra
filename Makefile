FC = gfortran
FCFLAGS = -pedantic-errors -std=f95 -Wall -Wextra -Werror -O2 -g
FLFLAGS =
SRC = prng.f90 linearalgebra.f90 driver.f90
OBJ = ${SRC:.f90=.o}
BIN = driver

$(BIN): $(OBJ)
	$(FC) -o $@ $^ $(FLFLAGS)

%.o: %.f90
	$(FC) -c $< $(FCFLAGS)

.PHONY: clean
clean:
	$(RM) driver *.o *.mod *~
