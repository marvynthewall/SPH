CC = gcc
CFLAGS = -Wall -Wextra -O3 -I./include

SRC = src/grid_system.c src/riemann.c src/hydro.c src/grid_init.c src/grid_io.c
OBJ = $(SRC:.c=.o)

all: bin/grid_sod bin/grid_kh

bin/grid_sod: examples/grid_sod.c $(OBJ)
	@mkdir -p bin
	$(CC) $(CFLAGS) $^ -o $@ -lm

bin/grid_kh: examples/grid_kh.c $(OBJ)
	@mkdir -p bin
	$(CC) $(CFLAGS) $^ -o $@ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)
	rm -rf bin
