all: matmul matmul_mpi

matmul:
	mpicc -Wall -pedantic matmul.c -o matmul.out

matmul_mpi:
	mpicc -Wall -pedantic matmul_mpi.c -o matmul_mpi.out

test: matmul matmul_mpi
	./test_matmul.sh

clean:
	rm -f matmul.out matmul_mpi.out