# Matrix Multiply (MPI Assignment)

Your task is to edit **`matmul_mpi.c`** so the matrix multiplication runs correctly on multiple MPI processes. Focus only on the region delimited by:

```
/* =========================== AREA OF INTEREST ============================ */

...

/* ========================= END AREA OF INTEREST ========================== */
```

MPI is already initialized for you, and a working sequential implementation is available in **`matmul.c`** for reference.

## Building

Use the provided `Makefile`:

- `make matmul`
- `make matmul_mpi`
- `make all`

## Testing

Run the full test suite with:

```
make test
```

This invokes `test_matmul.sh`, which executes both the sequential and MPI versions for several matrix sizes and verifies that they produce the expected output (`"Hey, it worked!"`).