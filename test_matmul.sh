#!/bin/bash

N_DIMS=(1000 1500 2000)
FAILED=0

run_and_check() {
    local label="$1"
    local cmd="$2"

    echo ">>> $label"
    if ! OUTPUT=$(eval "$cmd" 2>&1); then
        echo "$OUTPUT"
        echo "[ERROR] Command failed: $cmd"
        FAILED=1
        return
    fi

    echo "$OUTPUT"
    if ! grep -q "Hey, it worked!" <<< "$OUTPUT"; then
        echo "[ERROR] Validation string missing for: $label"
        FAILED=1
    fi
}

for size in "${N_DIMS[@]}"; do
    run_and_check "matmul_mpi (order $size, 4 ranks)" "mpirun -np 4 ./matmul_mpi.out $size"
done

if [[ $FAILED -ne 0 ]]; then
    echo "Some tests failed."
    exit 1
fi

echo "All tests passed."