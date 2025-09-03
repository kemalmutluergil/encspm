import numpy as np
import random
import numpy as np
import random
import numpy as np
import random

def generate_hs_matrix(
    N: int,
    num_diagonals: int,
    filename: str = "hs_matrix.mtx",
    density_low: float = 0.90,
    density_high: float = 1.00,
    seed: int | None = None,
):
    """
    Generate an NxN HS-diagonal matrix:
      - Choose `num_diagonals` HS diagonals (0..N-1), fill 90-100% of each.
      - Apply a random row permutation and a separate random column permutation.
      - Save in Matrix Market (.mtx), with comments describing diagonals and permutations.
    """
    if num_diagonals > N:
        raise ValueError("Number of diagonals cannot exceed N")
    if not (0.0 <= density_low <= density_high <= 1.0):
        raise ValueError("density_low/high must satisfy 0 <= low <= high <= 1")

    rng = np.random.default_rng(seed)

    # --- pick HS diagonals (0..N-1) ---
    diagonals = sorted(random.sample(range(N), num_diagonals))

    # --- build original COO lists without forming dense NxN ---
    rows, cols, vals = [], [], []
    for d in diagonals:
        i = np.arange(N)
        j = (i + d) % N
        dens = rng.uniform(density_low, density_high)
        mask = rng.random(N) < dens
        i, j = i[mask], j[mask]

        rows.extend(i.tolist())
        cols.extend(j.tolist())
        vals.extend([1] * len(i))

    rows = np.asarray(rows, dtype=np.int64)
    cols = np.asarray(cols, dtype=np.int64)
    vals = np.asarray(vals, dtype=np.int64)

    # --- separate random row and column permutations ---
    new_to_old_row = rng.permutation(N)
    old_to_new_row = np.empty(N, dtype=np.int64)
    old_to_new_row[new_to_old_row] = np.arange(N, dtype=np.int64)

    new_to_old_col = rng.permutation(N)
    old_to_new_col = np.empty(N, dtype=np.int64)
    old_to_new_col[new_to_old_col] = np.arange(N, dtype=np.int64)

    # --- apply permutations separately ---
    perm_rows = old_to_new_row[rows]
    perm_cols = old_to_new_col[cols]
    perm_vals = vals

    # sort by row,col for nicer Matrix Market file
    order = np.lexsort((perm_cols, perm_rows))
    perm_rows = perm_rows[order]
    perm_cols = perm_cols[order]
    perm_vals = perm_vals[order]

    # --- write Matrix Market with header & comments ---
    with open(filename, "w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"% HS diagonals (0..{N-1}) in original: {diagonals}\n")
        f.write(f"% Row permutation old->new: {old_to_new_row.tolist()}\n")
        f.write(f"% Col permutation old->new: {old_to_new_col.tolist()}\n")
        f.write(f"{N} {N} {len(perm_rows)}\n")
        for r, c, v in zip(perm_rows, perm_cols, perm_vals):
            f.write(f"{int(r)+1} {int(c)+1} {int(v)}\n")

    print(f"Matrix written to {filename}")



if __name__ == "__main__":
    sizes = [n for n in [256, 512, 1024, 2048, 4096] for _ in range(4)]
    num_diagonals = [10, 20, 40, 80, 20, 40, 80, 160, 40, 80, 160, 320, 80, 160, 320, 640, 160, 320, 640, 1280]

    assert len(sizes) == len(num_diagonals)

    for i in range(len(sizes)):
        filename = "./data/synth/permuted_" + str(sizes[i]) + "_" + str(num_diagonals[i]) + ".mtx"
        generate_hs_matrix(sizes[i], num_diagonals[i], filename)

# if __name__ == "__main__":
#     size = int(input("Enter matrix size: "))
#     num_diagonals = int(input("Enter number of diagonals: "))

#     filename = f"./data/synth/permuted_{size}_{num_diagonals}.mtx"
#     generate_hs_matrix(size, num_diagonals, filename)
