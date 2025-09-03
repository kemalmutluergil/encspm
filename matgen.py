import numpy as np
import random
from scipy.sparse import coo_matrix
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
      - Apply a random row/col permutation (same for rows and cols).
      - Save in Matrix Market (.mtx), with comments describing diagonals and permutations.

    Comments include:
      - HS diagonals (original)
      - Permutation old->new
      - Permutation new->old

    Matrix Market banner is written first (per spec), then comments.
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
        # positions on HS diagonal d are (i, (i+d) % N) for i=0..N-1
        i = np.arange(N)
        j = (i + d) % N
        # 90–100% density per diagonal
        dens = rng.uniform(density_low, density_high)
        mask = rng.random(N) < dens
        i, j = i[mask], j[mask]

        rows.extend(i.tolist())
        cols.extend(j.tolist())
        vals.extend([1] * len(i))  # use 1s; change if you want random weights

    rows = np.asarray(rows, dtype=np.int64)
    cols = np.asarray(cols, dtype=np.int64)
    vals = np.asarray(vals, dtype=np.int64)

    # --- build permutation ---
    # new_to_old[k] = old index that becomes row/col k in the permuted matrix
    new_to_old = rng.permutation(N)
    # old_to_new[old] = new index of that old row/col
    old_to_new = np.empty(N, dtype=np.int64)
    old_to_new[new_to_old] = np.arange(N, dtype=np.int64)

    # --- apply permutation as old->new (easy to reason about) ---
    # For every old coordinate (r, c), new coordinate is (old_to_new[r], old_to_new[c])
    perm_rows = old_to_new[rows]
    perm_cols = old_to_new[cols]
    perm_vals = vals

    # (optional) sort entries by row,col for nicer files
    order = np.lexsort((perm_cols, perm_rows))
    perm_rows = perm_rows[order]
    perm_cols = perm_cols[order]
    perm_vals = perm_vals[order]

    # --- write Matrix Market with correct header order ---
    with open(filename, "w") as f:
        # Banner FIRST (per Matrix Market spec)
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        # Comments (each line starting with %)
        f.write(f"% HS diagonals (0..{N-1}) in original: {diagonals}\n")
        f.write(f"% Permutation old->new: {old_to_new.tolist()}\n")
        # Problem size & nnz
        f.write(f"{N} {N} {len(perm_rows)}\n")
        # 1-based entries
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
