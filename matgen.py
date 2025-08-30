import numpy as np
import random
from scipy.sparse import coo_matrix

def generate_hs_matrix(N, num_diagonals, filename="matrix.mtx"):
    if num_diagonals > N:
        raise ValueError("Number of diagonals cannot exceed N")

    # pick HS diagonals randomly
    diagonals = random.sample(range(N), num_diagonals)

    # create NxN zero matrix
    mat = np.zeros((N, N), dtype=int)

    # fill selected diagonals
    for d in diagonals:
        for i in range(N):
            j = (i + d) % N
            mat[i, j] = 1  # you could put random values instead

    # create random permutation
    perm = list(range(N))
    random.shuffle(perm)

    # apply permutation (P * M * P^T)
    permuted = mat[np.ix_(perm, perm)]

    # convert to sparse coordinate format
    sparse_mat = coo_matrix(permuted)

    # write manually to Matrix Market file
    with open(filename, "w") as f:
        # custom comments
        f.write(f"% Generated HS diagonal matrix\n")
        f.write(f"% N = {N}\n")
        f.write(f"% Filled diagonals: {diagonals}\n")
        f.write(f"% Applied permutation: {perm}\n")
        # standard header
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{N} {N} {sparse_mat.nnz}\n")
        # write entries (Matrix Market is 1-based)
        for i, j, v in zip(sparse_mat.row, sparse_mat.col, sparse_mat.data):
            f.write(f"{i+1} {j+1} {v}\n")

    print(f"Matrix written to {filename}")
    print(f"Filled diagonals: {diagonals}")
    print(f"Permutation: {perm}")


if __name__ == "__main__":
    sizes = [256, 512, 1024, 1024, 2048, 2048, 2048, 2048, 4096, 4096]
    num_diagonals = [10, 20, 40, 80, 80, 160, 320, 400, 400, 800]
    
    for i in range(len(sizes)):
        filename = "./data/synth/permuted_" + str(sizes[i]) + "_" + str(num_diagonals[i]) + ".mtx"
        generate_hs_matrix(sizes[i], num_diagonals[i], filename)
