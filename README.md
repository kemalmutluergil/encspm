RCM and GOrder do roughly the same job when it comes to accomodating HS diagonals. 
Current dummy implementation takes 2 seconds for 256x256 matrix using N = 8192. Note that normally this would mean 32 computations at the same time.

# TODO:
- Fix the issues with square matrices of larger size
- Optimize `hs_mul` to use ciphertext slots efficiently
- Accomodate matrices of size larger than $N$
- Devise a way to tailor permutation specifically for HS diagonals

# Notes 

The RCM ordering does not necessarily help us reduce the number of occupied HS diagonals. Example:

$$
A = \begin{bmatrix}
10 & 0 & 2.5 & 0 & 0 \\
0 & 5 & 0 & 0 & 1.2 \\
3.3 & 0 & 7.1 & 0 & 0\\
0 & 0 & 0 & 4 & 0 \\
0 & 6.5 & 7.2 & 0 & 0 \\
\end{bmatrix}
$$
This matrix currently occupies 3 HS diagonals (0, 2, 3)

But the result of the RCM reordering of this matrix occupies 4 (0, 1, 3, 4):

$$
A' = \begin{bmatrix}
5 & 1.2 & 0 & 0 & 0 \\
6.5 & 0 & 7.2 & 0 & 0\\
3.3 & 0 & 7.1 & 3.3 & 0\\
0 & 0 & 2.5 & 10 & 0 \\
0 & 0 & 0 & 0 & 4 \\
\end{bmatrix}
$$

