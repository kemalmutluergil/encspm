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

