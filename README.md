This is an implementation of the algorithm presented in https://arxiv.org/abs/2410.22714 which computes the Selmer group of elliptic curves of the form $y^2 = x^3 + bx$ over $\mathbb{Q}(i)$ where $b\in\mathbb{Z}[i]$.  In particular, we define the class `GaussInt(x,y)` which represents a Gaussian integer $x + yi$.  Among others, this class has the method `.factor()` which provides the unique factorization of a Gaussian integer into primary primes, along with the exponent of $i$ and $1+i$ which appear in the unique factorization.

We list the main functions which go into creating the algorithm.
- `Laplacian(b)` which creates the Laplacian matrix associated with a Gaussian integer `b` over $\mathbb{F}_2$ along with the keys in `b.factor()` and the primary primes of `b` which index the matrix.
- `Laplacian_prime(b)` which modifies the Laplacian matrix according to the main algorithm.  The indexing set of primes is also modified accordingly.
- `partial_Selmer_group(b, s_d, t_d)` will output elements $d$ in the Selmer group with specified values of $s_d$ and $t_d$.  In particular, we take each vector in the solution set of $L'x = y_{(s_d,t_d)}$, translate the vector as an element of $\mathbb{Z}[i]$, and then determine if they satisfy any of the three congruence condition to demonstrate local solubility at $1+i$.
- `Selmer_group(b)` computes the entire Selmer group of $E_b/\mathbb{Q}(i)$ by concatenating the four calls of `partial_Selmer_group()` with $s_d\in\\{ 0,1 \\}$ and $t_d\in\\{0,1\\}$.

The Jupyter notebook `testing.ipynb` consists of an empirical verification to Corollary 6.1 and Theorem 6.8 as well as a verification of Example 6.1 to test the implementation of our algorithm.  

For convenience, one can readily use the algorithm at https://selmergraph.streamlit.app/ for any choice of $b$ or a congruence class of $b$.
