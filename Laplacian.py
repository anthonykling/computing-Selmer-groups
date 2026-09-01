from GaussianInt import *
import sympy as sp
import numpy as np

def degv(b, v, e=0):
    factors = b.factor()
    if v not in factors.keys():
        raise ValueError('p is not a primary prime factor of b')
    
    ve = factors[v]
    if e != 0 and (ve % 2 == e % 2):
        raise ValueError('The parity of the exponent of v and of e should be different according to the main algorithm.')
    
    oneV = [p for p in factors.keys() if factors[p] == 1]
    twoV = [p for p in factors.keys() if factors[p] == 2]
    threeV = [p for p in factors.keys() if factors[p] == 3]
    
    sum = 0
    if e == 1:
        for p in oneV:
             sum = sum + log_i(quartic_res(p,v))
        return sum
    if e == 2:
        for p in twoV:
            sum = sum + log_i(quartic_res(p,v))
        return sum
    if e == 3:
        for p in threeV:
            sum = sum + log_i(quartic_res(p,v))
        return sum
    if e == 0:
        for p in factors.keys():
            if p == v:
                continue
            sum = sum + log_i(quartic_res(p,v))
        return sum
    
def adj_matrix(b):
    '''
    Computes the adjacency matrix associated to b over F2
    '''
    factors = b.factor()
    primes = [p for p in factors if p != GaussInt(1,1) and p != GaussInt(0,1)]
    #oddV = [p for p in factors.keys() if factors[p] % 2 == 1]
    #evenV = [p for p in factors.keys() if factors[p] == 2]    
    #N = len(oddV) 
    #M = len(evenV)
    N = len(primes)    
    primes.sort(key = lambda p: factors[p] % 2 == 0)
    
    matrix = np.empty((N,N), dtype=int)
    for i, p_i in enumerate(primes):
        for j, p_j in enumerate(primes):
            if i == j:
                matrix[i,j] = 0
            else:
                matrix[i,j] = log_i(quartic_res(p_j, p_i)) % 2
                
    return factors, primes, matrix

def Laplacian(b):
    '''Computes the Laplacian associated to b over F2'''
    factors, primes, matrix = adj_matrix(b)
    for i, p_i in enumerate(primes):
        matrix[i,i] = degv(b, p_i, e=0) % 2
    return factors, primes, matrix

def Laplacian_prime(b):
    factors, primes, matrix = Laplacian(b)
    
    # Step 2
    remove = [q for q in primes 
              if factors[q] ==2 and (degv(b, q, e=1) + degv(b, q, e=3)) % 2 != (q.m()*b.t() + q.n()*b.s()) % 2
              ]
    indices = [primes.index(p) for p in remove]
    matrix = np.delete(matrix, indices, axis=0)
    matrix = np.delete(matrix, indices, axis=1)
    primes = [p for p in primes if p not in remove]
    
    for i, p_i in enumerate(primes):
        p = p_i
        # Step 3
        if (factors[p] % 2 == 1) and (degv(b, p, e=2) % 2 != (p.m()*b.t() + p.n()*b.s()) % 2):
            matrix[i,i] += 1
            matrix[i,i] = matrix[i,i] % 2
        # Step 4    
        if (factors[p] == 2) and ((degv(b,p,e=1) + 3*degv(b,p,e=3)) % 4 != (p.m()*b.t() + p.n()*b.s() +2*p.n()) % 4):
            matrix[i,i] += 1
            matrix[i,i] = matrix[i,i] % 2
    
    
    return factors, primes, matrix

def matrix_solver_F2(A, y):
    A = A.copy() % 2
    y = np.array(y, dtype=int) % 2

    m, n = A.shape
    M = np.column_stack((A, y))

    pivot_cols = []
    row = 0

    # Row reduction over F_2
    for col in range(n):
        pivot = np.where(M[row:, col] == 1)[0]

        if len(pivot) == 0:
            continue

        pivot = pivot[0] + row
        M[[row, pivot]] = M[[pivot, row]]

        for r in range(m):
            if r != row and M[r, col]:
                M[r] ^= M[row]

        pivot_cols.append(col)
        row += 1

        if row == m:
            break

    # Check for inconsistency
    for r in range(row, m):
        if np.all(M[r, :-1] == 0) and M[r, -1] == 1:
            return []

    # Free variables
    free_cols = [j for j in range(n) if j not in pivot_cols]

    solutions = []

    # Enumerate all choices of free variables
    for values in range(2 ** len(free_cols)):

        x = np.zeros(n, dtype=int)

        for j, col in enumerate(free_cols):
            x[col] = (values >> j) & 1

        # Back substitution
        for r, col in reversed(list(enumerate(pivot_cols))):
            x[col] = (
                M[r, -1] ^
                (np.dot(M[r, :-1], x) % 2)
            )

        solutions.append(x)

    return solutions    