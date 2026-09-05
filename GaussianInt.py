import sympy as sp
import math

def sum_two_squares(p):
    r = sp.sqrt_mod(-1, p, all_roots=True)[0]

    a, b = p, r

    while b * b > p:
        a, b = b, a % b

    c = b
    d = int(sp.sqrt(p - c * c))

    return c, d    

class GaussInt:
    def __init__(self, a, b=0):
        self.a = int(a)
        self.b = int(b)
        
    def __add__(self, other):
        return GaussInt(self.a + other.a,
                                self.b + other.b)
        
    def __sub__(self, other):
        return GaussInt(self.a - other.a,
                        self.b - other.b)

    def __mul__(self, other):
        return GaussInt(self.a*other.a - self.b*other.b,
                                self.a*other.b + self.b*other.a)
    
    def __neg__(self):
        return GaussInt(-self.a, -self.b)
    
    def __eq__(self, other):
        if isinstance(other, GaussInt):
            return self.a == other.a and self.b == other.b
        elif isinstance(other, int):
            return self.a == other and self.b == 0
        return NotImplemented  
    
          
    def __hash__(self):
        return hash((self.a, self.b))
    
    def __str__(self):
        if self.b == 0:
            return str(self.a)
        elif self.a == 0:
            return f"{self.b}i"
        elif self.b > 0:
            return f"{self.a} + {self.b}i"
        else:
            return f"{self.a} - {abs(self.b)}i"
        

    def __repr__(self):
        return f"GaussInt({self.a}, {self.b})"
    
    def __pow__(self, n):
        if not isinstance(n, int):
            return NotImplemented

        result = GaussInt(1)
        base = self

        if n < 0:
            raise ValueError('Exponent must be nonnegative')

        while n > 0:
            if n % 2 == 1:
                result = result * base
            base = base * base
            n //= 2

        return result

    def __mod__(self, other):
        denominator = other.norm()

        real_num = self.a * other.a + self.b * other.b
        imag_num = self.b * other.a - self.a * other.b

        # Round numerator/denominator to nearest integer
        def nearest_integer(num, den):
            if num >= 0:
                return (2 * num + den) // (2 * den)
            else:
                return -((2 * (-num) + den) // (2 * den))

        q_real = nearest_integer(real_num, denominator)
        q_imag = nearest_integer(imag_num, denominator)

        q = GaussInt(q_real, q_imag)

        return self - q * other
    
    def norm(self):
        return self.a**2 + self.b**2
    
    def conj(self):
        return GaussInt(self.a, -self.b)
    
    def divides(self, other):
        denominator = self.norm()

        real_num = self.a * other.a + self.b * other.b
        imag_num = self.a * other.b - self.b * other.a

        return (real_num % denominator == 0 and
                imag_num % denominator == 0)
        
    def quotient(self, other):
        denominator = other.norm()

        real_num = self.a * other.a + self.b * other.b
        imag_num = self.b * other.a - self.a * other.b

        if real_num % denominator != 0 or imag_num % denominator != 0:
            raise ValueError(f"{other} does not divide {self}")

        return GaussInt(real_num // denominator,
                        imag_num // denominator) 
    
    #def cong(self,other,mod):
    #    return mod.divides(self - other)
        
    def is_primary(self):
        t3 = GaussInt(1,1)**3
        z = self - GaussInt(1,0)
        return t3.divides(z)
    
    def primary_factor(self):
        for k in range(4):
            unit = GaussInt(0,1)**k
            z = unit * self
            if z.is_primary():
                return [(-k) % 4, z]
            
            
    def factor(self):
        z = self
        if z == GaussInt(0):
            raise ValueError('0 does not have a unique factorization')
        
        # initializing factors
        factors = {GaussInt(0,1):0, GaussInt(1,1):0} 
        # handling the 1+i factor
        while GaussInt(1, 1).divides(z):
            z = z.quotient(GaussInt(1, 1))
            factors[GaussInt(1, 1)] += 1
        # factoring the norm
        norm_factors = sp.factorint(z.norm())
        for p in norm_factors:
            if p % 4 == 3:
                prime = GaussInt(p)
                k, primary = prime.primary_factor()
                while prime.divides(z):
                    z = z.quotient(prime)
                    factors[GaussInt(0,1)] += k
                    factors[primary] = factors.get(primary, 0) + 1
            
            elif p % 4 == 1:
                a , b = sum_two_squares(p)
                pi = GaussInt(a,b)
                pi_bar = pi.conj()
 
                k1, pi_primary = pi.primary_factor()
                k2, pi_bar_primary = pi_bar.primary_factor()
                           
                while pi.divides(z):
                    z = z.quotient(pi)
                    factors[GaussInt(0, 1)] += k1
                    factors[pi_primary] = factors.get(pi_primary, 0) + 1

                while pi_bar.divides(z):
                    z = z.quotient(pi_bar)
                    factors[GaussInt(0, 1)] += k2
                    factors[pi_bar_primary] = factors.get(pi_bar_primary, 0) + 1                
        for k in range(4):
            if z == GaussInt(0, 1)**k:
                factors[GaussInt(0, 1)] += k
                break
        
        else:
                raise ValueError(f"Factorization failed; remaining factor: {z}")
        factors[GaussInt(0,1)] = factors[GaussInt(0,1)] % 4
        return factors                    
                    
    def s(self):
        # this is the power of i that appears in the primary factorization    
        f = self.factor()
        return f[GaussInt(0,1)]
    
    def t(self):
        # this is the power of 1+i that appears in the primary factorization
        count = 0
        z = self
        while GaussInt(1, 1).divides(z):
            count += 1
            z = z.quotient(GaussInt(1,1))
        return count
    
    def m(self):
        return log_i(quartic_res(GaussInt(1,1), self))
    
    def n(self):
        return log_i(quartic_res(GaussInt(0,1), self))
    
    def is_prime(self):
        a, b = abs(self.a), abs(self.b)

        if a == 0 and b == 0:
            return False

        if a == 0 or b == 0:
            n = max(a, b)
            return sp.isprime(n) and n % 4 == 3

        return sp.isprime(a*a + b*b)

def log_i(z):
    if z not in [GaussInt(0,1)**k for k in range(4)]:
        raise ValueError("input is not a power of i")
    for k in range(4):
        if z == GaussInt(0,1)**k:
            return int(k)
    

def quartic_res(z, p):
    if p == GaussInt(0):
        raise ValueError("Modulus cannot be zero")

    if not p.is_primary():
        raise ValueError("Modulus must be primary") 
    
    if p.divides(z):
        return GaussInt(0)
    
    exp = ((p.norm() - 1) // 4) 
    x = z**exp
    for k in range(4):
        if p.divides(x - GaussInt(0,1)**k):
            return GaussInt(0,1)**k
    else:
            raise ValueError("There was an error in calculated the residue symbol")
            



def gaussian_prime_congruent(norm_size, a, m):
    """
    Find a Gaussian prime p with
        N(p) < norm_size
        p == a (mod m)

    Parameters
    ----------
    norm_size : int
        Upper bound for the norm.
    a, m : GaussInt
        Congruence a modulo m.

    Returns
    -------
    GaussInt
        A Gaussian prime satisfying the conditions.

    Raises
    ------
    ValueError
        If no such prime is found below norm_size.
    """

    # We need p = a + q*m
    max_q_norm = norm_size // m.norm()

    max_q = math.isqrt(max_q_norm) + 1

    for x in range(-max_q, max_q + 1):
        for y in range(-max_q, max_q + 1):

            q = GaussInt(x, y)
            p = a + q * m

            if p.norm() >= norm_size:
                continue

            if p.is_prime():
                return p

    raise ValueError(
        f"No Gaussian prime congruent to {a} modulo {m} "
        f"with norm < {norm_size}"
    )

def gaussian_primes(norm_size, a, m):
    primes = []

    bound = math.isqrt(norm_size - 1)

    # Non-real/non-imaginary Gaussian primes
    for x in range(1, bound + 1):
        for y in range(1, bound + 1):

            n = x*x + y*y

            if n >= norm_size:
                continue

            # Gaussian prime iff norm is a rational prime
            if sp.isprime(n):
                p = GaussInt(x, y)

                if m.divides(p - a):
                    primes.append(p)

    # Rational Gaussian primes p == 3 mod 4
    for p in sp.primerange(3, norm_size):
        if p % 4 == 3:

            # Positive real representative
            z = GaussInt(p, 0)

            if m.divides(z - a):
                primes.append(z)

    # 1+i
    z = GaussInt(1, 1)

    if z.norm() < norm_size and m.divides(z - a):
        primes.append(z)

    return primes    
    
t = GaussInt(1,1)
t7 = GaussInt(1,1)**7
    
t7_residues = {
    GaussInt(a, b)
    for a in range(8)
    for b in range(16)
}

t7_squares = {z**2 % t7 for z in t7_residues}