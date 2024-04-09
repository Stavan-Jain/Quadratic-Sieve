import numpy as np
import math

class QuadraticSieve:
    
    def __init__(self, n):
        self.n = n
        self.matrix = np.array([])
    
    def gen_primes(self, limit):
        primes = [2]
        i=3
        while i < limit:
            is_prime = True
            for p in primes:
                if i % p == 0:
                    is_prime = False
            if is_prime:
                primes.append(i)
            i += 1
        return primes

    def factor_with_base(self, base, target):
        factors = [target] + [0] * len(base)
        for count, prime in enumerate(base):
            while (factors[0] % int(prime)) == 0:
                factors[0] = int(factors[0] / prime)
                factors[count+1] += 1
        if factors[0] != 1:
            factors[0] = -1
        return factors

    def find_bsmooth(self, B):
        primes = self.gen_primes(B)
        sq = int(math.sqrt(self.n))
        i = 1
        while len(self.matrix) <= len(primes):
            temp = sq + i
            current = ((temp)**2) % self.n
            factors = self.factor_with_base(primes, current)
            if factors[0] == 1:
                factors[0] = temp
                if len(self.matrix) == 0:
                    self.matrix = np.array([factors])
                else:
                    self.matrix = np.append(self.matrix, [factors], axis=0)
            i += 1
        return 
    
    #TODO: Prerana
    def get_B(self):
        return 
    
    #TODO: Stavan
    def find_nullspace(self):
        return 

    #TODO: Prerana
    def basic_principle(self):
        return False
    
x = QuadraticSieve(1649)
x.find_bsmooth(12)
print(x.matrix)