import numpy as np
from utils import *
import math 
import time

class QuadraticSieve:
    
    def __init__(self, n):
        self.n = n
        self.matrix = np.array([])
        self.bsmooth = np.array([])
        self.factor_base = np.array([])
        self.tonelli_relations = dict({})
    
    def gen_primes(self, limit):
        primes = [2]
        i=3
        while i <= limit:
            is_prime = True
            for p in primes:
                if i % p == 0:
                    is_prime = False
            if is_prime:
                primes.append(i)
            i += 1
        self.factor_base = np.array(primes)
        return primes

    def factor_with_base(self, base, target):
        temp = target
        factors = [0] * len(base)
        for count, prime in enumerate(base):
            while (temp % int(prime)) == 0:
                temp = int(temp / prime)
                factors[count] += 1
        if temp != 1:
            temp = -1
        return temp, factors

    def tonelli_2(self, square, prime, prime_power):
        if prime_power == 1:
            if square % prime == 1:
                return (1, 1)
            return (-1, -1)
        if prime_power == 2:
            if square % (prime ** prime_power) == 1:
                return (1, 3)
            return (-1, -1)
        if prime_power == 3:
            if square % 8 != 1:
                return (-1, -1)
            return (1, 3, 5, 7)
        new_relations = []
        mod = 2 << (prime_power - 1)
        for relation in self.tonelli_relations[prime][prime_power-1]:
            indicator = ((relation ** 2) - square) >> (prime_power - 1)
            if indicator % 2 == 0:
                new_relations.append(relation % mod)
            else:
                new_relations.append((relation + (2 << (prime_power - 3))) % mod)
        new_relations = new_relations + [mod - relation for relation in new_relations]
        new_relations = list(set(new_relations))
        new_relations.sort()
        return new_relations

    def tonelli_shanks(self, square, prime, prime_power=1):
        a = square
        b = 0
        power = 0
        k = 1
        m = int((prime - 1) / 2)
        if self.fast_powers(a, m, prime) != 1:
            return (-1, -1)
        if prime == 2:
            return self.tonelli_2(square, prime, prime_power)
        if square == 1:
            return (1, prime - 1)
        
        while True:
            legendre = self.fast_powers(a, m, prime)
            if legendre == 1:
                if m % 2 == 1:
                    break
                else:
                    m = m >> 1
                k = k + 1
            if legendre == prime - 1:
                if b == 0:
                    b = self.find_non_residue(prime)
                power = power + (2 ** (k - 1))
                a = (a * self.fast_powers(b, (2 ** (k-1)), prime)) % prime
        final_sqrt = self.fast_powers(a, int((m+1)/2), prime)
        if power > 0:
            b_inverse = self.find_inverse(b, prime)
            answer = final_sqrt * self.fast_powers(b_inverse, power // 2, prime) % prime
        else:
            answer = final_sqrt
        return (answer, (-answer) % prime)
                
    def fast_powers(self, base, exponent, mod):
        exponent = int(exponent)
        ret_val = 1
        while exponent > 0:
            if exponent % 2 == 1:
                ret_val = (ret_val * base) % mod
                exponent = exponent - 1
            base = (base * base) % mod
            exponent = exponent >> 1
        return ret_val

    def find_non_residue(self, mod):
        if mod == 2:
            return 1
        value = 2
        exponent = (mod - 1) / 2
        while True:
            current = self.fast_powers(value, exponent, mod)
            if current == mod - 1:
                return value
            value = value + 1

    def find_inverse(self, value, mod):
        # convention r = r_a * q_a + r_b
        coefficients = []
        r = mod
        r_a = value
        r_b = 1
        while r_b != 0:
            q_a = r // r_a
            r_b = r % r_a
            coefficients.append(q_a)
            r = r_a
            r_a = r_b
        top_row = [0, 1]
        bottom_row = [1, 0]
        for c in coefficients:
            top_row.append(c * top_row[-1] + top_row[-2])
            bottom_row.append(c * bottom_row[-1] + bottom_row[-2])
        if ((top_row[-2] * bottom_row[-1]) - (top_row[-1] * bottom_row[-2])) == -1:
            return (-1 * top_row[-2]) % mod
        return top_row[-2] % mod

    def tonelli_wrapper(self, candidate, value):
        current = math.log(candidate)
        factors = [0] * len(self.factor_base)
        for count, prime in enumerate(self.factor_base):
            if prime not in self.tonelli_relations:
                self.tonelli_relations.update({prime: [math.log(prime), self.tonelli_shanks(self.n % prime, prime, 1)]})
            power = 0
            while power >= 0 and prime == 2:
                if len(self.tonelli_relations[prime]) < (power + 2):
                    self.tonelli_relations[prime].append(self.tonelli_shanks(self.n % prime, prime, power + 1))
                if value % (prime ** (power + 1)) in self.tonelli_relations[prime][power + 1]:
                    print(f'found that %d ^ %d is a factor' % (prime, power + 1))
                    power += 1
                    current = current - self.tonelli_relations[prime][0]
                    factors[count] += 1
                else:
                    power = -1
                if current < 0.001:
                    power = -1
            if prime > 2: # will delete this clause when possible
                if value % prime in self.tonelli_relations[prime][1]:
                    current = current - self.tonelli_relations[prime][0]
                    factors[count] += 1
        if current < 0.001:
            return 1, factors
        return -1, factors

    def find_bsmooth(self, B):
        primes = self.gen_primes(B)
        sq = int(math.sqrt(self.n))
        i = 1
        while len(self.matrix) <= len(primes):
            temp = sq + i
            current = ((temp)**2) % self.n
            factored, factors = self.factor_with_base(primes, current)
            # factored, factors = self.tonelli_wrapper(current, i)
            if factored == 1:
                if len(self.matrix) == 0:
                    self.matrix = np.array([factors])
                    self.bsmooth = np.array(temp)
                else:
                    self.matrix = np.append(self.matrix, [factors], axis=0)
                    self.bsmooth = np.append(self.bsmooth, temp)
            i += 1
        return 
    
    def get_B(self):
        B = np.exp((1/2)*math.sqrt(math.log(self.n)*math.log(math.log(self.n))))
        return math.ceil(B)
    
    #returns two arrays A and B. A[i]**2 is congruent to B[i]**2 mod n for all i. 
    def find_linear_dependency_mod_2(self, M, v, factor_base):
        A = M % 2 != 0   # reducing the matrix mod 2 (True = 1, False = 0)
        m , n = np.shape(A)
        if m != len(v):
            raise Exception("length of vector != number of rows in matrix")
        if n != len(factor_base):
            raise Exception("length of factor base != number of columns in matrix")
        h, k = 0, 0                              #h : pivot row, k : pivot col
        linear_combinations = dict()    #dictionary to keep track of row operations
        for i in range(m):
            linear_combinations[i] = [v[i], [v[i]]]
        #this section computes the Row Echelon form of A while keeping track of row operations in `linear_combinations`
        while h < m and k < n:
            i_max = np.argmax(np.abs(A[h:, k])) + h      #getting the index of pivot column
            if not A[i_max][k]:                #if there is not pivot in this column, move to the next one
                k += 1
            else:
                r1, r2 = A[h].copy(), A[i_max].copy()
                A[h], A[i_max] = r2, r1         #swapping the current row with the row with pivot
                lc1, lc2 = linear_combinations[h].copy(), linear_combinations[i_max].copy()
                linear_combinations[h], linear_combinations[i_max] = lc2, lc1 #keeping track of index for each number in the vector v
                for i in range(h+1, m):
                    if A[i][k]:
                        linear_combinations[i][1] = np.append(linear_combinations[i][1], linear_combinations[h][1])
                        A[i] = np.logical_xor(A[i], A[h]) # updating rows so that there are 0s under the pivot entry 
                h += 1
                k += 1

        #in the Reduced Echelon Form A of M, find the 0 rows (i.e. True not in row A[i]) and for each such row, 
        #store the linear combination of rows that yield 0 in `linear_dependencies`
        linear_dependencies = []
        for i in range(m):
            if True not in A[i]:
                linear_dependencies.append(filter_array(linear_combinations[i][1]))

        B = []
        #for each linear dependency, compute the corresponding product of prime powers (mod n) and store in an array B
        for dependency in linear_dependencies:
            indices = []
            for k in dependency:
                indices.append(np.where(v == k)[0])  #finding indices corresponding to B-smooth numbers
            prime_exp = sum([M[i]for i in indices])[0]
            prime_exp = prime_exp // 2
            b = np.prod(factor_base**prime_exp) #computing product of prime powers
            B.append(b % self.n)

        #for each linear dependency, compute the product of the B-smooth numbers (mod n) and store in an array A
        A = [np.prod(arr) % self.n for arr in linear_dependencies] 
        
        return A, B

    def basic_principle(self, a, b):
        if((a-b)%self.n==0 or (a+b)%self.n==0):
            return False
        else: 
            return math.gcd(abs(a-b), self.n)
    
    #driver code
    def find_prime_factor(self):
        B = self.get_B()
        self.find_bsmooth(B)
        A, B = Sieve.find_linear_dependency_mod_2(Sieve.matrix, Sieve.bsmooth, Sieve.factor_base)
        ret = []
        for i in range(len(A)):
            ret.append(Sieve.basic_principle(A[i], B[i]))
        return ret

        
Sieve = QuadraticSieve(1073)
# print(Sieve.find_prime_factor())
Sieve.gen_primes(19)

current = time.time()
sol_A = Sieve.tonelli_wrapper(152, 35)
end = time.time()
print("Tonelli took " + str(end - current) + " seconds")
print(sol_A)
current = time.time()
sol_B = Sieve.factor_with_base(Sieve.factor_base, 152)
end = time.time()
print("Brute force took " + str(end - current) + " seconds")
print(sol_B)