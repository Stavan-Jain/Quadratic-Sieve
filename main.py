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
        self.old_rows = None
        self.old_matrix = None
        self.reduced_rows = None
        self.old_bsmooth = None
        self.lincombs = dict()
        self.tonelli_relations = dict({})


     
    def mod_exponentiation(self, y, p):
        z = 1 
        x = self.n % p
        
        while (y > 0):
            if (y & 1):
                z = (z * x) % p

            y = y >> 1 
            x = (x * x) % p
        return z

    def eulers_criterion(self, primes):
        for p in primes: 
            if((self.n**(int((p-1)/2)))%p == p-1):
            # if(mod_exponentiation((int((p-1)/2)),p)):
                primes.remove(p)
        return primes

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

    def tonelli_repeated(self, square, prime, prime_power):
        old_mod = prime ** (prime_power - 1)
        new_relations = set()
        for x_k in self.tonelli_relations[prime][prime_power - 1]:
            y_k = self.find_inverse(2 * x_k, old_mod)
            x_next = (x_k - ((x_k ** 2) - square) * y_k) % (old_mod * prime)
            new_relations.add(x_next)
        new_relations = list(new_relations)
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
        if prime_power > 1:
            return self.tonelli_repeated(square, prime, prime_power)
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
                self.tonelli_relations.update({prime: [math.log(prime), self.tonelli_shanks(self.n, prime, 1)]})
            power = 0
            while power >= 0:
                if len(self.tonelli_relations[prime]) < (power + 2):
                    self.tonelli_relations[prime].append(self.tonelli_shanks(self.n, prime, power + 1))
                if value % (prime ** (power + 1)) in self.tonelli_relations[prime][power + 1]:
                    # print(f'found that %d ^ %d is a factor' % (prime, power + 1))
                    power += 1
                    current = current - self.tonelli_relations[prime][0]
                    factors[count] += 1
                else:
                    power = -1
                if current < 0.001:
                    power = -1
        if current < 0.001:
            return 1, factors
        return -1, factors

    def find_bsmooth(self, B, tonelli=True):
        primes = self.gen_primes(B)
        print("Done generating primes")
        sq = int(math.sqrt(self.n))
        i = 1
        while len(self.matrix) <= len(primes):
            temp = sq + i
            current = ((temp)**2) % self.n
            
            if tonelli:
                factored, factors = self.tonelli_wrapper(current, temp)
            else:
                factored, factors = self.factor_with_base(primes, current)
            
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
    
    #row reduces M (mod 2) while keeping track of the linear combinations in each row in "linear_combinations"
    def row_reduce(self, M, linear_combs):
        linear_combinations = linear_combs.copy()
        A = M.copy()
        m, n = np.shape(A)
        h, k = 0, 0
        while h < m and k < n:
            i_max = np.argmax(np.abs(A[h:, k])) + h      #getting the index of pivot row in column k
            if not A[i_max][k]:                #if there is no pivot in this column, move to the next one
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
        return A, linear_combinations

    #in the Reduced Echelon Form A of M, find the 0 rows (i.e. True not in row A[i]) and for each such row, 
    #store the linear combination of rows that yield 0 in `linear_dependencies`.
    #given a dictionary of linear combinations for each row, find the linear dependencies
    def find_lin_dependencies(self, M, linear_combs):
        numrows, _ = np.shape(M)
        linear_dependencies = []
        for i in range(numrows):
            if True not in M[i]:
                linear_dependencies.append(filter_array(linear_combs[i][1]))
        return linear_dependencies

    # given new exponent vectors `newrows` and new bsmooth numbers `newbsmooth`, appropriately concatenate them with the old arrays
    # M: exponent vectors, A : exponent vectors reduced mod 2 (with old rows row-reduced), v: bsmooth numbers
    def initialize_objects(self, newrows, newbsmooth):
        nr = newrows % 2 != 0
        linear_combinations = self.lincombs
        if self.reduced_rows is not None: A = np.concatenate((self.reduced_rows, nr)) #A contains the exponent vectors mod 2 to be row reduced
        else: A = nr
        if self.old_matrix is not None : M = np.concatenate((self.old_matrix, newrows)) #M is the complete (unreduced) 'exponent vector' matrix
        else: M = newrows
        if self.old_bsmooth is not None : v = np.concatenate((self.old_bsmooth, newbsmooth)) #v is a vector that contains all the bsmooth numbers
        else: v = newbsmooth
        m , n = np.shape(A)
        l = len(linear_combinations)
        for i in range(len(newbsmooth)):
            linear_combinations[l + i] = [newbsmooth[i], [newbsmooth[i]]]
        return A, M, v, linear_combinations
    
    #newrows: new exponent vectors (2d array). dimensions (m x n)
    #newbsmooth: new bsmooth numbers corresponding to new rows. dimensions (m)
    #returns two arrays C and B. C[i]**2 is congruent to B[i]**2 mod n for all i. 
    def find_congruent_squares(self, newrows, newbsmooth, factorbase): 
        A, M, v, linear_combinations = self.initialize_objects(newrows, newbsmooth)
        A, linear_combinations = self.row_reduce(A, linear_combinations)
        linear_dependencies = self.find_lin_dependencies(A, linear_combinations)

        #tracks information for next call to `lin_dep_mod2`
        self.old_matrix = A 
        self.old_bsmooth = v
        self.lincombs = linear_combinations
        self.old_rows = M

        B = []
        #print(linear_dependencies)
        #for each linear dependency, compute the corresponding product of prime powers (mod n) and store in an array B
        for dependency in linear_dependencies:
            indices = []
            for k in dependency:
                indices.append(np.where(v == k)[0])  #finding indices corresponding to B-smooth numbers
            prime_exp = sum([M[i]for i in indices])[0]
            #print(prime_exp)
            prime_exp = prime_exp // 2

            b = 1
            for i in range(len(factorbase)):
                b = (b* pow(factorbase[i].item(),prime_exp[i].item(),self.n)) % self.n
            #b = np.prod(factorbase**prime_exp) #computing product of prime powers
            # .item() converts np.int64 to int
            B.append(b % self.n)

        #for each linear dependency, compute the product of the B-smooth numbers (mod n) and store in an array C
        C = []
        for dep in linear_dependencies:
            prod = 1
            for e in dep:
                prod = (prod* e.item()) % self.n
            C.append(prod)
        #C = [np.prod(arr).item() % self.n for arr in linear_dependencies]  

        return C, B 

    def basic_principle(self, a, b):
        #print(type(a), type(b))
        if((a-b)%self.n==0 or (a+b)%self.n==0):
            return False
        else: 
            return math.gcd(abs(a-b), self.n)
    
    #driver code
    def find_prime_factor(self, tonelli=True):
        B = self.get_B()
        self.find_bsmooth(B, tonelli)
        A, C = self.find_congruent_squares(self.matrix, self.bsmooth, self.factor_base)
        ret = []
        for i in range(len(A)):
            j = self.basic_principle(A[i], C[i])
            ret.append(j)
        return ret

#Sieve = QuadraticSieve(3837523)       
#Sieve = QuadraticSieve(77340247)
#Sieve = QuadraticSieve(100109*100271)
#Sieve = QuadraticSieve(100109* 386429)
#Sieve = QuadraticSieve(100271* 5009317 )
#Sieve = QuadraticSieve(10023234*12345679)
Sieve = QuadraticSieve(310248241 * 383838383)
#Sieve = QuadraticSieve(16921456439215439701)
#Sieve = QuadraticSieve(384869498225059)
print(Sieve.find_prime_factor())
