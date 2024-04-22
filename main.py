import numpy as np
from utils import *
import math 
import time

class QuadraticSieve:
    
    def __init__(self, n):
        self.n = int(n)
        self.matrix = np.array([])
        self.bsmooth = np.array([])
        self.factor_base = []
        self.old_rows = None
        self.old_matrix = None
        self.old_bsmooth = None
        self.lincombs = dict()
        self.tonelli_relations = dict({})
        self.tonelli = True
        self.i = 1
        self.tonelli_time = 0
        self.modular_plus = None
        self.modular_minus = None

    # use's euler's criterion to determine for which B-smooth primes is n a quadratic residue
    # if n is not a quadratic residue then the corresponding prime is removed from the factor base 
    def eulers_criterion(self, primes):
        new_primes = primes.copy()
        for p in primes: 
            if  self.fast_powers(self.n, int((p-1)/2), p) != 1:
                new_primes.remove(p)
        return new_primes

    # brute force generate primes under the B-smooth limit
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
        return primes

    # this is the brute force function to determine whether the target is b-smooth
    def factor_with_base(self, base, target):
        temp = target
        factors = [0] * len(base)
        for count, prime in enumerate(base):
            if prime == -1:
                continue
            while (temp % int(prime)) == 0:
                temp = int(temp / prime)
                factors[count] += 1
        if temp != 1:
            temp = -1
        return temp, factors

    # this function deals with the special case of p=2 for tonelli-shanks
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

    # this function uses the base tonelli-shanks relations to generate relations for powers of primes
    def tonelli_repeated(self, square, prime, prime_power):
        old_mod = prime ** (prime_power - 1)
        new_relations = set()
        for relation in self.tonelli_relations[prime][prime_power - 1]:
            x_k = int(relation)
            y_k = int(self.find_inverse(2 * x_k, old_mod))
            x_next = (x_k - ((x_k ** 2) - square) * y_k) % (old_mod * prime)
            new_relations.add(x_next)
        new_relations = list(new_relations)
        new_relations.sort()
        return new_relations

    # this function sends the case to tonelli_2 or tonelli_repeated if necessary
    # else, the function does the tonelli-shanks algorithm
    def tonelli_shanks(self, square, prime, prime_power=1):
        if prime == 2:
            return self.tonelli_2(square, prime, prime_power)
        if prime_power > 1:
            return self.tonelli_repeated(square, prime, prime_power)
        if square == 1:
            return (1, prime - 1)
        
        a = square
        b = 0
        power = 0
        k = 1
        m = int((prime - 1) / 2)

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
                a = (a * int(self.fast_powers(b, (2 ** (k-1)), prime))) % prime
        final_sqrt = self.fast_powers(a, int((m+1)/2), prime)
        if power > 0:
            b_inverse = self.find_inverse(b, prime)
            answer = final_sqrt * self.fast_powers(b_inverse, power // 2, prime) % prime
        else:
            answer = final_sqrt
        return (answer, (-answer) % prime)

    # this function computes powers modulo mod (quickly)            
    def fast_powers(self, base, exponent, mod):
        exponent = int(exponent)
        ret_val = 1
        while exponent > 0:
            if (exponent & 1):
                ret_val = (ret_val * base) % mod
                exponent = exponent - 1
            base = (base * base) % mod
            exponent = exponent >> 1
        return ret_val

    # this function finds a non-residue for use in the tonelli-shanks algorithm
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

    # this function uses the extended gcd algorithm to find modular inverses
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

    # this wrapper uses the generated tonelli-relations and generates more when needed
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
                    #print(f'found that %d ^ %d is a factor' % (prime, power + 1))
                    #print(f'%s %s^%s %s' % (value, prime, (power + 1), self.tonelli_relations[prime][power + 1]))
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

    # note tonelli only used for (sqrt(n)+k)^2 < 2n to keep consistent factor base
    # this function finds a specified number of b-smooth numbers
    def find_bsmooth(self, num_to_gen, tonelli=True, i=1):
        sq = int(math.sqrt(self.n))
        self.i = i
        current = 1
        self.matrix = np.array([])
        while len(self.matrix) <= num_to_gen:
            temp = sq + self.i
            current = ((temp)**2) % self.n
            if tonelli and (temp * temp) < 2 * self.n:
                factored, factors = self.tonelli_wrapper(current, temp)
            else:
                if self.tonelli:
                    print(f'Switching to Brute Force Solving %s' % temp)
                    self.tonelli = False
                factored, factors = self.factor_with_base(self.factor_base, current)
            if factored == 1:
                print(f'%d %d %d %s' % (temp, current, factored, factors))
                if len(self.matrix) == 0:
                    self.matrix = np.array([factors])
                    self.bsmooth = np.array(temp)
                else:
                    self.matrix = np.append(self.matrix, [factors], axis=0)
                    self.bsmooth = np.append(self.bsmooth, temp)
            self.i += 1
        return 
    
    def better_find_bsmooth(self, num_to_gen, tonelli=True, i=1):
        if tonelli == False:
            self.find_bsmooth(num_to_gen, tonelli, i)
            return
        
        self.i = i
        sq = int(math.sqrt(self.n))
        self.matrix = np.array([])
        limit = 10000

        while len(self.matrix) <= num_to_gen:
            starting_plus = sq + self.i
            starting_minus = sq - self.i
            self.modular_plus = dict({})
            self.modular_minus = dict({})
            for prime in self.tonelli_relations:
                self.modular_plus.update({prime: [starting_plus % (prime ** x) for x in range(1, len(self.tonelli_relations[prime]))]})
                self.modular_minus.update({prime: [starting_minus % (prime ** x) for x in range(1, len(self.tonelli_relations[prime]))]})
            possible_plus = [0] * limit
            possible_minus = [0] * limit
            for prime in self.tonelli_relations:
                current_log = self.tonelli_relations[prime][0]
                for count, relation in enumerate(self.tonelli_relations[prime][1:]):
                    modulo = prime ** (count + 1)
                    plus_start = self.modular_plus[prime][count]
                    minus_start = self.modular_minus[prime][count]
                    for rel in relation:
                        plus_offset = rel - plus_start
                        if plus_start > rel:
                            plus_offset += modulo

                        minus_offset = minus_start - rel
                        if minus_start < rel:
                            minus_offset += modulo

                        if plus_offset > minus_offset:
                            difference = plus_offset - minus_offset
                            while plus_offset < limit:
                                possible_plus[plus_offset] += current_log
                                possible_minus[plus_offset - difference] += current_log
                                plus_offset += modulo
                        else:
                            difference = minus_offset - plus_offset
                            while minus_offset < limit:
                                possible_minus[minus_offset] += current_log
                                possible_plus[minus_offset - difference] += current_log
                                minus_offset += modulo
            
            cutoff = 0.5 * math.log(self.n / 2) + math.log(starting_plus + limit - sq) - 1.6 * math.log(self.factor_base[-1])
            # From Silverman "The Multiple Polynomial Quadratic Sieve"
            candidates_p = []
            candidates_m = []
            for x in range(limit):
                if possible_plus[x] > cutoff:
                    temp_p = starting_plus + x
                    candidates_p.append(temp_p)
                if possible_minus[x] > cutoff:
                    temp_m = starting_minus - x
                    candidates_m.append(temp_m)   

            self.i = self.i + limit

            for candidate in candidates_p:
                factored, factors = self.factor_with_base(self.factor_base, candidate ** 2 - self.n)
                if factored == 1:
                    print(f'%d %d %d %s' % (candidate, candidate ** 2 - self.n, factored, factors))
                    if len(self.matrix) == 0:
                        self.matrix = np.array([factors])
                        self.bsmooth = np.array(candidate)
                    else:
                        self.matrix = np.append(self.matrix, [factors], axis=0)
                        self.bsmooth = np.append(self.bsmooth, candidate)
                    if len(self.matrix) > num_to_gen:
                        break

            for candidate in candidates_m:
                value = -1 * (candidate ** 2 - self.n)
                factored, factors = self.factor_with_base(self.factor_base, value)
                if factored == 1:
                    factors[0] = 1
                    print(f'%d %d %d %s' % (candidate, value, factored, factors))
                    if len(self.matrix) == 0:
                        self.matrix = np.array([factors])
                        self.bsmooth = np.array(candidate)
                    else:
                        self.matrix = np.append(self.matrix, [factors], axis=0)
                        self.bsmooth = np.append(self.bsmooth, candidate)

    # generates tonelli relations and the relations of powers of primes up to num_powerss
    def generate_tonelli(self, num_powers):
        for prime in self.factor_base:
            if prime == -1:
                continue
            self.tonelli_relations.update({prime: [math.log(prime), list(set(self.tonelli_shanks(self.n, prime, 1)))]})
            if prime == 2:
                for x in range(5 * num_powers):
                    self.tonelli_relations[prime].append(list(set(self.tonelli_shanks(self.n, prime, x + 2))))    
            else:
                for x in range(num_powers):
                    self.tonelli_relations[prime].append(list(set(self.tonelli_shanks(self.n, prime, x + 2))))

    # this function generates the limit B
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
                        linear_combinations[i][1] = filter_array(np.append(linear_combinations[i][1], linear_combinations[h][1]))
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
                linear_dependencies.append(filter_array(linear_combs[i][1])) # this step is taking a very long time
                # each array in linear_combs is like 3 million elements long
        return linear_dependencies

    # given new exponent vectors `newrows` and new bsmooth numbers `newbsmooth`, appropriately concatenate them with the old arrays
    # M: exponent vectors, A : exponent vectors reduced mod 2 (with old rows row-reduced), v: bsmooth numbers
    def initialize_objects(self, newrows, newbsmooth):
        nr = newrows % 2 != 0
        linear_combinations = self.lincombs
        
        #print("Initializing")
        if self.old_matrix is not None: 
            A = np.concatenate((self.old_matrix, nr)) #A contains the exponent vectors mod 2 to be row reduced
        else: 
            A = nr

        if self.old_rows is not None : 
            M = np.concatenate((self.old_rows, newrows)) #M is the complete (unreduced) 'exponent vector' matrix
        else: 
            M = newrows

        if self.old_bsmooth is not None : 
            v = np.concatenate((self.old_bsmooth, newbsmooth)) #v is a vector that contains all the bsmooth numbers
        else: 
            v = newbsmooth

        m , n = np.shape(A)
        l = len(linear_combinations)
        for i in range(len(newbsmooth)):
            linear_combinations[l + i] = [newbsmooth[i], [newbsmooth[i]]]
        return A, M, v, linear_combinations
    
    #newrows: new exponent vectors (2d array). dimensions (m x n)
    #newbsmooth: new bsmooth numbers corresponding to new rows. dimensions (m)
    #returns two arrays C and B. C[i]**2 is congruent to B[i]**2 mod n for all i. 
    def find_congruent_squares(self, newrows, newbsmooth): 
        A, M, v, linear_combinations = self.initialize_objects(newrows, newbsmooth)
        A, linear_combinations = self.row_reduce(A, linear_combinations)
        linear_dependencies = self.find_lin_dependencies(A, linear_combinations)
        #tracks information for next call to `lin_dep_mod2`
        self.old_matrix = A 
        self.old_bsmooth = v
        self.lincombs = linear_combinations
        self.old_rows = M

        B = []
        #for each linear dependency, compute the corresponding product of prime powers (mod n) and store in an array B
        for dependency in linear_dependencies:
            indices = []
            for k in dependency:
                indices.append(np.where(v == k)[0])  #finding indices corresponding to B-smooth numbers
            prime_exp = sum([M[i]for i in indices])[0]
            #print(prime_exp)
            prime_exp = prime_exp // 2

            b = 1
            for i in range(len(self.factor_base)):
                b = (b* pow(self.factor_base[i], prime_exp[i].item(), self.n)) % self.n
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
            return 1
        else: 
            return math.gcd(abs(a-b), self.n)
    
    #driver code
    def find_prime_factor(self, tonelli=True):
        ret = []
        B = self.get_B()
        if self.n < 100000000000 or not tonelli:
            tonelli = False
            self.factor_base = self.gen_primes(B)
        else:
            self.factor_base = [-1] + self.eulers_criterion(self.gen_primes(B))
            self.generate_tonelli(2)
        print(self.factor_base)
        num_to_gen = len(self.factor_base) + 5
        while len(ret) == 0:
            self.better_find_bsmooth(num_to_gen, tonelli, self.i)
            A, C = self.find_congruent_squares(self.matrix, self.bsmooth)
            for i in range(len(A)):
                j = self.basic_principle(A[i], C[i])
                if j > 1:
                    ret.append(j)
            num_to_gen = 1
        return ret

#Sieve = QuadraticSieve(101 * 109)
#Sieve = QuadraticSieve(1093 * 3511)       
#Sieve = QuadraticSieve(8101 * 9547)
#Sieve = QuadraticSieve(100109 * 100271)
#Sieve = QuadraticSieve(100109 * 386429)
#Sieve = QuadraticSieve(100271 * 5009317)
#Sieve = QuadraticSieve(10000019 * 1000003)
#Sieve = QuadraticSieve(310248241 * 383838383)
#Sieve = QuadraticSieve(2860486313 * 5915587277)            # first test case - 16921456439215439701
Sieve = QuadraticSieve(100123456789 * 1012346665879)
#Sieve = QuadraticSieve(46839566299936919234246726809)      # second test case

current = time.time()
res = Sieve.find_prime_factor(tonelli=True)
end = time.time()
print(res)
#print(f'Tonelli took %f seconds' % (Sieve.tonelli_time))
print(f'Total took %f seconds' % (end - current))
