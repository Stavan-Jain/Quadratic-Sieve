import numpy as np
from utils import *

class QuadraticSieve:
    
    def __init__(self, n):
        self.n = n
        self.matrix = np.array([[]])

    #TODO: Jack
    def find_bsmooth(self, B):
        return 
    
    #TODO: Prerana
    def get_B(self):
        return 
    
    #TODO: Stavan
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
        linear_dependencies = []
        for i in range(m):
            if True not in A[i]:
                linear_dependencies.append(filter_array(linear_combinations[i][1]))

        B = []
        for dependency in linear_dependencies:
            indices = []
            b = 1
            for k in dependency:
                indices.append(np.where(v == k)[0])
            prime_exp = sum([M[i]for i in indices])[0]
            prime_exp = np.array([int(prime / 2) for prime in prime_exp])
            print(prime_exp)
            for i in range(n):
                b *= factor_base[i]**prime_exp[i]

            print(b % self.n)
            B.append(b % self.n)
        A = [np.prod(arr) % self.n for arr in linear_dependencies]
        return A, B

    #TODO: Prerana
    def basic_principle(self):
        return False
    
Sieve = QuadraticSieve(3837523)
print(Sieve.find_linear_dependency_mod_2(M1, v1, base))

