import numpy as np
import math 

class QuadraticSieve:
    
    def __init__(self, n):
        self.n = n
        self.matrix = np.array([[]])

    #TODO: Jack
    def find_bsmooth(self, B):
        return 
    
    def get_B(self):
        B = np.exp((1/2)*math.sqrt(math.log(self.n)*math.log(math.log(self.n))))
        return math.ceil(B)
    
    #TODO: Stavan
    def find_nullspace(self):
        return 

    def basic_principle(self, a, b):
        if((a-b)%self.n==0 or (a+b)%self.n==0 or (b-a)%self.n==0):
            return False
        else: 
            return math.gcd(abs(a-b), self.n)