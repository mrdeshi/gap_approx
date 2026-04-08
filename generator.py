from covering import Covering
import numpy as np
import random

class Generator:
    def __init__(self, population=range(100)):
        self.population = population
        pass

    def pick(self,n):
        jobs = random.sample(self.population,n)
        return jobs

    def gen(self, w=100):

        m_step = 5
        m = 0
        
        ins = []
        for i in range(w):
            m= random.randint(2,5)
            n_jobs = random.randint(m_step,15)
            ins.append(Covering(m,n_jobs,self.pick(n_jobs)))
        return ins
    

    