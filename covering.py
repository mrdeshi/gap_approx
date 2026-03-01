from pyscipopt import Model
from itertools import combinations
import numpy as np
"""
We consider the case separately when there are just 2 machines. In that case, the problem is equivalent
with the covering polyhedra with an easier formulation. We can further assume that no machine-exclusive job exists.
"""

class Covering:
    def __init__(self, n_machines, n_jobs, p_times = None, p_max = 100, seed = 1):
        self.n_jobs = n_jobs
        self.n_machines = n_machines
        self.p_times = p_times if p_times else [np.random.randint(p_max) for _ in range(n_jobs)]
        self.p_max = p_max
        self.seed = seed

    def opt_LP(self, verbose=False, C_max = None):
        """
        Do a binary search to find the smallest integer for which is_feasible(T) is true.
        The initial guesses are as follows:
            l = highest processing time for any job - 1
            r = sum of all processing times
        We maintain the invariant that LB is in (l, r].

        Hint: as sometimes we compute the integer value, just to speed up the computation, we provide the C_max as "right" initial guess.
        """
        right = sum(self.p_times)
        left = sum(self.p_times)//self.n_machines - 1
        if C_max:
            right = min(right, C_max)
        _, x_keep = self.is_feasible(right)
        while right - left > 1:
            m = (left + right)//2
            is_feasible, x = self.is_feasible(m)
            if is_feasible:
                x_keep = x
                right = m
                if verbose:
                    print(f"LP feasible for T={m}")
            else:
                left = m

        return right, x_keep
    
    def half_integral_opt_LP(self, verbose=False, C_max = None):
        """
        Do a binary search to find the smallest integer for which is_half_integral_feasible(T) is true.
        The initial guesses are as follows:
            l = highest processing time for any job - 1
            r = sum of all processing times
        We maintain the invariant that LB is in (l, r].

        Hint: as sometimes we compute the integer value, just to speed up the computation, we provide the C_max as "right" initial guess.
        """
        right = sum(self.p_times)
        left = sum(self.p_times)//self.n_machines - 1
        if C_max:
            right = min(right, C_max)
        _, x_keep = self.is_half_integral_feasible(right)
        while right - left > 1:
            m = (left + right)//2
            is_half_integral_feasible, x = self.is_half_integral_feasible(m)
            if is_half_integral_feasible:
                x_keep = x
                right = m
                if verbose:
                    print(f"LP half integral feasible for T={m}")
            else:
                left = m

        return right, x_keep

    def is_feasible(self, T, job_pair = [-1, -1]):
        """
        :param T: integer
        :return: True if LP(T) is feasible, otherwise False
        """
        model = Model('Restricted assignment with 2 processing times')
        model.hideOutput()

        # Determining valid configurations (with makespan at most T) for each machine in the form of a dict
        # Values are lists of tuples, one tuple for each valid configuration
        # When a job pair is specified, we leave out all configs containing both jobs.
        configs = [c for length in range(1, self.n_jobs + 1) for c in combinations(range(self.n_jobs), length) if
                   sum(int(self.p_times[j]) for j in c) <= T and not (job_pair[0] in c and job_pair[1] in c)]

        # Decision variables
        x = {}
        for c in configs:
            x[c] = model.addVar(vtype="C", name=f"x({c})", lb=0.0)

        # The sum of the variables is at most 2.
        model.addCons(sum(x[c] for c in configs) <= self.n_machines)

        # Each job gets allocated at least once
        for j in range(self.n_jobs):
            model.addCons(sum(x[c] for c in configs if j in c) >= 1)

        model.optimize()

        if model.getStatus() == 'optimal':
            x_val = dict(zip(x.keys(), [model.getVal(x[e]) for e in x.keys()]))
            return True, x_val
        return False, {}
    
    def is_half_integral_feasible(self, T, job_pair = [-1, -1]):
        """
        :param T: integer
        :return: True if LP(T) is feasible, otherwise False
        """
        model = Model('Restricted assignment with 2 processing times')
        model.hideOutput()

        # Determining valid configurations (with makespan at most T) for each machine in the form of a dict
        # Values are lists of tuples, one tuple for each valid configuration
        # When a job pair is specified, we leave out all configs containing both jobs.
        configs = [c for length in range(1, self.n_jobs + 1) for c in combinations(range(self.n_jobs), length) if
                   sum(int(self.p_times[j]) for j in c) <= T and not (job_pair[0] in c and job_pair[1] in c)]

        # Decision variables
        x = {}
        w = {}
        for c in configs:
            x[c] = model.addVar(vtype="C", name=f"x({c})", lb=0.0)
            w[c] = model.addVar(vtype="I", name=f"w({c})", lb=0.0)
        

        # The sum of the variables is at most 2.
        model.addCons(sum(x[c] for c in configs) <= self.n_machines)
    
        # half integrality
        for c in configs:
            model.addCons(2*x[c]==w[c])
            model.addCons(w[c]<=2)

        # Each job gets allocated at least once
        for j in range(self.n_jobs):
            model.addCons(sum(x[c] for c in configs if j in c) >= 1)

        model.optimize()

        if model.getStatus() == 'optimal':
            x_val = dict(zip(x.keys(), [model.getVal(x[e]) for e in x.keys()]))
            return True, x_val
        return False, {}
    
    
    def solve(self, T, job_pair=[-1, -1]):
        print("to implement")



    def opt_IP(self, verbose=False):
        model = Model('Restricted assignment with 2 processing times')
        if not verbose:
            model.hideOutput()

        # Decision variables
        x = {}
        for i in range(self.n_machines):
            for j in range(self.n_jobs):
                x[i, j] = model.addVar(vtype="B", name=f"x({i},{j})", lb=0.0)

        # Makespan
        C_max = model.addVar(vtype="C", name="C_max", lb=0.0)

        # Objective function
        model.setObjective(C_max, "minimize")

        # Constraint 1. You have to allocate each job
        for j in range(self.n_jobs):
            model.addCons(sum(x[i, j] for i in range(self.n_machines)) == 1)

        # Constraint 2. The processing time on each machine must be at most C_max
        for i in range(self.n_machines):
            model.addCons(sum(x[i, j] * int(self.p_times[j]) for j in range(self.n_jobs)) <= C_max)

        # Print the model
        # model.writeProblem('model.lp')

        # Optimize the model
        model.optimize()
        # model.freeTransform()

        solution = list({j: i for j in range(self.n_jobs) for i in range(self.n_machines) if model.getVal(x[i, j]) > 0.5}.values())

        return solution, model.getObjVal()

    def gap(self):
        return self.opt_IP()[1] / self.opt_LP()

