import random
from pyscipopt import Model
from itertools import combinations
from itertools import product
from pyscipopt import quicksum
import numpy as np
import time
from tqdm import tqdm 
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

    def is_feasible(self, T, job_pair = [-1, -1], columnGeneration=False):
        
        """
        :param T: integer
        :return: True if LP(T) is feasible, otherwise False
        """
        model = Model('CONFIGURATION LP')
        model.hideOutput()
        # Determining valid configurations (with makespan at most T) for each machine in the form of a dict
        # Values are lists of tuples, one tuple for each valid configuration
        # When a job pair is specified, we leave out all configs containing both jobs.
        configs = [c for length in range(1, self.n_jobs + 1) for c in combinations(range(self.n_jobs), length) if
                sum(int(self.p_times[j]) for j in c) <= T and not (job_pair[0] in c and job_pair[1] in c)]

        

        if (not columnGeneration):
           
            
            # Decision variables
            x = {}
            for c in configs:
                x[c] = model.addVar(vtype="C", name=f"x({c})", lb=0.0)
            # The sum of the variables is at most 2.
            model.addCons(quicksum(x[c] for c in configs) <= self.n_machines)

            # Each job gets allocated at least once
            for j in range(self.n_jobs):
                model.addCons(quicksum(x[c] for c in configs if j in c) >= 1)

            model.optimize()

            if model.getStatus() == 'optimal':
                x_val = dict(zip(x.keys(), [model.getVal(x[e]) for e in x.keys()]))
                return True, x_val
            return False, {}


        # Decision variables
        x = {}
        #constraint
        s = []

        group = []
        for j in range(self.n_jobs):
            s.append(0)
            group.append(j)

        #return self.columnGen(x,s,model,T,group)

        total_machines = 0
        random.shuffle(configs)
        #configs.reverse()
        for c in tqdm(configs):
            
            model.freeTransform()
            x[c] = model.addVar(vtype="C", name=f"x({c})", lb=0.0)


            # The sum of the variables is at most 2.
            if not total_machines:
                total_machines = model.addCons((x[c]) <= self.n_machines)
            else:
                model.addConsCoeff(total_machines,x[c],1)


            # Each job gets allocated at least once
            check = False
            for j in range(self.n_jobs):
                if j in c:
                    if s[j] == 0:
                        s[j] = model.addCons(x[c]>=1)
                        check = True
                    else: 
                        model.addConsCoeff(s[j], x[c],1)
            for j in range(self.n_jobs):
                if s[j] == 0:
                    check = True

            if check: continue

            model.optimize()

            if model.getStatus() == 'optimal':
                x_val = dict(zip(x.keys(), [model.getVal(x[e]) for e in x.keys()]))
                return True, x_val
        return False, {}
    

    #deprecated, nt
    def columnGen(self,x,s, model,T,total_machines, group = []):

        if group:
            result = self.columnGen(x,s,model,T,group[:-1])
            if result and len(result)>0:
                for partial in result:
                    if len(partial)==0: continue
                    print("AAAAAAAAAA",partial)
                    model.freeTransform()
                    x[partial] = model.addVar(vtype="C", name=f"x({partial})", lb=0.0)


                    # The sum of the variables is at most 2.
                    if not total_machines:
                        total_machines = model.addCons((x[partial]) <= self.n_machines)
                    else:
                        model.addConsCoeff(total_machines,x[group],1)


                    # Each job gets allocated at least once
                    check = False
                    for j in range(self.n_jobs):
                        if j in partial:
                            if s[j] == 0:
                                s[j] = model.addCons(x[partial]>=1)
                                check = True
                            else: 
                                model.addConsCoeff(s[j], x[partial],1)
                    for j in range(self.n_jobs):
                        if s[j] == 0:
                            check = True

                    if not check: return result
                    model.optimize()

                    if model.getStatus() == 'optimal':
                        x_val = dict(zip(x.keys(), [model.getVal(x[e]) for e in x.keys()]))
                        return True, x_val
            return result + [c + [group[-1]] for c in result]
            
            
        else:
            return [[]]

    
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
        #w = {}
        for c in configs:
            x[c] = model.addVar(vtype="I", name=f"2x({c})", lb=0.0)
            #w[c] = model.addVar(vtype="I", name=f"w({c})", lb=0.0)
        

        # The sum of the variables is at most 2.
        model.addCons(quicksum(x[c]/2 for c in configs) <= self.n_machines)
    
        # half integrality
        for c in configs:
            #model.addCons(2*x[c]==w[c])
            model.addCons(x[c]<=2)

        # Each job gets allocated at least once
        for j in range(self.n_jobs):
            model.addCons(quicksum(x[c]/2 for c in configs if j in c) == 1)

        model.optimize()

        if model.getStatus() == 'optimal':
            x_val = dict(zip(x.keys(), [model.getVal(x[e]/2) for e in x.keys()]))
            return True, x_val
        return False, {}
    
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
            model.addCons(quicksum(x[i, j] for i in range(self.n_machines)) == 1)

        # Constraint 2. The processing time on each machine must be at most C_max
        for i in range(self.n_machines):
            model.addCons(quicksum(x[i, j] * int(self.p_times[j]) for j in range(self.n_jobs)) <= C_max)

        # Print the model
        # model.writeProblem('model.lp')

        # Optimize the model
        model.optimize()
        # model.freeTransform()

        solution = list({j: i for j in range(self.n_jobs) for i in range(self.n_machines) if model.getVal(x[i, j]) > 0.5}.values())

        return solution, model.getObjVal()

    
    
    def solve(self, T, job_pair=[-1, -1]):
        print("to implement")

    def find_half_integral_persistence(self,dictionary, column_generation=True, verbose=False):
        """
        :param dictionary: the feasible solution
        :return: True if there is any arbitrary reformulation of y that is half-integral given a persistent, feasble configuration, otherwise False
        """
        clean = {} #dictionary with only configurations that has positive non-zero value 
        if(not dictionary):
            print("not feasible")
            return False
        if verbose :print("original y= \n")
        for k in dictionary.keys():
            if (dictionary.get(k)>0.0000000000000001):
                clean[k] = dictionary.get(k)
                if verbose:print("subset=",k,"ibrid machine=",dictionary.get(k))

        complete = {} #dictionary for each job the % of completeness in {0,0.5,1}
    

        for set in clean.keys():
            clean[set] = 0
        
        for job in range(self.n_jobs):
            complete[job] = 0
        
        n_subsets= len(clean.keys())

        if column_generation:
            return self.config_find(0,clean,complete,0,0)

        for solution in product([0,0.5,1],repeat=n_subsets):
            i = 0
            for  key in clean.keys():
                clean[key] = solution[i]
                i+=1
            
            for s in clean.keys():
                value = clean[s]
                for j in s:
                    complete[j] += value

    
            check = True
            total_power=0
            for s in clean.keys():
                total_power += clean[s]
            
            if (total_power> self.n_machines):
                check = False
            for job in range(self.n_jobs):
                if (complete[job]<1):
                    check = False
            
            if (check):
                return clean
            
            for job in range(self.n_jobs):
                complete[job] = 0
            
        return False


    #with column generation
    def config_find(self, v, clean, complete, power, step):

        if v>=1.5:
            return False

        check = True
        total_power=0
        for s in clean.keys():
            total_power += clean[s]
        
        if (total_power> self.n_machines):
            return False
        
        for job in range(self.n_jobs):
            if (complete[job]<1):
                check = False
        
        if (check):
            return clean
    
        elif(step==len(list(clean.keys()))):
            return False
        

        p = power
        clean_cop = clean.copy()
        complete_cop = complete.copy()
        subset = list(clean_cop.keys())[step] #
        clean_cop[subset] = v
        p+=v
        for j in subset:
            complete_cop[j] += v

        

        return self.config_find(v+0.5,clean,complete,power, step) or self.config_find(0,clean_cop,complete_cop,p, step+1)


    
    def gap(self):
        return self.opt_IP()[1] / self.opt_LP()

