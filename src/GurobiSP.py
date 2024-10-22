from util import *
import engine
import numpy as np

from itertools import combinations
import gurobipy as gp
from gurobipy import GRB

INT_INF = int(1e9)

class GUROBI_SOLVER:    
    def __init__(self, initial_sol, helper_matrix, riders_available, bundle_costs, bundle_riders, timelimit, heur_time, verbose):

        self.verbose = verbose
        self.timelimit = timelimit
        self.heur_time = heur_time

        self.initial_sol = initial_sol
        self.helper_matrix = helper_matrix
        self.riders_available = riders_available
        self.bundle_costs = bundle_costs
        self.bundle_riders = bundle_riders

        self.n = len(bundle_costs) # number of bundles
        self.K = len(helper_matrix) # number of orders
        self.r = len(riders_available) # number of riders
        
        self.prepare()
        
        
    def prepare(self):
        start_time = time.time()
        
        env = gp.Env(empty=True)
        env.setParam('OutputFlag', self.verbose==2)
        env.start()
        
        self.model = gp.Model(env=env)
        if self.verbose == 2: print("License check done!\n")
        
        ################## Set Variables ##################
        x = self.model.addVars(self.n, vtype=GRB.BINARY, name='x')
        for i in range(self.n): x[i].start = self.initial_sol[i]
        self.sol = x

        ################## Set Objective ##################
        self.model.setObjective(gp.quicksum(self.bundle_costs[i] * x[i] / self.K for i in range(self.n)), GRB.MINIMIZE)
        
        ################## Set Constraints ##################
        # (2)
        for idx, L in enumerate(self.helper_matrix):
            self.model.addConstr(gp.quicksum(x[i] for i in L) >= 1, name=f"constraint_M_{idx}")
        # (3)
        self.model.addConstr(gp.quicksum(x[i] for i in range(self.n) if self.bundle_riders[i] == 0) <= self.riders_available[0], name="constraint_type_1")
        self.model.addConstr(gp.quicksum(x[i] for i in range(self.n) if self.bundle_riders[i] == 1) <= self.riders_available[1], name="constraint_type_2")
        self.model.addConstr(gp.quicksum(x[i] for i in range(self.n) if self.bundle_riders[i] == 2) <= self.riders_available[2], name="constraint_type_3")
        
        spent_time = time.time() - start_time
        ################## set model Params ##################
        self.model.Params.MIPGap = 0                                  # Relative MIP optimality gap
        self.model.Params.TimeLimit = self.timelimit - spent_time     # Time limit
        self.model.Params.Threads = 4                                 # Number of parallel threads to use
        self.model.Params.Seed = 42                                   # Random number seed
        self.model.Params.LogToConsole = (self.verbose==2)            # Console logging
        self.model.Params.MIPFocus = 1
        self.model.Params.NoRelHeurTime = self.heur_time if self.heur_time != -1 else self.model.Params.TimeLimit
        self.model.Params.SimplexPricing = 2
        #self.model.Params.Heuristics = 1 if self.heur_time == -1 else 0.05
        #self.model.Params.ImproveStartGap = 1
        #self.model.Params.ImproveStartGap = 1e+30  # Delay bound improvement efforts
        #self.model.Params.Cuts = 1             # Disable cuts to save time (optional)
        #self.model.Params.Presolve = 1  
        ######################################################
            
    
    def solve(self):
        self.model.optimize()
        if self.verbose == 2: print("Objective value:", self.model.ObjVal)
        if self.verbose >= 1: print(f"Gurobi gap : {self.model.MIPGap}")
        
        sol = np.array([self.sol[i].X > 0.5 for i in range(self.n)])
        return sol
  
    
def gurobi_sol_to_answer(gurobi_sol, bundle_vectors, dist_mat, riders):
    
    K = dist_mat.shape[0] // 2
    set_cover_sol = [[bundle_vectors[i][0], bundle_vectors[i][2], bundle_vectors[i][3]] for i in np.nonzero(gurobi_sol)[0]]
    
    cnt = np.zeros(K, dtype=int)
    for i in range(len(set_cover_sol)):
        source = set_cover_sol[i][1]
        for x in source: cnt[x] += 1

    order_pool = np.array([i for i, x in enumerate(cnt) if x > 1])
    np.random.shuffle(order_pool)
    
    for order in order_pool:
        best_cost = np.inf
        where_to_keep = 0
        
        for bundle_id, bundle in enumerate(set_cover_sol):
            rider_type, source, dest = bundle
            if order not in source: continue
            capa, fixed_cost, var_cost = riders[rider_type].extract()[1]
            
            n = len(source)
            if n == 1:
                dist_incremental = dist_mat[order, order+K]
                cost_incremental = fixed_cost + dist_incremental/100 * var_cost
            
            elif source[n-1] == order and dest[0] == order:
                bef, aft = source[n-2], dest[1]+K
                dist_incremental = dist_mat[bef,order] + dist_mat[order, order+K] + dist_mat[order+K, aft] - dist_mat[bef,aft]
                cost_incremental = dist_incremental/100 * var_cost
            
            else:
                dist_incremental = 0
                
                i = source.index(order)
                if i == 0:
                    dist_incremental += dist_mat[source[i], source[i+1]]
                elif i == n-1:
                    dist_incremental += dist_mat[source[i-1], source[i]] + dist_mat[source[i], dest[0]+K] - dist_mat[source[i-1], dest[0]+K]
                else:
                    dist_incremental += dist_mat[source[i-1], source[i]] + dist_mat[source[i], source[i+1]] - dist_mat[source[i-1], source[i+1]]
                
                i = dest.index(order)
                if i == n-1:
                    dist_incremental += dist_mat[dest[i-1], dest[i]]
                elif i == 0:
                    dist_incremental += dist_mat[source[n-1], dest[i]+K] + dist_mat[dest[i]+K, dest[i+1]+K] - dist_mat[source[n-1], dest[i+1]+K]
                else:
                    dist_incremental += dist_mat[dest[i-1]+K, dest[i]+K] + dist_mat[dest[i]+K, dest[i+1]+K] - dist_mat[dest[i-1]+K, dest[i+1]+K]
                    
                cost_incremental = dist_incremental/100 * var_cost
                
            if cost_incremental < best_cost:
                best_cost = cost_incremental
                where_to_keep = bundle_id
                
            
        for bundle_id, bundle in enumerate(set_cover_sol):
            rider_type, source, dest = bundle
            if order not in source: continue
            if where_to_keep == bundle_id: continue
            
            source.remove(order)
            dest.remove(order)

            set_cover_sol[bundle_id] = [rider_type, source, dest]
        
    best_sol = []
    for rider_type, source, dest in set_cover_sol:
        if len(source) == 0: continue
        best_sol.append([rider_type, source, dest])
    return best_sol
    
    

def solveSP(initial_sol_for_gurobi, constraint_helper_for_gurobi, riders_available, bundle_vectors_for_gurobi, timelimit, verbose, dist_mat, riders, heur_time)->list[str, list[int], list[int]]:
    start_time = time.time()
    K = len(dist_mat)//2
    initial_sol_for_gurobi_processed = np.zeros(len(bundle_vectors_for_gurobi), dtype=bool)
    initial_sol_for_gurobi_processed[initial_sol_for_gurobi] = True
    
    bundle_costs = np.array([x[1] for x in bundle_vectors_for_gurobi])
    
    rider_dict = { "WALK" : 0, "BIKE" : 1, "CAR" : 2}
    def convert(s): return rider_dict[s]

    bundle_riders = np.array([convert(x[0]) for x in bundle_vectors_for_gurobi])
    riders_available_processed = {convert(k) : v for k,v in riders_available.items()}

    solver = GUROBI_SOLVER(initial_sol_for_gurobi_processed, constraint_helper_for_gurobi, riders_available_processed, bundle_costs, bundle_riders, timelimit, heur_time, verbose)
    gurobi_sol = solver.solve()
    
    best_sol = gurobi_sol_to_answer(gurobi_sol, bundle_vectors_for_gurobi, dist_mat, riders)
    
    ## prepare for next iter ##
    sol_edge = np.zeros((2*K, 2*K), dtype=int)
    sol_rider = np.zeros(K, dtype=int)
    init_bundle = bundle_vectors_for_gurobi
    
    for rider_type, source, dest in best_sol:
        n = len(source)
        
        for i in range(n-1): sol_edge[source[i], source[i+1]] = 1
        sol_edge[source[n-1], dest[0]+K] = 1
        for i in range(n-1): sol_edge[dest[i]+K, dest[i+1]+K] = 1
        
        for x in source: sol_rider[x] = rider_dict[rider_type]

    return best_sol, sol_edge, sol_rider, init_bundle, solver.model.ObjVal