from util import *
import engine
from GurobiSP import solveSP
import time
from extra import get_general_timelist_and_anneal

def get_prob_num(K, timelimit):
    if K == 2000 and timelimit == 300: return 1
    if K == 2000 and timelimit == 480: return 2
    if K == 1000 and timelimit ==  60: return 3
    if K == 1000 and timelimit == 180: return 4
    if K == 1000 and timelimit == 300: return 5
    if K ==  750 and timelimit ==  30: return 6
    if K ==  750 and timelimit ==  60: return 7
    if K ==  500 and timelimit ==  15: return 8
    if K ==  500 and timelimit ==  30: return 9
    if K ==  300 and timelimit ==  15: return 10
    return 0

def get_timelist_and_anneal(prob_num, K, timelimit):
    if prob_num ==  0 :
        return get_general_timelist_and_anneal(K, timelimit)
    data_dict = {
        1: [[[195, 40], [10, 15], [15, 15]], 1.003, 3.5],
        2: [[[240, 60], [25, 30], [25, 30], [30, 30]], 1.03, 3.5],
        3: [[[25, 15], [10, 10]], 1.12, 1],
        4: [[[15, 15], [15, 15], [15, 15], [15, 15], [15, 15], [15, 15]], 1.04, 3],
        5: [[[180, 60], [15, 15], [15, 15]], 1.03, 3],
        6: [[[7.5, 7.5], [6.5, 8.5]], 1.07, 0.5],
        7: [[[20, 20], [8, 12]], 1.1, 0.5],
        8: [[[6.2, 8.8]], 1.1, 0],
        9: [[[7.5, 7.5], [6.5, 8.5]], 1.12, 0],
        10: [[[9.5, 5.5]], 1.1, 0.0]
    }
    return data_dict[prob_num][0], data_dict[prob_num][1], data_dict[prob_num][2]

def get_use_power(prob_num):
    return prob_num not in [3,4]

def get_wt_tw(prob_num):
    if prob_num in [3,10]:
        return [0.75, 0.75]
    else:
        return [0.2, 1.0]

def get_use_worst(prob_num):
    return (prob_num != 5)

def get_first_heur_time(prob_num):
    if prob_num ==5:
        return 3
    return 0

def get_use_old_version(prob_num):
    return prob_num in [9]

def algorithm(K:int, all_orders:list[Order], all_riders:list[Rider], dist_mat:np.ndarray[np.int32], timelimit=60, hparam=None)->list[str, list[int], list[int]]:
    # 0 : no logging , 1 : compact logging , 2 : all logging
    verbose = 0
    seed = 1
    
    start_time = time.time()
    ##################### hparam setting for problems #####################
    if not hparam: hparam = dict()

    prob_num = get_prob_num(K, timelimit)

    timelist, annealer_modify_ratio, SPARE_TIME = get_timelist_and_anneal(prob_num, K, timelimit)
    if K >= 2000: timelimit -= 10 # 서버 cpu timelimit 계속 나옴?

    use_power = get_use_power(prob_num)

    wt_tw = get_wt_tw(prob_num)
    hparam["wt_weight"] = wt_tw[0]
    hparam["tw_weight"] = wt_tw[1]
    
    use_worst = get_use_worst(prob_num)
    hparam["use_worst"] = use_worst

    first_heur_time = get_first_heur_time(prob_num)

    use_old = get_use_old_version(prob_num)
    hparam["use_old"]=use_old

    ##################### Data preparation #####################

    orders, riders, riders_available = preprocess(all_orders, all_riders, dist_mat)
    walk_T, walk_info = riders["WALK"].extract()
    bike_T, bike_info = riders["BIKE"].extract()
    car_T, car_info = riders["CAR"].extract()
    
    sol_edge = np.zeros((2*K, 2*K), dtype=int)
    sol_rider = np.zeros(K, dtype=int)
    init_bundle = []

    ##################### Start solving #####################

    if verbose >= 1: print("############################# Start #############################")
    if verbose >= 1: print(f"timelist={timelist}, modify_ratio : {annealer_modify_ratio}")
    
    for i in range(len(timelist)):
        if verbose >= 1: print(f"################## iter {i+1} #######################")
        ########### Time info ####################################
        ALNS_TIME = timelist[i][0]
        SP_TIME = timelist[i][1]
        
        if i==len(timelist)-1:
            lefttime = timelimit - (time.time() - start_time)
            SP_TIME = timelist[i][1]
            ALNS_TIME = lefttime - SP_TIME - SPARE_TIME # 사실 운영진이 여유시간 1초 주는거 같은디
            
        ########### ALNS ####################################
        if verbose >= 1: print("entering ALNS, spent time:", time.time() - start_time)
        result = engine.run(ALNS_TIME, annealer_modify_ratio, seed,
                            sol_edge, sol_rider, init_bundle,
                            orders, walk_T, walk_info, bike_T, bike_info, car_T, car_info, riders_available, dist_mat, hparam, use_power,
                            verbose==2)
        
        best_sol_without_SP_score = result[0][0]
        best_sol_without_SP = result[0][1]

        initial_sol_for_gurobi = result[1][0]
        constraint_helper_for_gurobi = result[1][1]
        bundle_vectors_for_gurobi = result[1][2]
        if verbose >= 1: print("ALNS result : ", best_sol_without_SP_score/K)
        if verbose >= 1: print("bundles generated in total : ", len(bundle_vectors_for_gurobi))

        ########### SP ####################################

        if verbose >= 1: print("entering SP spent time", time.time() - start_time)
        heur_time = first_heur_time if i == 0 else -1
        best_sol, sol_edge, sol_rider, init_bundle, gurobi_score = solveSP(initial_sol_for_gurobi, constraint_helper_for_gurobi, riders_available, bundle_vectors_for_gurobi, SP_TIME, verbose, dist_mat, riders, heur_time)
        if verbose >= 1: print("GUROBI result(without SC->SP conversion) : ", gurobi_score)
    
    return best_sol


  
# Order : [.id, .order_time, .cook_time, .volume, .deadline]
# Rider : [.type, .speed, .capa, .var_cost, .fixed_cost, .service_time, .available_number]
# dist_mat : 2K * 2K int matrix
# Every type is an Integer
class Order_Info():
    def __init__(self, order:Order):
        self.start = order.order_time + order.cook_time # int
        self.end = order.deadline # int
        self.volume = order.volume # int
        
class Rider_Info():
    def __init__(self, rider:Rider, dist_mat):
        self.T = np.round(dist_mat / rider.speed + rider.service_time).astype(int)
        self.capa = rider.capa # int

        self.fixed_cost = rider.fixed_cost # int 
        self.var_cost = rider.var_cost # int

    def cost(self, dist:int)->float:
        return self.fixed_cost + dist/100 * self.var_cost
    
    def extract(self)->tuple[np.ndarray, list[int]]:
        return self.T, [self.capa, self.fixed_cost, self.var_cost]

def preprocess(all_orders, all_riders, dist_mat)->tuple[dict[int, Order_Info], dict[str, Rider_Info], dict[str, int]]:
    orders = dict()
    orders_array = []
    for order in all_orders:
        info = Order_Info(order)
        orders[order.id] = info
        orders_array.append([info.start, info.end, info.volume])
    riders = dict()
    riders_available = dict()
    for rider in all_riders:
        riders[rider.type] = Rider_Info(rider, dist_mat)
        riders_available[rider.type] = rider.available_number

    return np.array(orders_array), riders, riders_available