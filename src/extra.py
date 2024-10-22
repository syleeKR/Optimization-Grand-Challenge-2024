
# customize this as you want(depending on your problem size and timelimit). 
# Important : Just a simple example(still.. decently good)!!!
def get_general_timelist_and_anneal(K, timelimit):
    timelist = [[timelimit/2, timelimit/2]]
    if 120 <= timelimit and timelimit <= 180:
        timelist = [[timelimit - 60, 30], [15,15]]
    elif 180 < timelimit:
        timelist = [[timelimit - 105, 45], [15,15],[15,15]]
    anneal = min(1.12, (6 ** 1/timelist[0][0]))
    # spare time : for timelimt, depends on the cpu performance
    spare_time = 1
    if K>=1000 and timelimit >= 180:
        spare_time =  3.5
    if K>=2000:
        timelist[0][0] -= 5
        timelist[0][1] -= 5


    return timelist, anneal, spare_time
