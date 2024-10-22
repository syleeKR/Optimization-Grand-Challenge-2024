#ifndef ANNEALER_HPP
#define ANNEALER_HPP

#include <bits/stdc++.h>
using namespace std;

struct SimulatedAnnealing{
    float T;
    float modify_ratio;
    int modify_cnt = 0;
    default_random_engine rng;
    int iter_per_second;

    SimulatedAnnealing(float initial_T, float modify_ratio, int K, int seed = 42):T(initial_T), modify_ratio(modify_ratio){
        iter_per_second = 200000 / K;
        rng.seed(seed);
    }

    bool accept(float cost_new, float cost_old, int iter) // need fix if time limit is larger
    {
        if (modify_cnt + 1 < iter/iter_per_second ) {
            modify_cnt++;
            modify_temp();
        }
        
        if(cost_new < cost_old)return true;
        else
        {
            float acceptance_prob = exp((cost_old - cost_new) / T);
            std::uniform_real_distribution<> dis(0.0, 1.0);
            float random_val = dis(rng);
            if (random_val < acceptance_prob)return true;
            else return false;
        }
    }
    void modify_temp()
    {
        T /= modify_ratio;
    }
};

# endif