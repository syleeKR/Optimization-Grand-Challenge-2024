#ifndef REPAIR_HPP
#define REPAIR_HPP

#include "helper_classes.hpp"
#include "helper_functions.hpp"
#include "repair_tools/optimize_rider.hpp"
#include "repair_tools/random_repair.hpp"
#include "repair_tools/little_random_repair.hpp"
#include "repair_tools/regret_repair.hpp"

#include <vector>
#include <bits/stdc++.h>
#include <omp.h>
namespace py = pybind11;
using namespace std;
#define REP0(i,n) for(int i=0; i<n; i++)
#define REP1(i,n) for(int i=1; i<=n; i++)
#define REP(i,a,b) for(int i=a; i<=b; i++)
#define sz(v) ((int)(v).size())
#define all(v) (v).begin(), (v).end()
#define compress(v) sort(all(v)); v.erase(    unique(all(v)) , v.end()   )
#define reset(X) memset(X, 0, sizeof(X))
#define pb push_back
#define fi first
#define se second
#define pii pair<int, int>
#define pll pair<ll, ll>
#define vint vector<int>
#define vll vector<ll>
#define vpii vector<pair<int, int>>
#define vpll vector<pair<ll,ll>>


struct Builder{

    Hyperparameter hparam;
    
    int K;
    int l;
    default_random_engine rng;
    bool use_power;
    
    int * orders_ptr; //orders array 1D
    Rider_Info rider_info; // info about each riders
    int * dist_mat_ptr; // distance matrix array 1D

    bool use_old; //use old version?

    Builder(
        int K, 
        int * orders_ptr,
        Rider_Info & rider_info,
        int * dist_mat_ptr_input,
        bool use_power,
        Hyperparameter & hparam,
        int seed=42
        ): hparam(hparam), K(K), use_power(use_power), orders_ptr(orders_ptr), rider_info(rider_info), dist_mat_ptr(dist_mat_ptr_input)
    {  
        rng.seed(seed);
        l = K *2;
        use_old = hparam.use_old;

    }

    void build(
        vint & ids_to_build, 
        Solution & solution, 
        map<string, int> & riders_available
        )
    {
        int n = sz(solution.solutions);
        
        std::vector<int> indices(n);
        for (int i = 0; i < n; ++i)indices[i] = i;
        std::shuffle(indices.begin(), indices.end(), rng);

        vint insert_order;
        for(int currsize = sz(ids_to_build); currsize>=1; currsize--)
        {
            std::uniform_int_distribution<std::size_t> distuniform(0, currsize - 1);
            int random_index = distuniform(rng);
            insert_order.pb(random_index);
        }
        
        optimize_rider_type(indices,solution, riders_available,rider_info,orders_ptr, dist_mat_ptr,l);
        if(!use_old)repair_little_random(insert_order, ids_to_build, solution, orders_ptr, rider_info, riders_available, dist_mat_ptr, l, use_power, rng, hparam.consider_size);
        if(use_old)repair_little_random_old(insert_order, ids_to_build, solution, orders_ptr, rider_info, riders_available, dist_mat_ptr, l, use_power, rng, hparam.consider_size);
    }
};

#endif