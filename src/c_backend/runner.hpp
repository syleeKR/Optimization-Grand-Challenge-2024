
#ifndef RUNNER_HPP
#define RUNNER_HPP

#include "helper_classes.hpp"
#include "helper_functions.hpp"
#include "destroy_solution.hpp"
#include "annealer.hpp"
#include "utils.hpp"
#include "hyperparam.hpp"
#include "main_iteration.hpp"

#include <vector>
#include <bits/stdc++.h>
#include <omp.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
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



pair<SOLUTION_FORMAT, GUROBI_INPUT> run(
    double total_timelimit,
    float annealer_modify_ratio,
    int seed,
    py::array_t<int> & sol_edge, // (2*K, 2*K)
    py::array_t<int> & sol_rider, // (K,)
    vector<Bundle> init_bundle,
    py::array_t<int> & orders, // K by 3
    py::array_t<int> & walk_T, // 2K by 2K
    vint & walk_info,
    py::array_t<int> & bike_T,
    vint & bike_info,
    py::array_t<int> & car_T,
    vint & car_info,
    map<string, int> & riders_available, 
    py::array_t<int> & dist_mat,
    py::dict py_hparam,
    bool use_power,
    bool verbose = false)
{    
    /* ##################################### read in inputs ##################################### */
    auto start_time = chrono::high_resolution_clock::now();
    omp_set_num_threads(4);

    py::buffer_info walk_T_buffer = walk_T.request();
    py::buffer_info bike_T_buffer = bike_T.request();
    py::buffer_info car_T_buffer = car_T.request();

    py::buffer_info dist_mat_buffer = dist_mat.request();
    py::buffer_info orders_buffer = orders.request();

    int* walk_T_ptr = static_cast<int*>(walk_T_buffer.ptr);
    int* bike_T_ptr = static_cast<int*>(bike_T_buffer.ptr);
    int* car_T_ptr = static_cast<int*>(car_T_buffer.ptr);

    int* dist_mat_ptr = static_cast<int*>(dist_mat_buffer.ptr);
    int* orders_ptr = static_cast<int*>(orders_buffer.ptr);
    size_t l = car_T_buffer.shape[0]; // this is equal to 2K;
    int K = l/2;
    
    py::buffer_info sol_edge_buffer = sol_edge.request();
    py::buffer_info sol_rider_buffer = sol_rider.request();
    
    int* sol_edge_ptr = static_cast<int*>(sol_edge_buffer.ptr);
    int* sol_rider_ptr = static_cast<int*>(sol_rider_buffer.ptr);
    
    auto cpp_hparam = py_hparam.cast<map<string,float>>();
    
    /* ##################################### initialize ##################################### */
    
    Hyperparameter hparam = Hyperparameter(cpp_hparam);
    
    bool first_alns = (sz(init_bundle)==0) ? true : false;
    Solution bestsol = Solution();
    Rider_Info rider_info = Rider_Info(walk_T_ptr, walk_info, bike_T_ptr, bike_info,  car_T_ptr, car_info);

    BundleStorage bundle_storage = BundleStorage(init_bundle);
    
    little_random_get_ready(orders_ptr, rider_info, dist_mat_ptr, l, hparam);

    AdaptiveDestroyer destroyer(K, orders_ptr,rider_info, dist_mat_ptr, hparam, seed);
    Builder builder(K, orders_ptr,rider_info, dist_mat_ptr, use_power, hparam, seed);
    
    /* ##################################### build initial solution ##################################### */
    if (!first_alns)
    {
        bestsol = input_to_solution(sol_edge_ptr, sol_rider_ptr, rider_info, dist_mat_ptr, riders_available, K);
    }
    else
    {
        vint allids; REP0(i, K) allids.pb(i);
        builder.build(allids, bestsol,riders_available);
    }
    bundle_storage.append(bestsol);
    
    if (verbose) {
        cout<<"time took for initial build : "<<get_time(start_time)<<endl;
        cout<<"initial cost : "<< bestsol.cost / K<<endl;}

    /* ##################################### initialize SA ##################################### */
    
    float temperature = (first_alns) ? bestsol.cost * 0.003 : 0;
    SimulatedAnnealing annealer(temperature, annealer_modify_ratio, K, seed);
    
    /* ##################################### Iterative update ##################################### */

    int iter = sequential_destroy_insert(
        total_timelimit,
        start_time, 
        bestsol, 
        riders_available,
        annealer,
        destroyer,
        builder,
        bundle_storage,
        hparam,
        verbose 
    );

    if(verbose) {
        cout<<"size of problem: "<<K<<", total iter: "<<iter<<endl;
    }
    // cout<<"size of problem: "<<K<<", total iter: "<<iter<<endl;
    // for(auto x: destroyer.probability_distribution)cout<<x<<" ";
    // cout<<endl;

    return {bestsol.extract(), bundle_storage.extract(bestsol, K, rider_info, dist_mat_ptr, l)};
}


#endif