#ifndef MAIN_ITERATION_HPP
#define MAIN_ITERATION_HPP


#include "helper_classes.hpp"
#include "helper_functions.hpp"
#include "destroy_solution.hpp"
#include "annealer.hpp"
#include "utils.hpp"
#include "repair.hpp"

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



int sequential_destroy_insert(
    double total_timelimit,
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time, 
    Solution & bestsol, 
    map<string, int> & riders_available,
    SimulatedAnnealing & annealer,
    AdaptiveDestroyer & destroyer,
    Builder & builder,
    BundleStorage & bundle_storage,
    Hyperparameter & hparam,
    bool verbose )
{
    Solution curr = bestsol;
    int iter = 0;
    int K = builder.K;

    while(true){
        double spent_time = get_time(start_time);
        if (spent_time > total_timelimit - hparam.time_margin)break;
        iter++;

        Solution newsol = curr;
        map<string, int> new_riders_available = riders_available;

        vint ids_destroyed = destroyer.destroy(newsol, new_riders_available);
        
        builder.build( ids_destroyed, newsol,new_riders_available);
        bundle_storage.append(newsol);

        int update_score = 0;
        if (newsol.cost < bestsol.cost)
        {            
            update_score = hparam.score1;
            bestsol = newsol;
            if (verbose) cout<<"current best : "<<bestsol.cost/K<<" on iter "<<iter<< " with method : "<<destroyer.current_method<<endl;
        }
        else if (newsol.cost < curr.cost)
            update_score = hparam.score2;
        else if (annealer.accept(newsol.cost, curr.cost, spent_time))
            update_score = hparam.score3;

        if (update_score != 0) {
            curr = newsol;
            riders_available = new_riders_available;
        }
        destroyer.update(update_score);
    }

    return iter;
}






#endif