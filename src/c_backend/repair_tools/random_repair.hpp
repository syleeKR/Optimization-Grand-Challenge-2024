#ifndef RANDOM_REPAIR_HPP
#define RANDOM_REPAIR_HPP

#include "../helper_classes.hpp"
#include "../helper_functions.hpp"
#include "../hyperparam.hpp"

#include <vector>
#include <bits/stdc++.h>
#include <omp.h>
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


void repair_random(
    default_random_engine &rng,
    vint & ids_to_build, 
    Solution & solution, 
    int * orders_ptr,
    Rider_Info & rider_info,
    map<string, int> & riders_available, 
    int * dist_mat_ptr,
    int l
    )
{
    Hyperparameter hparam;
    
    while ( !ids_to_build.empty() )
    {
        std::uniform_int_distribution<std::size_t> distuniform(0, sz(ids_to_build) - 1);
        int random_index = distuniform(rng);
        int order_id_to_append  = ids_to_build[random_index];

        vector<tuple<float, int, string, float, vint, vint>> feasible_solutions(sz(solution.solutions)+1, {INF, -1, "", 0.0, {}, {}});
        
        // Case 1. to bundle
        # pragma omp parallel for
        for(int bundle_id = 0; bundle_id < sz(solution.solutions); bundle_id++)
        {
            float cost_before = solution.get_cost(bundle_id);

            Res res = investigate(order_id_to_append, solution.get_source(bundle_id), solution.get_dest(bundle_id), rider_info, orders_ptr, dist_mat_ptr, l);
            
            for(string rider_type : res.optimal_order) {
                if(riders_available[rider_type] > 0 || solution.get_rider_type(bundle_id) == rider_type) {
                    if (res.get_feasibility(rider_type)) {
                        float cost_after = res.get_cost(rider_type);
                        float cost_incremental = cost_after - cost_before;
                        
                        // 걍 대강 튜닝해둠. 수정 필요.
                        float capa = rider_info.prepare(rider_type).second[0];
                        for (int id : res.get_source(rider_type)) capa -= orders_ptr[f(id,2,3)];
                        float C = hparam.alpha1 * cost_incremental / 50 + hparam.alpha2 * capa;
                        
                        feasible_solutions[bundle_id] = make_tuple(C, bundle_id, rider_type, cost_after, res.get_source(rider_type), res.get_dest(rider_type));
                        break;
                    }
                }
            }
        }
            
        // Case 2. to empty
        Res res = investigate(order_id_to_append, {}, {}, rider_info, orders_ptr, dist_mat_ptr, l);
        
        for (string rider_type: res.optimal_order) {
            if (riders_available[rider_type] >0 && res.get_feasibility(rider_type)) {
                float cost_incremental = res.get_cost(rider_type);
                    
                // 걍 대강 튜닝해둠. 수정 필요.
                float capa = rider_info.prepare(rider_type).second[0];
                for (int id : res.get_source(rider_type)) capa -= orders_ptr[f(id,2,3)];
                float C = hparam.alpha1 * cost_incremental / 50 + hparam.alpha2 * capa;

                                                
                feasible_solutions[sz(solution.solutions)] = make_tuple(C, -1, rider_type,cost_incremental, res.get_source(rider_type), res.get_dest(rider_type));
                break;
            }
        }

        // remove it from ids_to_build
        auto it = std::find(all(ids_to_build), order_id_to_append); ids_to_build.erase(it);

        float cost_incre; int bundle_id; string rider_type; float newcost; vint newsource; vint newdest;
        sort(all(feasible_solutions));
        tie(cost_incre, bundle_id, rider_type, newcost, newsource, newdest) = feasible_solutions[0];

        // remove original bundle
        if (bundle_id != -1) {
            riders_available[solution.get_rider_type(bundle_id)]++;
            solution.remove(bundle_id);
        }

        // append inserted bundle
        solution.append(make_tuple(rider_type, newcost, newsource, newdest));
        riders_available[rider_type]--;
    }
}

#endif