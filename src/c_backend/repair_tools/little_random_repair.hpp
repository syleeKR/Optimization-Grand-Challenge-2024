#ifndef LITTLE_RANDOM_REPAIR
#define LITTLE_RANDOM_REPAIR

#include "../helper_classes.hpp"
#include "../helper_functions.hpp"
#include "../hyperparam.hpp"
#include "../utils.hpp"
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

map<string, vector<vector<int>>> promising;

void little_random_get_ready(int * orders_ptr, Rider_Info & rider_info, int * dist_mat_ptr, int l, Hyperparameter & hparam)
{
    int K = l/2;
    for (string rider_type : RIDER) {
        int* T = rider_info.get_T(rider_type);
        
        vector<vector<int>> rider_promising(K, vector<int>(K));
        for (int i = 0; i < K; ++i) {
            for (int j = 0; j < K; ++j) {
                if (i == j) continue;
                int dd1 = dist_mat_ptr[f(i,j,l)];
                int wt1 = max(orders_ptr[f(j,0,3)] - T[f(i,j,l)] - orders_ptr[f(i,1,3)], 0);
                int tw1 = max(orders_ptr[f(i,0,3)] + T[f(i,j,l)] - orders_ptr[f(j,1,3)], 0);
                int score1 = dd1 + (hparam.wt_weight * wt1 + hparam.tw_weight * tw1);
                
                int dd2 = dist_mat_ptr[f(i+K,j+K,l)];
                int wt2 = max(orders_ptr[f(j,0,3)] - T[f(i+K,j+K,l)] - orders_ptr[f(i,1,3)], 0);
                int tw2 = max(orders_ptr[f(i,0,3)] + T[f(i+K,j+K,l)] - orders_ptr[f(j,1,3)], 0);
                int score2 = dd2 +  (hparam.wt_weight * wt2 + hparam.tw_weight * tw2);
                
                rider_promising[i][j] = score1 + score2;
            }
        }
        promising[rider_type] = rider_promising;
    }
}

void repair_little_random(
    vint & insert_order,
    vint & ids_to_build, 
    Solution & solution, 
    int * orders_ptr,
    Rider_Info & rider_info,
    map<string, int> & riders_available, 
    int * dist_mat_ptr,
    int l,
    bool use_power,
    default_random_engine & rng,
    int consider_size
    )
{
    int n_iter = sz(ids_to_build);
    int maxcapacity = max(max(rider_info.walk_info[0], rider_info.bike_info[0]), rider_info.car_info[0]);
    for(int rep=0; rep<n_iter; rep++){
        int order_id_to_append = ids_to_build[insert_order[rep]];
    
        vector<tuple<float, int, string, float, vint, vint>> feasible_solutions(sz(solution.solutions)+1, {INF, -1, "", 0.0, {}, {}});
        
        int n = sz(solution.solutions);
        vector<pair<float,int>> vv(n, {INF, -1});
        
        # pragma omp parallel for
        for (int bundle_id = 0; bundle_id < n; bundle_id++) {
            string rider_type = solution.get_rider_type(bundle_id);
            int tmp = 0;
            bool ok =true;
            int filled =0 ;
            for (auto x : solution.get_source(bundle_id))
            {
                tmp += promising[rider_type][order_id_to_append][x];
                filled += orders_ptr[f(x,2,3)];
                if(orders_ptr[f(order_id_to_append,1,3)] < orders_ptr[f(x,0,3)]) {ok = false;break;}
                if(orders_ptr[f(x,1,3)] < orders_ptr[f(order_id_to_append,0,3)]) {ok = false;break;}
            }
            if (filled + orders_ptr[f(order_id_to_append, 2,3)] > maxcapacity)ok=false;
            if(ok)vv[bundle_id] = {(float)tmp / sz(solution.get_source(bundle_id)), bundle_id};
        }
        
        vint bundle_indices;
        
        if (use_power) {
            vector<pair<float,int>> vvv;
            REP0(i, n) if (vv[i].se != -1) vvv.pb(vv[i]);
            sort(all(vvv));

            int n_consideration = min(max(consider_size, l/100), sz(vvv));
            int nn = min(sz(vvv), 5*n_consideration);

            vint candidate(nn); REP0(i, nn) candidate[i] = vvv[i].se;

            bundle_indices.resize(n_consideration);

            std::uniform_real_distribution<> dis(0.0, 1.0);
            for (int i = 0; i < n_consideration; ++i) {
                int noise_idx = floor(pow(dis(rng), 5) * sz(candidate));
                int bundle_id = candidate[noise_idx];
                bundle_indices[i] = bundle_id;
                candidate.erase(candidate.begin() + noise_idx);
            }
        }
        else {
            vector<pair<float,int>> candidates;
            for (auto pair : vv) {
                if (pair.se == -1) continue;
                candidates.pb(pair);
            }
            sort(candidates.begin(), candidates.end());

            int n_consideration = min(max(consider_size, l/100), sz(candidates));

            // 4 core optimization for 750 : problem 6 and 7
            n_consideration = min((n_consideration + 3) / 4 * 4, sz(candidates));

            bundle_indices.resize(n_consideration);
            REP0(i,n_consideration)bundle_indices[i] = candidates[i].second;
        }

        // Case 1. to bundle
        # pragma omp parallel for
        for(int bundle_id : bundle_indices) 
        {
            float cost_before = solution.get_cost(bundle_id);

            Res res = investigate(order_id_to_append, solution.get_source(bundle_id), solution.get_dest(bundle_id), rider_info, orders_ptr, dist_mat_ptr, l);
            
            for(string rider_type : res.optimal_order) {
                if(riders_available[rider_type] > 0 || solution.get_rider_type(bundle_id) == rider_type) {
                    if (res.get_feasibility(rider_type)) {
                        float cost_after = res.get_cost(rider_type);
                        float cost_incremental = cost_after - cost_before;

                        feasible_solutions[bundle_id] = make_tuple(cost_incremental, bundle_id, rider_type, cost_after, res.get_source(rider_type), res.get_dest(rider_type));
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
                feasible_solutions[sz(solution.solutions)] = make_tuple(cost_incremental, -1, rider_type,cost_incremental, res.get_source(rider_type), res.get_dest(rider_type));
                break;
            }
        }

        // remove it from ids_to_build
        auto it = std::find(all(ids_to_build), order_id_to_append); ids_to_build.erase(it);

        float cost_incre; int bundle_id; string rider_type; float newcost; vint newsource; vint newdest;

        int minindex = 0; float minval = get<0>(feasible_solutions[0]);
        for(int i=1 ; i<sz(feasible_solutions); i++)
        {
            if ( get<0>(feasible_solutions[i]) < minval)
            {
                minindex = i;
                minval = get<0>(feasible_solutions[i]);
            }
        }
        tie(cost_incre, bundle_id, rider_type, newcost, newsource, newdest) = feasible_solutions[minindex];

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




void repair_little_random_old(
    vint & insert_order,
    vint & ids_to_build, 
    Solution & solution, 
    int * orders_ptr,
    Rider_Info & rider_info,
    map<string, int> & riders_available, 
    int * dist_mat_ptr,
    int l,
    bool use_power,
    default_random_engine & rng,
    int consider_size
    )
{
    int n_iter = sz(ids_to_build);
    int maxcapacity = max(max(rider_info.walk_info[0], rider_info.bike_info[0]), rider_info.car_info[0]);
    for(int rep=0; rep<n_iter; rep++){
        int order_id_to_append = ids_to_build[insert_order[rep]];
    
        vector<tuple<float, int, string, float, vint, vint>> feasible_solutions(sz(solution.solutions)+1, {INF, -1, "", 0.0, {}, {}});
        
        int n = sz(solution.solutions);
        vector<pair<float,int>> vv(n, {INF, -1});
        
        # pragma omp parallel for
        for (int bundle_id = 0; bundle_id < n; bundle_id++) {
            string rider_type = solution.get_rider_type(bundle_id);
            int tmp = 0;
            bool ok =true;
            int filled =0 ;
            for (auto x : solution.get_source(bundle_id))
            {
                tmp += promising[rider_type][order_id_to_append][x];
                filled += orders_ptr[f(x,2,3)];
                if(orders_ptr[f(order_id_to_append,1,3)] < orders_ptr[f(x,0,3)]) {ok = false;break;}
                if(orders_ptr[f(x,1,3)] < orders_ptr[f(order_id_to_append,0,3)]) {ok = false;break;}
            }
            if (filled + orders_ptr[f(order_id_to_append, 2,3)] > maxcapacity)ok=false;
            if(ok)vv[bundle_id] = {(float)tmp / sz(solution.get_source(bundle_id)), bundle_id};
        }
        
        vint bundle_indices;
        
        if (use_power) {
            vector<pair<float,int>> vvv;
            REP0(i, n) if (vv[i].se != -1) vvv.pb(vv[i]);
            sort(all(vvv));

            int n_consideration = min(max(consider_size, l/100), sz(vvv));
            int nn = min(sz(vvv), 5*n_consideration);

            vint candidate(nn); REP0(i, nn) candidate[i] = vvv[i].se;

            bundle_indices.resize(n_consideration);

            std::uniform_real_distribution<> dis(0.0, 1.0);
            for (int i = 0; i < n_consideration; ++i) {
                int noise_idx = floor(pow(dis(rng), 5) * sz(candidate));
                int bundle_id = candidate[noise_idx];
                bundle_indices[i] = bundle_id;
                candidate.erase(candidate.begin() + noise_idx);
            }
        }
        else {
            vector<pair<float,int>> candidates;
            for (auto pair : vv) {
                if (pair.se == -1) continue;
                candidates.pb(pair);
            }
            sort(candidates.begin(), candidates.end());

            int n_consideration = min(max(consider_size, l/100), sz(candidates));

            // 4 core optimization for 750 : problem 6 and 7
            n_consideration = min((n_consideration + 3) / 4 * 4, sz(candidates));

            bundle_indices.resize(n_consideration);
            REP0(i,n_consideration)bundle_indices[i] = candidates[i].second;
        }

        // Case 1. to bundle
        # pragma omp parallel for
        for(int bundle_id : bundle_indices) 
        {
            float cost_before = solution.get_cost(bundle_id);

            Res res = investigate_old(order_id_to_append, solution.get_source(bundle_id), solution.get_dest(bundle_id), rider_info, orders_ptr, dist_mat_ptr, l);
            
            for(string rider_type : res.optimal_order) {
                if(riders_available[rider_type] > 0 || solution.get_rider_type(bundle_id) == rider_type) {
                    if (res.get_feasibility(rider_type)) {
                        float cost_after = res.get_cost(rider_type);
                        float cost_incremental = cost_after - cost_before;

                        feasible_solutions[bundle_id] = make_tuple(cost_incremental, bundle_id, rider_type, cost_after, res.get_source(rider_type), res.get_dest(rider_type));
                        break;
                    }
                }
            }
        }
            
        // Case 2. to empty
        Res res = investigate_old(order_id_to_append, {}, {}, rider_info, orders_ptr, dist_mat_ptr, l);
        
        for (string rider_type: res.optimal_order) {
            if (riders_available[rider_type] >0 && res.get_feasibility(rider_type)) {
                float cost_incremental = res.get_cost(rider_type);                          
                feasible_solutions[sz(solution.solutions)] = make_tuple(cost_incremental, -1, rider_type,cost_incremental, res.get_source(rider_type), res.get_dest(rider_type));
                break;
            }
        }

        // remove it from ids_to_build
        auto it = std::find(all(ids_to_build), order_id_to_append); ids_to_build.erase(it);

        float cost_incre; int bundle_id; string rider_type; float newcost; vint newsource; vint newdest;

        int minindex = 0; float minval = get<0>(feasible_solutions[0]);
        for(int i=1 ; i<sz(feasible_solutions); i++)
        {
            if ( get<0>(feasible_solutions[i]) < minval)
            {
                minindex = i;
                minval = get<0>(feasible_solutions[i]);
            }
        }
        tie(cost_incre, bundle_id, rider_type, newcost, newsource, newdest) = feasible_solutions[minindex];

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