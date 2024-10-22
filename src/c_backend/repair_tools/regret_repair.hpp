#ifndef REGRET_REPAIR_HPP
#define REGRET_REPAIR_HPP

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

map<int, vector<tuple<float, int, string>>> get_sorted_feasible_solutions(vector<vector<tuple<float, int, string>>> & feasible_solutions)
{
    map<int, vector<tuple<float, int, string>>> temp;
    for (int i=0; i<sz(feasible_solutions); i++)
    {
        vector<tuple<float, int, string>> x = feasible_solutions[i];
        if (sz(x)>0){
            for(auto tup : x)temp[i].push_back(tup);
        }
    }
    for (auto& x : temp) {
        // Sort the vector based on the float element of the tuple
        std::sort(x.second.begin(), x.second.end(), [](const std::tuple<float, int, std::string>& a, const std::tuple<float, int, std::string>& b) {
            return std::get<0>(a) < std::get<0>(b);
        });
    }
    return temp;
}

vector<float> probability_distribution = {0.8, 0.15, 0.05};
std::discrete_distribution<int> distribution(probability_distribution.begin(), probability_distribution.end());
std::uniform_real_distribution<> real_distribution(0.29999, 0.3);

int get_order_id_to_append(default_random_engine &rng, map<int, vector<tuple<float, int, string>>> & feasible_solutions, int regret_option = 1)
{
    int how_fool = distribution(rng);
    
    int order_id_to_append = -1;
    if (regret_option>0)
    {
    /* ######################### Use k-regret ############################### */
        // Calculate regret values
        vector<pair<float, int>> vv;
        for (const auto& pair : feasible_solutions) {
            int key = pair.first; float regret = 0;
            const vector<tuple<float, int, string>>& vec = pair.second;
            for(int i = 1; i<=min(regret_option, sz(vec)-1); i++) regret += (std::get<0>(vec[i]) - std::get<0>(vec[0]));
            for(int i=sz(vec); i<=regret_option; i++) regret += INT_MAX; //add infinite value. we prefer ones with little possible options left
            
            vv.push_back({regret, key});
        }

        sort(all(vv));
        int n = max(0, sz(vv)-1 - how_fool);
        order_id_to_append = vv[n].second;
    }
    else if(regret_option==0)
    {
        /*################Use Greedy#########################*/
        vector<pair<float, int>> vv;
        for(const auto & pair : feasible_solutions)
        {
            int i = pair.first;
            const vector<tuple<float, int, string>>& vec = pair.second;
            if (!vec.empty()) {
                vv.push_back({get<0>(vec[0]), i});
            }
            
            sort(all(vv));
            int n = min(sz(vv)-1, how_fool);
            order_id_to_append = vv[n].second;
        }
    }
    return order_id_to_append;
}


void repair_regret(
    default_random_engine &rng,
    vint & ids_to_build, 
    Solution & solution, 
    int * orders_ptr,
    Rider_Info & rider_info,
    map<string, int> & riders_available, 
    int * dist_mat_ptr,
    int l,
    int regret_option=1)
{    
    Hyperparameter hparam;

    //// 1. Let's go
    int K = l/2;
    Cache cache = Cache(); Res res;
    while ( !ids_to_build.empty() )
    {
        vector<vector<tuple<float, int, string>>> feasible_solutions(K);
        vector<vector<pair<vint, Res>>> to_append_to_cache(K);

        # pragma omp parallel for
        for (int order_id : ids_to_build)
        {
            Res res;
            
            for(int bundle_id = 0; bundle_id < sz(solution.solutions); bundle_id++)
            {
                vint merged_id_list = solution.get_dest(bundle_id);
                float cost_before = solution.get_cost(bundle_id);
                merged_id_list.push_back(order_id);

                if (cache.check(merged_id_list)) {
                    res = cache.retrieve(merged_id_list);
                } 
                else {
                    res = investigate(order_id, solution.get_source(bundle_id),solution.get_dest(bundle_id),rider_info, orders_ptr, dist_mat_ptr, l);
                    to_append_to_cache[order_id].push_back({merged_id_list, res});
                }
                
                for(string rider_type : res.optimal_order) {
                    if(riders_available[rider_type] > 0 || solution.get_rider_type(bundle_id) == rider_type) {
                        if (res.get_feasibility(rider_type)) {
                            float cost_after = res.get_cost(rider_type);
                            float cost_incremental = cost_after - cost_before;
                            
                            float capa = rider_info.prepare(rider_type).second[0];
                            for (int id : res.get_source(rider_type)) capa -= orders_ptr[f(id,2,3)];
                            
                            float alpha2 = real_distribution(rng);
                            float C = (1.0 - alpha2) * cost_incremental + alpha2 * capa;
                            
                            feasible_solutions[order_id].push_back(make_tuple(C, bundle_id, rider_type));
                            break;
                        }
                    }
                }
            }
            
            vint one_element_vector = {order_id};
            if (cache.check(one_element_vector)) res = cache.retrieve(one_element_vector);
            else
            {
                vint empty_source, empty_dest;
                res = investigate(order_id, empty_source,empty_dest,rider_info, orders_ptr, dist_mat_ptr, l);
                to_append_to_cache[order_id].push_back({one_element_vector, res});
            }

            for (string rider_type: res.optimal_order)
            {
                if (riders_available[rider_type] >0 && res.get_feasibility(rider_type))
                {
                    float cost_incremental = res.get_cost(rider_type);
                    
                    float capa = rider_info.prepare(rider_type).second[0];
                    for (int id : res.get_source(rider_type)) capa -= orders_ptr[f(id,2,3)];
                    
                    float C = hparam.alpha1 * cost_incremental + hparam.alpha2 * capa;
                    
                    feasible_solutions[order_id].push_back(make_tuple(C, -1, rider_type));
                    break;
                }
            }
        }
        
        // to_append_to_cache to cache
        for(int i = 0; i < K; i++) {
            vector<pair<vint, Res>> vector_of_values_to_append = to_append_to_cache[i];
            if (sz(vector_of_values_to_append) > 0){
                for (auto v_res_pair : vector_of_values_to_append) {
                    cache.append(v_res_pair.first, v_res_pair.second);
                }
            }
        }
        // get best order_id
        map<int, vector<tuple<float, int, string>>> sorted_feasible_solutions = get_sorted_feasible_solutions(feasible_solutions);
        int order_id_to_append = get_order_id_to_append(rng, sorted_feasible_solutions, regret_option);
        if (order_id_to_append == -1) break;

        // remove it from ids_to_build
        auto it = std::find(all(ids_to_build), order_id_to_append);
        if (it != ids_to_build.end()) ids_to_build.erase(it);

        // TODO : what if not a single id can be appended??
        // get best inserted bundle
        float cost_incre; int bundle_id; string rider_type;
        tie(cost_incre, bundle_id, rider_type) = sorted_feasible_solutions[order_id_to_append][0];

        // remove original bundle
        vint merged_id_list;
        if (bundle_id != -1) {
            merged_id_list = solution.get_source(bundle_id);
            riders_available[solution.get_rider_type(bundle_id)]++;
            solution.remove(bundle_id);
        }
        merged_id_list.push_back(order_id_to_append);

        // append inserted bundle
        Res cached_res = cache.retrieve(merged_id_list);
        solution.append(make_tuple(rider_type, cached_res.get_cost(rider_type),cached_res.get_source(rider_type), cached_res.get_dest(rider_type)));
        riders_available[rider_type]--;
    }
}

#endif