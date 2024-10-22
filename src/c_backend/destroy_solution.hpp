#ifndef DESTROY_SOLUTION_HPP
#define DESTROY_SOLUTION_HPP

#include "helper_classes.hpp"
#include "helper_functions.hpp"
#include "./destroy_tools/methods.hpp"
#include "./destroy_tools/WorstRemoval.hpp"
#include "./destroy_tools/HARemoval.hpp"
#include "./destroy_tools/semiWorstRemoval.hpp"


#include <vector>
#include <bits/stdc++.h>
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


struct AdaptiveDestroyer{
    int n_sampling_method = 9; // method of different sampling methods

    Hyperparameter hparam;
    float rho = hparam.rho; //updating parameter
    float smoothing_ratio = hparam.smoothing_ratio; // smoothing for insertion operators
    float min_prob = hparam.min_prob; //destroy prob min
    float max_prob = hparam.max_prob; // destroy prob max
    float route_removal_ratio = hparam.route_removal_ratio; // route removal ratio
    int update_period = hparam.update_period; // period for updating prob
    float worst_noise = hparam.worst_noise;
    bool use_worst = hparam.use_worst;

    vint num_called; //number selected for each method
    vint num_success; // success number for each method

    int current_method; 
    int iter = 0;
    default_random_engine rng;
    int K; // number of orders
    int l; // 2 * number of orders. this is needed because we have 1D array
    vector<float> probability_distribution; // probability to choose each method

    int * orders_ptr; //orders array 1D
    Rider_Info rider_info; // info about each riders
    int * dist_mat_ptr; // distance matrix array 1D
    
    vector<vector<float>> dist_rel, time_rel, load_rel, shop_dist_rel, dlv_dist_rel;
    vector<vector<float>> history;
    
    AdaptiveDestroyer(
                    int K, 
                    int * orders_ptr,
                    Rider_Info & rider_info,
                    int * dist_mat_ptr_input,
                    Hyperparameter & hparam,
                    int seed=42
                    ):hparam(hparam), K(K), orders_ptr(orders_ptr), rider_info(rider_info), dist_mat_ptr(dist_mat_ptr_input)
    {     
        // Hyperparameter setting
        rho = hparam.rho; //updating parameter
        smoothing_ratio = hparam.smoothing_ratio; // smoothing for insertion operators
        min_prob = hparam.min_prob; //destroy prob min
        max_prob = hparam.max_prob; // destroy prob max
        route_removal_ratio = hparam.route_removal_ratio; // route removal ratio
        update_period = hparam.update_period; // period for updating prob
        worst_noise = hparam.worst_noise;
        use_worst = hparam.use_worst;
        
        if (!use_worst) n_sampling_method--;
        
        rng.seed(seed);
        l = 2*K;
        num_called.resize(n_sampling_method);
        num_success.resize(n_sampling_method);
        
        // dist_rel : a distance relativity measure
        dist_rel.resize(K, vector<float>(K));
        float dmax = 0; REP0(i, l) REP0(j, l) if (i != j) dmax = max(dmax, (float)dist_mat_ptr[f(i,j,l)]);
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                float dmean = (float)(dist_mat_ptr[f(i,j,l)] + dist_mat_ptr[f(i,j+K,l)] + dist_mat_ptr[f(i+K,j,l)] + dist_mat_ptr[f(i+K,j+K,l)]) / 4;
                dist_rel[i][j] = dmean / dmax;
            }
        }
        shop_dist_rel.resize(K, vector<float>(K));
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                shop_dist_rel[i][j] = (dist_mat_ptr[f(i,j,l)]);
            }
        }
        dlv_dist_rel.resize(K, vector<float>(K));
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                dlv_dist_rel[i][j] = (dist_mat_ptr[f(i+K,j+K,l)]);
            }
        }
        
        vector<float> window_center(K);
        REP0(i, K) window_center[i] = (float)(orders_ptr[f(i,0,3)] + orders_ptr[f(i,1,3)]) / 2;
        time_rel.resize(K, vector<float>(K));
        float tmax = 0;
        REP0(i, K) REP0(j, K) { time_rel[i][j] = abs(window_center[i] - window_center[j]); tmax = max(tmax, time_rel[i][j]); }
        REP0(i, K) REP0(j, K) { time_rel[i][j] /= tmax; }
        
        load_rel.resize(K, vector<float>(K));
        float lmax = 0;
        REP0(i, K) REP0(j, K) { load_rel[i][j] = abs(orders_ptr[f(i,2,3)] - orders_ptr[f(j,2,3)]); lmax = max(lmax, load_rel[i][j]); }
        REP0(i, K) REP0(j, K) { load_rel[i][j] /= lmax; }

        history.resize(l, vector<float>(l, INF));
        
        // make probability_distribution
        probability_distribution.resize(n_sampling_method);
        for(int i=0; i<n_sampling_method; i++)probability_distribution[i] = 1.0/(float)n_sampling_method;
    }

    /* ########## preparation functions ################*/

    int get_number_to_destroy()
    {
        std::uniform_real_distribution<> dis(min_prob, max_prob);
        float prob = dis(rng);
        int n_destroy = floor(prob * K);
        return min(hparam.max_destroy, n_destroy);
    }

    int get_sampling_method()
    {
        std::discrete_distribution<int> distribution(probability_distribution.begin(), probability_distribution.end());
        return distribution(rng);
    }

    vector<float> get_samples(int n_destroy)
    {
        vector<float> samples;
        std::uniform_real_distribution<> realunidist(0.0, 1.0);
        REP0(rep, n_destroy){samples.pb(realunidist(rng));}
        return samples;
    }

    

    /* ########## destroy function ################*/
    void destroy_given_ids(vint & ids_to_destroy, Solution & sol, map<string, int> & riders_available)
    {   
        map<int, int> id_to_bundle;
        for (int bundle_id = 0; bundle_id < sz(sol.solutions); bundle_id++) {
            for (int id : sol.get_source(bundle_id))
                id_to_bundle[id] = bundle_id;
        }

        vint deleteid;
        for(auto id : ids_to_destroy){
            int bundleid = id_to_bundle[id];
            Bundle bundle = sol.solutions[bundleid];
            
            vint new_source_order = get<2>(bundle);
            new_source_order.erase(find(all(new_source_order), id));
            vint new_dest_order = get<3>(bundle);
            new_dest_order.erase(find(all(new_dest_order), id));

            string ridertype = get<0>(bundle);
            pair<int *, vint> info = rider_info.prepare(ridertype);
            float new_cost =  get_bundle_cost(new_source_order, new_dest_order, info.second, dist_mat_ptr,l);

            if (new_source_order.empty()) {
                riders_available[ridertype]++;
                deleteid.pb(bundleid);
            }
            else {
                Bundle solinstance = make_tuple(ridertype, new_cost, new_source_order, new_dest_order);
                sol.solutions[bundleid] = solinstance;
            }
        }
        sol.remove(deleteid);
    }

    /* ########## store results for updating probabilities ################*/

    void update_prob()
    {
        vector<float> new_prob(n_sampling_method);
        REP0(i, n_sampling_method) if (num_called[i] > 0) new_prob[i] = (float)num_success[i] / num_called[i];
        float sum = std::accumulate(all(new_prob), 0.0f);
        REP0(i, n_sampling_method) if (sum > 0) new_prob[i] = new_prob[i] / sum;
        
        for (int i = 0; i < n_sampling_method; ++i) {
            probability_distribution[i] = (1 - rho) * probability_distribution[i] + rho * new_prob[i];
        }
        // Normalize the probability distribution
        sum = std::accumulate(probability_distribution.begin(), probability_distribution.end(), 0.0f);
        for (int i = 0; i < n_sampling_method; ++i) {
            probability_distribution[i] = smoothing_ratio + (1 - smoothing_ratio * n_sampling_method) * probability_distribution[i] / sum;
        }
        // Reset num_called and num_success
        std::fill(num_called.begin(), num_called.end(), 0);
        std::fill(num_success.begin(), num_success.end(), 0);
    }

    void update(int update_score)
    {
        num_called[current_method] ++;
        num_success[current_method] += update_score;
        if(iter % update_period == 0) update_prob();
    }

    /* ############################################### Main function ############################################### */
    vint destroy(Solution & sol, map<string, int> & riders_available){
        iter++; 
        current_method = get_sampling_method();
        int n_destroy = get_number_to_destroy();

        history_update(sol,K, history);
        
        vint ids_destroyed;
        if (current_method == 0)ids_destroyed = random_removal(n_destroy, K, rng);
        else if (current_method ==1)ids_destroyed =  distance_oriented_removal(n_destroy, K, rng, dist_rel);
        else if (current_method ==2)ids_destroyed = route_removal(n_destroy, sol, route_removal_ratio, rng);
        else if (current_method ==3)ids_destroyed = shaw_removal(n_destroy, sol, K, hparam.shaw_d, hparam.shaw_t, hparam.shaw_l, hparam.shaw_noise, dist_rel, time_rel, load_rel, rng);
        else if (current_method ==4)ids_destroyed = historical_action_pair_removal(n_destroy, sol, K, history);
        else if(current_method ==5)ids_destroyed=shop_oriented_removal(n_destroy, K, rng, shop_dist_rel);
        else if(current_method ==6)ids_destroyed=dlv_oriented_removal(n_destroy, K, rng, dlv_dist_rel);
        else if (current_method ==7)
        {
            vector<float> samples = get_samples(n_destroy);
            ids_destroyed = semiworst_removal(n_destroy, sol, samples, worst_noise,K, rider_info, dist_mat_ptr, l);
        }
        else if (current_method ==8)
        {
            vector<float> samples = get_samples(n_destroy);
            ids_destroyed = worst_removal(n_destroy, sol, riders_available, samples, worst_noise,K, rider_info, dist_mat_ptr, l);
        }
        if(current_method!=8)destroy_given_ids(ids_destroyed, sol, riders_available);
        return ids_destroyed;
    }
};

#endif