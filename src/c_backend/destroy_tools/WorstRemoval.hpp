#ifndef WORSTREMOVAL_HPP
#define WORSTREMOVAL_HPP

#include "../helper_classes.hpp"
#include "../helper_functions.hpp"
#include "../annealer.hpp"

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

vint update_bundle(int bundle_id, vector<pair<int,float>> & cost_of_order, Solution & sol, int K, Rider_Info & rider_info, int * dist_mat_ptr, int l )
{
    float org_cost = sol.get_cost(bundle_id);
    vint source = sol.get_source(bundle_id);
    vint dest = sol.get_dest(bundle_id);
    vint rider = rider_info.prepare(sol.get_rider_type(bundle_id)).se;
    
    int m = sz(source);
    if (m == 1) {
        int order = source[0];
        cost_of_order[order] = {round(100 * org_cost), 0.0};
        return source;
    }
    
    // source update
    for (int i = 0; i < m; ++i) {
        int order = source[i];
        int next = (i == m-1) ? dest[0]+K : source[i+1];
        int dist_increment = -dist_mat_ptr[f(order,next,l)];
        if (i != 0) dist_increment += dist_mat_ptr[f(source[i-1],next,l)] - dist_mat_ptr[f(source[i-1],order,l)];   
        cost_of_order[order].fi = -dist_increment * rider[2];
    }
    
    // dest update
    for (int i = 0; i < m; ++i) {
        int order = dest[i];
        int befo = (i == 0) ? source[m-1] : dest[i]+K;
        int dist_increment = -dist_mat_ptr[f(befo,order+K,l)];
        if (i != m-1) dist_increment += dist_mat_ptr[f(befo,dest[i+1]+K,l)] - dist_mat_ptr[f(order+K,dest[i+1]+K,l)];
        cost_of_order[order].fi -= dist_increment * rider[2];
    }
    
    // exception
    if (source[m-1] == dest[0]) {
        int order = dest[0];
        int befo = source[m-2], next = dest[1]+K;                
        int dist_increment = dist_mat_ptr[f(befo,next,l)] - dist_mat_ptr[f(befo,order,l)] - dist_mat_ptr[f(order,order+K,l)] - dist_mat_ptr[f(order+K,next,l)];     
        
        cost_of_order[order].fi = - dist_increment * rider[2] ;
    }
    
    for (int order : source) cost_of_order[order].se = org_cost - (float)cost_of_order[order].fi / 100;
    return source;
}

std::vector<std::pair<int, int>> removeElements(std::vector<std::pair<int, int>>& vv, const std::vector<int>& eraseplease) {
    // Convert toerase to an unordered_set for O(1) lookups
    std::unordered_set<int> toerase_set(eraseplease.begin(), eraseplease.end());

    // Use the erase-remove idiom with remove_if and check the unordered_set
    vv.erase(
        std::remove_if(vv.begin(), vv.end(),
            [&toerase_set](const std::pair<int, int>& p) {
                // Check if p.second is in the set
                return toerase_set.find(p.second) != toerase_set.end();
            }),
        vv.end());

    return vv; // Return the modified vv
}

std::vector<std::pair<int, int>> merge(std::vector<std::pair<int, int>>& vv, std::vector<std::pair<int, int>>& tmp) {
    std::vector<std::pair<int, int>> result;
    result.reserve(vv.size() + tmp.size());  // Reserve space for the result

    // Indices for iterating through vv and tmp
    size_t i = 0, j = 0;

    // Merge the two vectors while maintaining the decreasing order
    while (i < vv.size() && j < tmp.size()) {
        if (vv[i] > tmp[j]) {
            result.push_back(vv[i]);
            ++i;
        } else {
            result.push_back(tmp[j]);
            ++j;
        }
    }

    // Add remaining elements from vv, if any
    while (i < vv.size()) {
        result.push_back(vv[i]);
        ++i;
    }

    // Add remaining elements from tmp, if any
    while (j < tmp.size()) {
        result.push_back(tmp[j]);
        ++j;
    }

    return result;
}

vector<pair<int,int>> merge(vector<pair<int,int>> & vv, vint & considerplease, vector<pair<int,float>> & cost_of_order)
{
    vector<pair<int, int>> tmp;
    for(auto x : considerplease)if(cost_of_order[x].fi > 0)tmp.pb({cost_of_order[x].fi, x});
    sort(tmp.begin(), tmp.end(), greater<pair<int,int>>());

    return merge(vv, tmp);
}



vint worst_removal(int n_destroy, Solution & sol, map<string, int> & riders_available, vector<float> & samples, float worst_noise, int K, Rider_Info & rider_info, int * dist_mat_ptr, int l)
{
    vint to_remove;
            
    vector<pair<int,float>> cost_of_order(K, {0, 0.0});

    # pragma omp parallel for
    for (int bundle_id = 0; bundle_id < sz(sol.solutions); bundle_id++)update_bundle(bundle_id, cost_of_order, sol, K, rider_info, dist_mat_ptr, l);
    
    
    bool need_update = false;
    vint eraseplease;

    vector<pair<int,int>> vv;
    for (int i = 0; i < K; i++)if (cost_of_order[i].fi > 0) vv.pb({cost_of_order[i].fi, i});
    sort(vv.begin(), vv.end(), greater<pair<int,int>>());


    for (int rep = 0; rep < n_destroy; rep++)
    {
        if (need_update)vv  = merge(vv, eraseplease, cost_of_order);
        
        int noise_idx = floor(pow(samples[rep], worst_noise) * sz(vv));
        int worst_order = vv[noise_idx].se;
        to_remove.pb(worst_order);

        int worst_bundle_id = sol.get_bundle_id(worst_order);
        cost_of_order[worst_order].fi = -1;

        // update sol and riders_available.
        pair<bool, string> res = sol.remove_order(worst_order, worst_bundle_id, cost_of_order[worst_order].se);
        if(res.first==true)
        {
            riders_available[res.second]++;
            vv.erase(
            std::remove_if(vv.begin(), vv.end(),
                [worst_order](const std::pair<int, int>& p) {
                    return p.second == worst_order;
                }),
            vv.end());
            need_update = false;
        }
        else 
        {
            eraseplease = update_bundle(worst_bundle_id, cost_of_order, sol, K, rider_info, dist_mat_ptr, l);
            eraseplease.pb(worst_order);
            removeElements(vv, eraseplease);
            need_update = true;
        }

    }
    return to_remove;
}


#endif