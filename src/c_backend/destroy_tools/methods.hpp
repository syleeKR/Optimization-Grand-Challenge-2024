#ifndef METHODS_HPP
#define METHODS_HPP

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


vint random_removal(int n_destroy, int K, default_random_engine & rng){
        std::vector<int> indices(K);
        for (int i = 0; i < K; ++i) indices[i] = i;

        std::shuffle(indices.begin(), indices.end(), rng);
        std::vector<int> to_remove(indices.begin(), indices.begin() + n_destroy);
        return to_remove;
    }

vint distance_oriented_removal(int n_destroy, int K,default_random_engine & rng , vector<vector<float>> & dist_rel)
{
    // sample first order randomly
    std::uniform_int_distribution<std::size_t> dist(0, K-1);
    int id1 = dist(rng);
    
    vector<pair<float,int>> vv;
    REP0(id2, K) vv.push_back({dist_rel[id1][id2], id2});
    sort(all(vv));
    
    vint to_remove;
    REP0(i, n_destroy) to_remove.push_back(vv[i].se);
    return to_remove;
}
vint shop_oriented_removal(int n_destroy, int K, default_random_engine & rng, vector<vector<float>> & shop_dist_rel)
{
    // sample first order randomly
    std::uniform_int_distribution<std::size_t> dist(0, K-1);
    int id1 = dist(rng);
    
    vector<pair<float,int>> vv;
    REP0(id2, K) vv.push_back({shop_dist_rel[id1][id2], id2});
    sort(all(vv));
    
    vint to_remove;
    REP0(i, n_destroy) to_remove.push_back(vv[i].se);
    return to_remove;
}
vint dlv_oriented_removal(int n_destroy, int K, default_random_engine & rng, vector<vector<float>> & dlv_dist_rel)
{
    // sample first order randomly
    std::uniform_int_distribution<std::size_t> dist(0, K-1);
    int id1 = dist(rng);
    
    vector<pair<float,int>> vv;
    REP0(id2, K) vv.push_back({dlv_dist_rel[id1][id2], id2});
    sort(all(vv));
    
    vint to_remove;
    REP0(i, n_destroy) to_remove.push_back(vv[i].se);
    return to_remove;
}

vint route_removal(int n_destroy, Solution & sol, float route_removal_ratio, default_random_engine & rng)
{
    int n_bundle = sz(sol.solutions);
    std::vector<int> indices(n_bundle);
    for (int i = 0; i < n_bundle; ++i)indices[i] = i;
    std::shuffle(indices.begin(), indices.end(), rng);

    int n_destroy_bundle = floor(n_bundle * route_removal_ratio);
    vint to_remove;
    for (int i=0; i<n_destroy_bundle; ++i) {
        for (int id : sol.get_source(indices[i]))
            to_remove.push_back(id);
    }
    return to_remove;
}

vint shaw_removal(int n_destroy, Solution & sol, int K, float shaw_d , float shaw_t, float shaw_l, float shaw_noise, vector<vector<float>> & dist_rel, vector<vector<float>> & time_rel, vector<vector<float>> & load_rel, default_random_engine & rng)
{
    std::uniform_int_distribution<std::size_t> init_uniform(0, K-1);
    int init_order = init_uniform(rng);
    vint to_remove = { init_order };
    
    vector<bool> chk(K);
    chk[init_order] = true;
    
    std::uniform_real_distribution<> dis(0.0, 1.0);
    while (sz(to_remove) < n_destroy)
    {
        std::uniform_int_distribution<std::size_t> random_choose(0, sz(to_remove) - 1);
        int id1 = to_remove[ random_choose(rng) ];
        
        vector<pair<float,int>> vv;
        REP0(id2, K) if (!chk[id2]) {
            float rel = shaw_d * dist_rel[id1][id2] + shaw_t * time_rel[id1][id2] + shaw_l * load_rel[id1][id2];
            vv.push_back({rel, id2});
        }
        sort(all(vv));
        int noise_idx = floor(pow(dis(rng), shaw_noise) * sz(vv));
        int shaw_order = vv[noise_idx].se;
        
        to_remove.pb(shaw_order);
        chk[shaw_order] = true;
    }
    return to_remove;
}


void history_update(Solution & sol, int K, vector<vector<float>> & history)
{
    float cost = sol.cost / K;
    for (int bundle_id = 0; bundle_id < sz(sol.solutions); bundle_id++) {
        const vint& s = sol.get_source(bundle_id);
        const vint& d = sol.get_dest(bundle_id);
        
        int n = sz(s);
        for (int i = 0; i < n-1; ++i) history[s[i]][s[i+1]] = min(history[s[i]][s[i+1]], cost);
        history[s[n-1]][d[0]+K] = min(history[s[n-1]][d[0]+K], cost);
        for (int i = 0; i < n-1; ++i) history[d[i]+K][d[i+1]+K] = min(history[d[i]+K][d[i+1]+K], cost);
    }
}


#endif