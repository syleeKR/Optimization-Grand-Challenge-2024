#ifndef SEMIWORSTREMOVAL_HPP
#define SEMIWORSTREMOVAL_HPP

#include "../helper_classes.hpp"
#include "../helper_functions.hpp"
#include "../annealer.hpp"
#include "WorstRemoval.hpp"

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

vint update_bundle(int bundle_id, vector<int> & cost_of_order, Solution & sol, int K, Rider_Info & rider_info, int * dist_mat_ptr, int l )
{
    vint source = sol.get_source(bundle_id);
    vint dest = sol.get_dest(bundle_id);
    vint rider = rider_info.prepare(sol.get_rider_type(bundle_id)).se;
    
    int m = sz(source);
    if (m == 1) {
        // 여기가 ㅅㅂ 예측이 안됨 sol.cost가져다 쓰는거여서 완전 랜덤이었을텐데.
        int order = source[0];
        float one_source_cost = rider[1] + rider[2] * dist_mat_ptr[f(order, order+K, l)]/100.0;
        cost_of_order[order] = round(100 * one_source_cost);
        return source;
    }
    
    // source update
    for (int i = 0; i < m; ++i) {
        int order = source[i];
        int next = (i == m-1) ? dest[0]+K : source[i+1];
        int dist_increment = -dist_mat_ptr[f(order,next,l)];
        if (i != 0) dist_increment += dist_mat_ptr[f(source[i-1],next,l)] - dist_mat_ptr[f(source[i-1],order,l)];   
        cost_of_order[order] -= dist_increment * rider[2];
    }
    
    // dest update
    for (int i = 0; i < m; ++i) {
        int order = dest[i];
        int befo = (i == 0) ? source[m-1] : dest[i]+K;
        int dist_increment = -dist_mat_ptr[f(befo,order+K,l)];
        if (i != m-1) dist_increment += dist_mat_ptr[f(befo,dest[i+1]+K,l)] - dist_mat_ptr[f(order+K,dest[i+1]+K,l)];
        cost_of_order[order] -= dist_increment * rider[2];
    }
    
    // exception
    if (source[m-1] == dest[0]) {
        int order = dest[0];
        int befo = source[m-2], next = dest[1]+K;                
        int dist_increment = dist_mat_ptr[f(befo,next,l)] - dist_mat_ptr[f(befo,order,l)] - dist_mat_ptr[f(order,order+K,l)] - dist_mat_ptr[f(order+K,next,l)];     
        
        cost_of_order[order] = - dist_increment * rider[2] ;
    }
    
    return source;
}


vector<pair<int,int>> merge(vector<pair<int,int>> & vv, vint & considerplease, vector<int> & cost_of_order)
{
    vector<pair<int, int>> tmp;
    for(auto x : considerplease)if(cost_of_order[x] > 0)tmp.pb({cost_of_order[x], x});
    sort(tmp.begin(), tmp.end(), greater<pair<int,int>>());

    return merge(vv, tmp);
}


// ha sibal 그때 그시절이랑 똑같은진 모르겠다 일단 helper class쪽 오류 고친거 때문에 org_cost 쓰던 m=1일때 로직 좀 다른거 빼고는 같은듯한데
vint semiworst_removal(int n_destroy, Solution & originalsol, vector<float> & samples, float worst_noise, int K, Rider_Info & rider_info, int * dist_mat_ptr, int l)
{
    Solution sol = originalsol; 
    vint to_remove;
            
    vector<int> cost_of_order(K, 0);

    # pragma omp parallel for
    for (int bundle_id = 0; bundle_id < sz(sol.solutions); bundle_id++)update_bundle(bundle_id, cost_of_order, sol, K, rider_info, dist_mat_ptr, l);
    
    vector<pair<int,int>> vv;
    for (int i = 0; i < K; i++)if (cost_of_order[i] > 0) vv.pb({cost_of_order[i], i});
    sort(vv.begin(), vv.end(), greater<pair<int,int>>());
    
    vint eraseplease;
    int worst_order;
    for (int rep = 0; rep < n_destroy; rep++)
    {
        if(rep>=1)
        {
            eraseplease.pb(worst_order);
            removeElements(vv, eraseplease);
            merge(vv, eraseplease, cost_of_order);
        }

        int noise_idx = floor(pow(samples[rep], worst_noise) * sz(vv));
        worst_order = vv[noise_idx].se;
        
        int worst_bundle_id = sol.get_bundle_id(worst_order);
        cost_of_order[worst_order] = -1;
        
        to_remove.pb(worst_order);

        // update sol and riders_available.
        float dummy = 0.0;
        pair<bool, string> res = sol.remove_order(worst_order, worst_bundle_id, dummy);
        if(worst_bundle_id<sz(sol.solutions))eraseplease = update_bundle(worst_bundle_id, cost_of_order, sol, K, rider_info, dist_mat_ptr, l);
        else eraseplease.clear();
    }
    return to_remove;
}


#endif