
#ifndef HAREMOVAL_HPP
#define HAREMOVAL_HPP

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

vint historical_action_pair_removal(int n_destroy, Solution & sol, int K, vector<vector<float>>& history)
{
    vector<pair<float,int>> vv(K);
    REP0(i, K) vv[i].se = i;
    
    # pragma omp parallel for
    for (int bundle_id = 0; bundle_id < sz(sol.solutions); bundle_id++) {
        const vint& s = sol.get_source(bundle_id);
        const vint& d = sol.get_dest(bundle_id);
        
        int n = sz(s);
        for (int i = 0; i < n-1; ++i) {
            vv[s[i]].fi += history[s[i]][s[i+1]];
            vv[s[i+1]].fi += history[s[i]][s[i+1]];
        }
        vv[s[n-1]].fi += history[s[n-1]][d[0]+K];
        vv[d[0]].fi += history[s[n-1]][d[0]+K];
        for (int i = 0; i < n-1; ++i) {
            vv[d[i]].fi += history[d[i]+K][d[i+1]+K];
            vv[d[i+1]].fi += history[d[i]+K][d[i+1]+K];
        }
    }
    sort(all(vv), greater<pair<float,int>>());
    
    vint to_remove;
    for (int i = 0; i < n_destroy; ++i) to_remove.pb(vv[i].se);
    return to_remove;
}


#endif