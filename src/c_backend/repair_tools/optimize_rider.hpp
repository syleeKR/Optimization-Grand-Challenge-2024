#ifndef OPTIMIZE_RIDER_HPP
#define OPTIMIZE_RIDER_HPP

#include "../helper_classes.hpp"
#include "../helper_functions.hpp"
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



void optimize_rider_type(vint & indices, Solution & sol, map<string, int> & riders_available, Rider_Info & rider_info, int * orders_ptr, int * dist_mat_ptr, int l)
{
    // always three things to update. sol.solutions, sol.cost, riders_available
    for(auto i : indices)
    {
        string rider_type = sol.get_rider_type(i);
        vint source = sol.get_source(i);
        vint dest = sol.get_dest(i);

        vector<pair<string, float>> res = get_optimal_rider_type(source, dest, rider_info, orders_ptr, dist_mat_ptr, l);

        for(auto x: res)
        {
            if (x.first == rider_type)break;
            if(riders_available[x.first] >0){
                string new_rider_type = x.first;
                float newcost = x.second;
                sol.solutions[i] = make_tuple(new_rider_type, newcost, source, dest);
                riders_available[rider_type]++;
                riders_available[new_rider_type] --;
                break;
            }
        }
    }
    sol.update_cost();
}

#endif