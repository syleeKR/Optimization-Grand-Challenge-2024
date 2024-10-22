#ifndef POLISH_HPP
#define POLISH_HPP


#include <bits/stdc++.h>
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

/* ################################################################################################ */
#include "helper_classes.hpp" 
#include "utils.hpp"



void polish_bundle(tuple<string, float, vint, vint> & bundle,
                    int * T_ptr,
                    vint & rider_info,
                    int * orders_ptr,
                    int * dist_mat_ptr,
                    int l)
{
    vint source = get<2>(bundle);
    vint dest = get<2>(bundle);
    int n = sz(source);
    

    int permutation_count = 120;
    if (n==1)permutation_count=1;
    else if(n==2)permutation_count=2;
    else if(n==3)permutation_count=6;
    else if(n==4)permutation_count=24;

    vector<vector<int>> permutations_source;
    vector<vector<int>> permutations_dest;

    int current_count = 0;
    do {
        permutations_source.push_back(source);
        current_count++;
    } while (std::next_permutation(source.begin(), source.end()) && current_count < permutation_count);


    current_count = 0;
    do {
        permutations_dest.push_back(dest);
        current_count++;
    } while (std::next_permutation(dest.begin(), dest.end()) && current_count < permutation_count);
    
    

    vector< pair<pair<bool, float>, pair<vint, vint>>> storage;
    omp_set_num_threads(4);
    #pragma omp parallel for collapse(2)
    for(auto s: permutations_source)
    {
        for(auto d: permutations_dest)
        {
            pair<bool, float> res = check_path(s,d,T_ptr, rider_info, orders_ptr, dist_mat_ptr, l);

            #pragma omp critical
            {
                storage.emplace_back(std::make_pair(res, std::make_pair(s, d)));
            }
            
        }
    }
    float min_cost = get<1>(bundle);
    for(auto x: storage)
    {
        if (x.first.first ==true && x.first.second < min_cost)
        {
            min_cost = x.first.second;
            get<2>(bundle) = x.second.first;
            get<3>(bundle) = x.second.second;
            get<1>(bundle) = min_cost;
        }
    }
}

void polish(double lefttime, Solution & sol,
    int * orders_ptr,
    int * walk_T_ptr,
    vint & walk_info,
    int * bike_T_ptr,
    vint & bike_info,
    int * car_T_ptr,
    vint & car_info,
    int * dist_mat_ptr,
    int l)
{
    auto start_time = chrono::high_resolution_clock::now();

    int n = sz(sol.solutions);

    for(int i=0; i<n; i++)
    {
        if(get_time(start_time) > lefttime-0.5)break;
        if(get<0>(sol.solutions[i]) == "WALK")polish_bundle(sol.solutions[i], walk_T_ptr, walk_info, orders_ptr, dist_mat_ptr, l);
        else if(get<0>(sol.solutions[i]) == "BIKE")polish_bundle(sol.solutions[i], bike_T_ptr, bike_info, orders_ptr, dist_mat_ptr, l);
        else polish_bundle(sol.solutions[i], car_T_ptr, car_info, orders_ptr, dist_mat_ptr, l);
    }
    sol.update_cost();
}




    /* ##################################### Polishing final result ##################################### */

    // devote last 2 sec to polishing
    // local_search.search(bestsol, riders_available, 5);
    // double lefttime = hparam.polish_time;
    // polish(lefttime, bestsol,
    //     orders_ptr, 
    //     walk_T_ptr, walk_info,
    //     bike_T_ptr, bike_info,
    //     car_T_ptr, car_info,
    //     dist_mat_ptr,
    //     l);






# endif