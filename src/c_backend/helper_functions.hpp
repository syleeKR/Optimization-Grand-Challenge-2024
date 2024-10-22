
#ifndef HELPER_FUNCTIONS_HPP
#define HELPER_FUNCTIONS_HPP

#include "helper_classes.hpp"
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
pair<bool, float> check_path(const vector<int> & source, 
                             const vector<int> & dest, 
                             int * T_ptr, 
                             vector<int> & rider, 
                             int * orders_ptr, 
                             int * dist_mat_ptr,
                             int l);
// ################################################################################################### //

int f(int i, int j, int l){return i*l+j;}

float get_bundle_cost(vint & source, vint & dest, vint & rider, int * dist_mat_ptr, int l){
    // calculate the cost of a bundle
    int K= l/2;
    int n =sz(source);
    if (n==0)return 0;
    int dist = 0;
    REP(i,1,n-1){dist += dist_mat_ptr[f(source[i-1], source[i], l)];}
    dist += dist_mat_ptr[f(source[n-1], dest[0] + K, l)];
    REP(i,1,n-1){ dist += dist_mat_ptr[f(dest[i-1] + K, dest[i] + K, l)] ; }

    return rider[1] + rider[2] * (float)dist/100;
}

vector<pair<string, float>> get_optimal_rider_type(vint & source, vint & dest, Rider_Info & rider, int * orders_ptr, int * dist_mat_ptr, int l)
{
    // given source, dest , find all feasible rider types. Return the feasible types with their costs, sorted.

    int total_volume = 0;
    for (int id : source)total_volume += orders_ptr[f(id,2,3)];
    pair<int * ,vint> walk = rider.prepare("WALK");
    pair<int * ,vint> bike = rider.prepare("BIKE");
    pair<int * ,vint> car = rider.prepare("CAR");

    pair<bool, float> res_walk = check_path(source, dest, walk.first, walk.second, orders_ptr, dist_mat_ptr, l);
    pair<bool, float> res_bike = check_path(source, dest, bike.first, bike.second, orders_ptr, dist_mat_ptr, l);
    pair<bool, float> res_car = check_path(source, dest, car.first, car.second, orders_ptr, dist_mat_ptr, l);

    vector<pair<float, string>> feasible_types;
    if (res_walk.first && total_volume <= walk.second[0])feasible_types.pb({res_walk.second, "WALK"});
    if (res_bike.first && total_volume <= bike.second[0])feasible_types.pb({res_bike.second, "BIKE"});
    if (res_car.first && total_volume <= car.second[0])feasible_types.pb({res_car.second, "CAR"});

    sort(all(feasible_types));

    vector<pair<string, float>> sorted_rider_type;
    for(auto x : feasible_types)sorted_rider_type.pb({x.se, x.fi});

    return sorted_rider_type;
}

/* ######################################################################################################### */
pair<bool, float> check_path(const vector<int> & source, 
                             const vector<int> & dest, 
                             int * T_ptr, 
                             vector<int> & rider, 
                             int * orders_ptr, 
                             int * dist_mat_ptr,
                             int l)
{
    int K = l/2, n = sz(source);
    
    int dist = 0, t_curr = orders_ptr[f(source[0],0,3)];
    for(int i = 1; i < n; i++)
    {
        dist += dist_mat_ptr[f(source[i-1], source[i], l)];
        t_curr =  max(t_curr + T_ptr[f(source[i-1] , source[i], l)], orders_ptr[f(source[i],0,3)]);
    }

    int bef = source[n-1];
    for (int i = 0; i < n; i++) {
        int aft = dest[i]+K;
        dist += dist_mat_ptr[f(bef, aft, l)];
        t_curr += T_ptr[f(bef, aft, l)];
        bef = aft;
        
        if (t_curr > orders_ptr[f(dest[i],1,3)]) return {false, INF};
    }
    return {true, rider[1] + rider[2] * (float)dist / 100};
}
    
/* ######################################################################################################### */

vector<vector<int>> insertAtAllPositions(int x, const vint & v) {
    int n = v.size();
    std::vector<std::vector<int>> res;

    for (int i = 0; i <= n; ++i) {
        std::vector<int> temp;
        for (int j = 0; j < i; ++j) temp.push_back(v[j]);
        temp.push_back(x);
        for (int j = i; j < n; ++j) temp.push_back(v[j]);
        res.push_back(temp);
    }
    return res;
}

tuple<bool, float, vector<int>, vector<int>> find_best_path_old(int x,
                                                            const vector<int> & curr_source,
                                                            const vector<int> & curr_dest, 
                                                            int * T_ptr, 
                                                            vector<int> & rider, 
                                                            int * orders_ptr, // n by 3 ndarray
                                                            int * dist_mat_ptr, 
                                                            int l)
{
    // 1. volume constraint
    int total_volume = orders_ptr[f(x,2,3)];
    for (int id : curr_source) total_volume += orders_ptr[f(id,2,3)];
    if (total_volume > rider[0]) return make_tuple(false, INT_MAX, std::vector<int>{}, std::vector<int>{});

    // 2. calculate dest deadline
    int K = l/2, n = sz(curr_source);
    vint deadline(n);
    for (int i = n-1; i >= 0; --i) {
        deadline[i] = orders_ptr[f(curr_dest[i], 1, 3)];
        if (i != n-1) deadline[i] = min(deadline[i], deadline[i+1] - T_ptr[f(curr_dest[i]+K, curr_dest[i+1]+K, l)]);
    }

    // 3. Now we check the timing constraint
    bool feasible = false;
    float best_cost = INF;
    int best_source_idx, best_dest_idx;

    vector<vint> source_sampled = insertAtAllPositions(x, curr_source);   
    vector<vint> dest_sampled = insertAtAllPositions(x, curr_dest);
    
    const vint& d = curr_dest;
    for (int kk = 0; kk <= n; ++kk)
    {
        const vint& s = source_sampled[kk];
        
        // calculate dist and source time
        int dist = 0, time = orders_ptr[f(s[0],0,3)];
        for (int i = 1; i <= n; ++i) {
            dist += dist_mat_ptr[f(s[i-1], s[i], l)];
            time = max(time + T_ptr[f(s[i-1], s[i], l)], orders_ptr[f(s[i],0,3)]);
        }
        int bef = s[n];
        for (int i = 0; i < n; ++i) {
            dist += dist_mat_ptr[f(bef, d[i]+K, l)];
            bef = d[i]+K;
        }

        // insert and calculate new dest
        bef = s[n];
        int deadline_x = orders_ptr[f(x, 1, 3)];
        for (int i = 0; i < n; ++i) {
            if (i != 0 && time > orders_ptr[f(bef-K, 1, 3)]) break;
            int aft = d[i]+K;

            bool chk1 = (time + T_ptr[f(bef, x+K, l)] <= deadline_x);
            bool chk2 = (time + T_ptr[f(bef, x+K, l)] + T_ptr[f(x+K, aft, l)] <= deadline[i]);
            if (chk1 && chk2) {
                int new_dist = dist + dist_mat_ptr[f(bef, x+K, l)] + dist_mat_ptr[f(x+K, aft, l)] - dist_mat_ptr[f(bef, aft, l)];
                float new_cost = rider[1] + rider[2] * (float)new_dist / 100;

                if (new_cost < best_cost) {
                    best_cost = new_cost;
                    best_source_idx = kk;
                    best_dest_idx = i;
                    feasible = true;
                }
            }
            time += T_ptr[f(bef, aft, l)];
            bef = aft;
        }

        if (n != 0 && time > orders_ptr[f(bef-K, 1, 3)]) continue;
        if (time + T_ptr[f(bef, x+K, l)] <= deadline_x) {
            int new_dist = dist + dist_mat_ptr[f(bef, x+K, l)];
            float new_cost = rider[1] + rider[2] * (float)new_dist / 100;

            if (new_cost < best_cost) {
                best_cost = new_cost;
                best_source_idx = kk;
                best_dest_idx = n;
                feasible = true;
            }
        }
    }
    if (!feasible) return make_tuple(false, INF, vint{}, vint{});
    return make_tuple(true, best_cost, source_sampled[best_source_idx], dest_sampled[best_dest_idx]);
}


Res investigate_old(
    int id_to_insert, 
    const vint & source, 
    const vint & dest,
    Rider_Info & rider_info,
    int * orders,
    int * dist_mat, 
    int l
    )
{
    bool w_f, b_f, c_f;
    float w_c, b_c, c_c;
    vint w_s, w_d, b_s, b_d, c_s, c_d;

    pair<int * , vint> walkinfo = rider_info.prepare("WALK");
    pair<int * , vint> bikeinfo = rider_info.prepare("BIKE");
    pair<int * , vint> carinfo = rider_info.prepare("CAR");

    tie(w_f, w_c, w_s, w_d) = find_best_path_old(id_to_insert, source, dest, walkinfo.first, walkinfo.second, orders, dist_mat, l);
    tie(b_f, b_c, b_s, b_d) = find_best_path_old(id_to_insert, source, dest, bikeinfo.first, bikeinfo.second, orders, dist_mat, l);
    tie(c_f, c_c, c_s, c_d) = find_best_path_old(id_to_insert, source, dest, carinfo.first, carinfo.second, orders, dist_mat, l);
    bool feasible = (w_f || b_f || c_f);

    // capa : 100 정도 scale
    // rider : 10000 정도 scale
    //  float w_C = 0.7 * w_c/6 + 0.3 * walkinfo.second[0];
    //  float b_C = 0.7 * b_c/6 + 0.3 * bikeinfo.second[0];
    //  float c_C = 0.7 * c_c/6 + 0.3 * carinfo.second[0];
    float w_C = w_c;
    float b_C = b_c;
    float c_C = c_c;

    
    vector<string> opt_order;
    if (w_C <= b_C && w_C < c_C)
    {
        if(b_C <= c_C)opt_order = {"WALK", "BIKE","CAR"};
        else opt_order = {"WALK", "CAR","BIKE"};
    }
    else if (b_C <= w_C && b_C < c_C)
    {
        if(w_C <= c_C)opt_order = {"BIKE", "WALK","CAR"};
        else opt_order = {"BIKE", "CAR","WALK"};
    }
    else
    {
        if(w_C <= b_C)opt_order = {"CAR", "WALK","BIKE"};
        else opt_order = {"CAR", "BIKE","WALK"};
    }
    tuple<bool, float, vint, vint> walk = tie(w_f, w_c, w_s, w_d);
    tuple<bool, float, vint, vint> bike = tie(b_f, b_c, b_s, b_d);
    tuple<bool, float, vint, vint> car = tie(c_f, c_c, c_s, c_d);

    return Res(walk, bike, car, opt_order, feasible);
}

/* ######################################################################################################### */

tuple<bool, float, vector<int>, vector<int>> find_best_path(int x,
                                                            const vector<int> & curr_source,
                                                            const vector<int> & curr_dest, 
                                                            int * T_ptr, 
                                                            vector<int> & rider, 
                                                            int * orders_ptr, // n by 3 ndarray
                                                            int * dist_mat_ptr, 
                                                            int l)
{
    // 1. volume constraint check
    int total_volume = orders_ptr[f(x,2,3)];
    for (int id : curr_source) total_volume += orders_ptr[f(id,2,3)];
    if (total_volume > rider[0]) return make_tuple(false, INT_MAX, std::vector<int>{}, std::vector<int>{});

    int K = l/2, n = sz(curr_source);
    vint s(n+1), d(n+1);
    s[0] = x, d[0] = x;
    REP0(i, n) s[i+1] = curr_source[i], d[i+1] = curr_dest[i];
    
    if (n == 0) {
        if (orders_ptr[f(x,0,3)] + T_ptr[f(x,x+K,l)] > orders_ptr[f(x,1,3)])
            return make_tuple(false, INT_MAX, std::vector<int>{}, std::vector<int>{});
        
        int dist = dist_mat_ptr[f(x,x+K,l)];
        float cost = rider[1] + rider[2] * (float)dist / 100;
        return make_tuple(true, cost, s, d);
    }
    
    // 하려는 것.
    // source_time : source의 각 자리에 x를 넣었을 때, source의 마지막 order를 가장 빠르게 pickup할 수 있는 시각
    // dest_time : dest의 각 자리에 x를 넣었을 때, dest의 처음 order를 가장 늦게 delivery할 수 있는 시각
    // source_time + (마지막 source order -> 처음 dest order) <= dest_time 일 때만 insertion이 feasible.
    
    // 2. get source_time
    vint source_time(n+1);
    for (int i = 0; i <= n; ++i) {
        int time = orders_ptr[f(s[0],0,3)];
        for (int j = 1; j <= n; ++j) time = max(time + T_ptr[f(s[j-1], s[j], l)], orders_ptr[f(s[j],0,3)]);
        source_time[i] = time;
        if (i != n) swap(s[i], s[i+1]);
    }
    
    // 3. get dest_time
    vint dest_time(n+1);
    for (int i = 0; i <= n; ++i) {
        int time = orders_ptr[f(d[n],1,3)];
        for (int j = n-1; j >= 0; --j) time = min(time - T_ptr[f(d[j]+K, d[j+1]+K, l)], orders_ptr[f(d[j],1,3)]);
        dest_time[i] = time;
        if (i != n) swap(d[i], d[i+1]);
    }
    
    // 4. check feasibility & find best insertion
    int best_incre = INT_MAX, best_source_idx = -1, best_dest_idx = -1;
    for (int i = 0; i <= n; ++i) {
        // source에서 dist 증가량
        int source_incre = 0;
        if (i != 0) source_incre += dist_mat_ptr[f(s[i-1], x, l)];
        if (i != n) source_incre += dist_mat_ptr[f(x, s[i], l)];
        if (i != 0 && i != n) source_incre -= dist_mat_ptr[f(s[i-1], s[i], l)];
            
        for (int j = 0; j <= n; ++j) {
            // time constraint
            int bef = (i == n ? x : s[n-1]), aft = (j == 0 ? x+K : d[0]+K);
            if (source_time[i] + T_ptr[f(bef, aft, l)] > dest_time[j]) continue;
            
            // dest에서 dist 증가량
            int dest_incre = 0;
            if (j != 0) dest_incre += dist_mat_ptr[f(d[j-1]+K, x+K, l)];
            if (j != n) dest_incre += dist_mat_ptr[f(x+K, d[j]+K, l)];
            if (j != 0 && j != n) dest_incre -= dist_mat_ptr[f(d[j-1]+K, d[j]+K, l)];
            
            // source, dest 사이에서 dist 증가량
            int btw_incre = dist_mat_ptr[f(bef, aft, l)] - dist_mat_ptr[f(s[n-1], d[0]+K, l)];
            
            int total_incre = source_incre + dest_incre + btw_incre;
            if (total_incre < best_incre) {
                best_incre = total_incre;
                best_source_idx = i;
                best_dest_idx = j;
            }
        }
    }
    if (best_incre == INT_MAX) return make_tuple(false, INF, vint{}, vint{});
    
    // get inserted source & dest
    for (int i = n-1; i >= best_source_idx; --i) swap(s[i], s[i+1]);
    for (int i = n-1; i >= best_dest_idx; --i) swap(d[i], d[i+1]);
    
    // 5. get total distance with best insertion
    int dist = 0;
    for (int i = 1; i <= n; ++i) dist += dist_mat_ptr[f(s[i-1], s[i], l)];
    dist += dist_mat_ptr[f(s[n], d[0]+K, l)];
    for (int i = 1; i <= n; ++i) dist += dist_mat_ptr[f(d[i-1]+K, d[i]+K, l)];
    
    // 6. get cost
    float cost = rider[1] + rider[2] * (float)dist / 100;
    
    return make_tuple(true, cost, s, d);
}


Res investigate(
    int id_to_insert, 
    const vint & source, 
    const vint & dest,
    Rider_Info & rider_info,
    int * orders,
    int * dist_mat, 
    int l
    )
{
    bool w_f, b_f, c_f;
    float w_c, b_c, c_c;
    vint w_s, w_d, b_s, b_d, c_s, c_d;

    pair<int * , vint> walkinfo = rider_info.prepare("WALK");
    pair<int * , vint> bikeinfo = rider_info.prepare("BIKE");
    pair<int * , vint> carinfo = rider_info.prepare("CAR");

    tie(w_f, w_c, w_s, w_d) = find_best_path(id_to_insert, source, dest, walkinfo.first, walkinfo.second, orders, dist_mat, l);
    tie(b_f, b_c, b_s, b_d) = find_best_path(id_to_insert, source, dest, bikeinfo.first, bikeinfo.second, orders, dist_mat, l);
    tie(c_f, c_c, c_s, c_d) = find_best_path(id_to_insert, source, dest, carinfo.first, carinfo.second, orders, dist_mat, l);
    bool feasible = (w_f || b_f || c_f);

    // capa : 100 정도 scale
    // rider : 10000 정도 scale
    //  float w_C = 0.7 * w_c/6 + 0.3 * walkinfo.second[0];
    //  float b_C = 0.7 * b_c/6 + 0.3 * bikeinfo.second[0];
    //  float c_C = 0.7 * c_c/6 + 0.3 * carinfo.second[0];
    float w_C = w_c;
    float b_C = b_c;
    float c_C = c_c;

    
    vector<string> opt_order;
    if (w_C <= b_C && w_C < c_C)
    {
        if(b_C <= c_C)opt_order = {"WALK", "BIKE","CAR"};
        else opt_order = {"WALK", "CAR","BIKE"};
    }
    else if (b_C <= w_C && b_C < c_C)
    {
        if(w_C <= c_C)opt_order = {"BIKE", "WALK","CAR"};
        else opt_order = {"BIKE", "CAR","WALK"};
    }
    else
    {
        if(w_C <= b_C)opt_order = {"CAR", "WALK","BIKE"};
        else opt_order = {"CAR", "BIKE","WALK"};
    }
    tuple<bool, float, vint, vint> walk = tie(w_f, w_c, w_s, w_d);
    tuple<bool, float, vint, vint> bike = tie(b_f, b_c, b_s, b_d);
    tuple<bool, float, vint, vint> car = tie(c_f, c_c, c_s, c_d);

    return Res(walk, bike, car, opt_order, feasible);
}


/* ######################################################################################################### */

Solution input_to_solution(int * sol_edge_ptr,
                           int * sol_rider_ptr,
                           Rider_Info & rider_info,
                           int * dist_mat_ptr,
                           map<string, int> & riders_available,
                           int K)
{
    int l = 2*K;
    Solution sol = Solution();
    for (int i = 0; i < K; ++i)
    {
        bool ok = true;
        REP0(j, K) if (sol_edge_ptr[f(j,i,l)]) { ok = false; break; }
        if (!ok) continue;
        
        string rider_type = RIDER[sol_rider_ptr[i]];
        vint rider = rider_info.get_info(rider_type);
        
        vint merged = {i};
        while (true) {
            int x = merged.back(), y = -1;
            REP0(j, l) if (sol_edge_ptr[f(x,j,l)]) { y = j; break; }
            if (y == -1) break;
            merged.push_back(y);
        }
        
        int n = sz(merged) / 2;
        vint source(n), dest(n);
        REP0(j, n) source[j] = merged[j];
        REP0(j, n) dest[j] = merged[j+n] - K;
        
        float cost = get_bundle_cost(source, dest, rider, dist_mat_ptr, l);
        sol.append({rider_type, cost, source, dest});
        riders_available[rider_type]--;
    }
    return sol;
}

#endif