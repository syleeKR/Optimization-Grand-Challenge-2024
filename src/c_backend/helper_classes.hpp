#ifndef HELPER_CLASSES_HPP
#define HELPER_CLASSES_HPP


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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

/* ################################################################################################ */

const float INF = 1e9;
const vector<string> RIDER = {"WALK", "BIKE", "CAR"};

typedef pair<string,vint> KEY;
typedef tuple<string, float, vector<int>, vector<int>> Bundle;
typedef pair<float, vector<tuple<string, vint, vint>>> SOLUTION_FORMAT;
typedef tuple<vector<int>, vector<vector<int>>, vector<Bundle>> GUROBI_INPUT;

struct Rider_Info
{
    int * walk_T_ptr;
    vint walk_info;
    int * bike_T_ptr;
    vint bike_info;
    int * car_T_ptr;
    vint car_info;
    Rider_Info(int * walk_T_ptr,
                    vint & walk_info,
                    int * bike_T_ptr,
                    vint & bike_info,
                    int * car_T_ptr,
                    vint & car_info):walk_T_ptr(walk_T_ptr),walk_info(walk_info), bike_T_ptr(bike_T_ptr),bike_info(bike_info),car_T_ptr(car_T_ptr),car_info(car_info)
    {
    }

    pair<int * , vint> prepare(string s)
    {
        if(s=="WALK") return {walk_T_ptr, walk_info};
        else if(s=="BIKE") return {bike_T_ptr, bike_info};
        else return {car_T_ptr, car_info};
    }
    
    int* get_T(string s)
    {
        if(s=="WALK") return walk_T_ptr;
        else if(s=="BIKE") return bike_T_ptr;
        else return car_T_ptr;
    }
    
    vint get_info(string s)
    {
        if(s=="WALK") return walk_info;
        else if(s=="BIKE") return bike_info;
        else return car_info;
    }
};

struct Res {
    tuple<bool, float, vint, vint> walk;
    tuple<bool, float, vint, vint> bike;
    tuple<bool, float, vint, vint> car;
    vector<string> optimal_order;
    bool feasible;

    Res(){}

    Res(tuple<bool, float, vint, vint> & w, tuple<bool, float, vint, vint> & b, tuple<bool, float, vint, vint> &c, vector<string> & opt_order, bool f)
    {
        walk = w;
        bike = b;
        car = c;
        optimal_order = opt_order;
        feasible = f;
    }
    tuple<bool, float, vint, vint> get_all(string s)
    {
        if (s == "WALK")return walk;
        else if(s=="BIKE")return bike;
        else return car;
    }
    bool get_feasibility(string s)
    {
        if (s == "WALK")return get<0>(walk);
        else if (s=="BIKE")return get<0>(bike);
        else return get<0>(car);
    }
    float get_cost(string s)
    {
        if (s == "WALK")return get<1>(walk);
        else if (s=="BIKE")return get<1>(bike);
        else return get<1>(car);
    }
    vint get_source(string s)
    {
        if (s == "WALK")return get<2>(walk);
        else if (s=="BIKE")return get<2>(bike);
        else return get<2>(car);
    }
    vint get_dest(string s)
    {
        if (s == "WALK")return get<3>(walk);
        else if (s=="BIKE")return get<3>(bike);
        else return get<3>(car);
    }
};


struct Cache {
    // Member variables
    map<vint, Res> storage;

    Cache(){}
    bool check(vint v){sort(all(v));return storage.count(v);}
    void append(vint v, Res r){sort(all(v));storage[v] = r;}
    Res retrieve(vint v){sort(all(v)); return storage[v];}
};


struct Solution
{
    float cost;
    vector<Bundle> solutions;

    Solution():cost(0.0){}
    Solution(const pair<vector<Bundle>, float> & sol_raw)
    {
        cost = sol_raw.second;
        solutions = sol_raw.first;
    }
    string get_rider_type(int i){return get<0>(solutions[i]);}
    float get_cost(int i){return get<1>(solutions[i]);}
    vint get_source(int i){return get<2>(solutions[i]);}
    vint get_dest(int i){return get<3>(solutions[i]);}

    int get_bundle_id(int order_id) // if no bundle return -1
    {
        REP0(i, sz(solutions))
        {
            vint orders = get_source(i);
            for(auto o : orders)
            {
                if ( o == order_id)return i;
            }
        }
        return -1;
    }

    
    void update_cost()
    {
        cost = 0;
        for (auto x : solutions)
        {
            cost += get<1>(x);
        }
    }
    void remove(vint & ids) //remove bundle!
    {
        set<int> tmp; for(auto x : ids)tmp.insert(x);
        vector<tuple<string, float, vint, vint>> new_solutions;
        REP0(i, sz(solutions))
        {
            if(tmp.count(i)==0)new_solutions.pb(solutions[i]);
        }
        solutions = new_solutions;
        update_cost();
    }
    void remove(int i) //remove bundle!
    {
        cost -= get<1>(solutions[i]);
        solutions.erase(solutions.begin() + i);
    }
    pair<bool, string> remove_order(int order_id, int bundleid, float newcost)
    {  
        string rider_type = get_rider_type(bundleid);
        vint source = get_source(bundleid);

        if (sz(source)==1){ remove(bundleid); return {true,rider_type}; }  

        vint dest = get_dest(bundleid);
        source.erase(find(all(source), order_id));
        dest.erase(find(all(dest), order_id));
        cost += (newcost - get<1>(solutions[bundleid]));
        solutions[bundleid] = make_tuple(rider_type, newcost, source, dest);
        return {false,rider_type};
    }
    void append(Bundle x)
    {
        solutions.push_back(x);
        cost += get<1>(x);
    }


 /*####################################################################*/
    pair<vector<Bundle>, float> get_sol()
    {
        return {solutions,cost};
    }
    pair<float, vector<tuple<string, vint, vint>>> extract()
    {
        vector<tuple<string, vint, vint>> ans;

        for(auto x  : solutions)
        {
            ans.pb(make_tuple(get<0>(x), get<2>(x),get<3>(x)));
        }
        return {cost, ans};
    }
};



struct BundleStorage {

    map<KEY, Bundle> storage;

    BundleStorage(){}

    BundleStorage(vector<Bundle> & init_storage)
    {
        for (Bundle& bundle : init_storage) {
            KEY key = get_key(bundle);
            storage[key] = bundle;
        }
    }
    
    KEY get_key(Bundle & bundle)
    {
        string rider_type = get<0>(bundle);
        vint orders = get<2>(bundle);
        sort(all(orders));

        KEY key = {rider_type, orders};
        return key;
    }
    
    void append(Bundle bundle)
    {
        KEY key = get_key(bundle);
        if(storage.count(key) ==0) storage[key] = bundle;
        else
        {
            float cost_of_bundle = get<1>(bundle);
            float cost_of_stored_bundle = get<1>(storage[key]);
            if(cost_of_bundle < cost_of_stored_bundle)storage[key] = bundle;
        }
    }

    void append(Solution & sol, float prob = 1.0)
    {
        // best sol이나 accept된 애들 prob는 1.0으로 나머지는 작게?
        // best sol일때는 반드시 1.0이어야 함. gurobi에 initial sol 줘야돼서
        // 나중을 위해서 일단 prob도 넣어놓음.
        // 2가지 고쳐야함
        // 1.unordered_map 으로 근데 custom hash 짜야됨
        // 2.py::array_t로 c++ -> python 넘기기 근데 충분히 빨라서 안해도될듯
        for(Bundle bundle : sol.solutions)
        {
            KEY key = get_key(bundle);

            if(storage.count(key) ==0)storage[key] = bundle;
            else
            {
                float cost_of_bundle = get<1>(bundle);
                float cost_of_stored_bundle = get<1>(storage[key]);
                if(cost_of_bundle < cost_of_stored_bundle)storage[key] = bundle;
            }
            
        }
    }
    
    /*
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
    */

    tuple<vector<int>, vector<vector<int>> , vector<Bundle>> extract(Solution & bestsol, int K, Rider_Info & rider_info, int * dist_mat_ptr, int l)
    {
        /* 
        n : number of bundles in storage
        K : number of orders
        returns 3 things
        #1. the best solution that will be provided to gurobi, 
            vector<int> containing the bundle indexes that constitute best sol
        #2. the helper matrix that will be used in making constraints
            vector<vector<int>> of size K, each vector containing the bundle indexes that contain the order
        #3. the vector of bundles
        */
        
        /*
        auto start_time = chrono::high_resolution_clock::now();
        
        map<KEY, Bundle> tmp;
        for (const auto& pair : storage) {
            Bundle bundle = pair.second;
            string rider_type; float cost; vint source, dest;
            tie(rider_type, cost, source, dest) = bundle;
            
            vint rider = rider_info.get_info(rider_type);
            
            int n = sz(source);
            for (int i = 0; i < (1<<n); ++i) {
                vint new_source = source, new_dest = dest;
                
                vint val_to_del;
                for (int j = 0; j < n; ++j) {
                    if ((i >> j) & 1) val_to_del.pb(new_source[j]);
                }
                
                for (int x : val_to_del) new_source.erase(find(all(new_source), x));
                for (int x : val_to_del) new_dest.erase(find(all(new_dest), x));
                
                float new_cost = get_bundle_cost(source, dest, rider, dist_mat_ptr, l);
                
                Bundle new_bundle = {rider_type, new_cost, new_source, new_dest};
                
                KEY key = get_key(new_bundle);
                if(tmp.count(key) == 0) tmp[key] = new_bundle;
                else
                {
                    float cost_of_bundle = get<1>(new_bundle);
                    float cost_of_stored_bundle = get<1>(tmp[key]);
                    if(cost_of_bundle < cost_of_stored_bundle) tmp[key] = new_bundle;
                }
            }
        }
        storage = tmp;
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start_time);
        cout << "time: " << duration.count() / 1000000.0 << endl;
        
        */
        
        // vector<tuple<float,KEY,Bundle>> vv;
        // for (const auto& pair : storage) {
        //     Bundle bundle = pair.second;
        //     float cost = get<1>(bundle);
        //     int n = sz(get<2>(bundle));
        //     vv.pb({cost/n, pair.first, pair.second});
        // }
        // sort(all(vv));
        
        // map<KEY, Bundle> tmp;
        // for (int i=0; i<min(50000, sz(vv)); ++i) {
        //     tmp[get<1>(vv[i])] = get<2>(vv[i]);
        // }
        // storage = tmp;
        // append(bestsol);

        /* #3 */
        map<KEY, int> keys_to_index;
        vector<Bundle> bundles;
        int idx = 0;
        for (const auto& pair : storage) {
            bundles.pb(pair.second);
            keys_to_index[pair.first] = idx;
            idx++;
        }

        /* #1 */
        vint gurobi_initial_solution;
        for(Bundle bundle : bestsol.solutions)
        {
            KEY key = get_key(bundle);
            int keyindex = keys_to_index[key];
            gurobi_initial_solution.pb(keyindex);
        }

        /* #2 */
        vector<vector<int>> constraint_helper_matrix(K);
        for(int i = 0; i < sz(bundles); i++)
        {
            Bundle bundle = bundles[i];
            vint orders = get<2>(bundle);
            for(int order : orders) constraint_helper_matrix[order].pb(i);
        }
        return make_tuple(gurobi_initial_solution, constraint_helper_matrix, bundles);
    }
};


#endif