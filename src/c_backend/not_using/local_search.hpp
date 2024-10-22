#ifndef LOCAL_SEARCH_HPP
#define LOCAL_SEARCH_HPP

#include "helper_classes.hpp"
#include "helper_functions.hpp"
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

struct Seq {
     // Duration, Time warp, Earlist visit, Latest visit, Cumulated distance, Load
    int D, TW, E, L, C, Q;
    string rider_type;
    
    Seq() {}
    Seq(int d, int tw, int e, int l, int c, int q, string r): D(d), TW(tw), E(e), L(l), C(c), Q(q), rider_type(r) {}
};

struct SequenceInfo {
    vint source, dest;
    Seq source_seq, dest_seq;
    
    SequenceInfo(vint & source, vint & dest, Seq & source_seq, Seq & dest_seq):
        source(source), dest(dest), source_seq(source_seq), dest_seq(dest_seq)
    {}
    
    tuple<vint,vint,Seq,Seq> extract() { return {source, dest, source_seq, dest_seq}; }
};

struct LocalSearch {
    int K, l;
    int * orders_ptr; //orders array 1D
    Rider_Info rider_info; // info about each riders
    int * dist_mat_ptr; // distance matrix array 1D
    default_random_engine rng;
    
    map<string,int> rider_to_num = {{"CAR", 0}, {"BIKE", 1}, {"WALK", 2}};
    vector<bool> overlap;
    
    LocalSearch(
        int K, 
        int * orders_ptr,
        Rider_Info & rider_info,
        int * dist_mat_ptr,
        int promising_size = 20,
        int seed = 42
    ): K(K), orders_ptr(orders_ptr), rider_info(rider_info), dist_mat_ptr(dist_mat_ptr), rng(seed)
    {
        l = 2*K;
        overlap.resize(K, false);
    }
    
    //////////////////////////////////////////////////////
    
    Seq concat_seq(const Seq & s1, const Seq & s2, int from, int to)
    {
        string rider_type = s1.rider_type;
        int* T = rider_info.get_T(rider_type);
        
        int delta = s1.D - s1.TW + T[f(from,to,l)];
        int delta_WT = max(s2.E - delta - s1.L, 0);
        int delta_TW = max(s1.E + delta - s2.L, 0);
        
        int D = s1.D + s2.D + T[f(from,to,l)] + delta_WT;
        int TW = s1.TW + s2.TW + delta_TW;
        int E = max(s2.E - delta, s1.E) - delta_WT;
        int L = min(s2.L - delta, s1.L) + delta_TW;
        int C = s1.C + s2.C + dist_mat_ptr[f(from,to,l)];
        int Q = s1.Q + s2.Q;
        return Seq(D, TW, E, L, C, Q, s1.rider_type);
    }
    
    bool is_valid(const Seq & s)
    {
        int capa = rider_info.get_info(s.rider_type)[0];
        return (s.TW == 0 && s.Q <= capa);
    }
    
    Seq source_to_seq(vector<int> & vec, string rider_type)
    {
        Seq s = Seq(0, 0, orders_ptr[f(vec[0],0,3)], orders_ptr[f(vec[0],1,3)], 0, orders_ptr[f(vec[0],2,3)], rider_type);
        for (int i = 1; i < sz(vec); ++i) {
            Seq ss = Seq(0, 0, orders_ptr[f(vec[i],0,3)], orders_ptr[f(vec[i],1,3)], 0, orders_ptr[f(vec[i],2,3)], rider_type);
            s = concat_seq(s, ss, vec[i-1], vec[i]);
        }
        return s;
    }
    
    Seq dest_to_seq(vector<int> & vec, string rider_type)
    {
        Seq s = Seq(0, 0, orders_ptr[f(vec[0],0,3)], orders_ptr[f(vec[0],1,3)], 0, 0, rider_type);
        for (int i = 1; i < sz(vec); ++i) {
            Seq ss = Seq(0, 0, orders_ptr[f(vec[i],0,3)], orders_ptr[f(vec[i],1,3)], 0, 0, rider_type);
            s = concat_seq(s, ss, vec[i-1]+K, vec[i]+K);
        }
        return s;
    }
    
    //////////////////////////////////////////////////////
    
    void cut(Bundle & bundle, vector<vector<SequenceInfo>> & bundle_cuts)
    {
        string rider_type; float cost; vint source, dest;
        tie(rider_type, cost, source, dest) = bundle;
        int rider_num = rider_to_num[rider_type];
        int n = sz(source);
        
        std::uniform_real_distribution<> dis(0, 1);
        bool pick = (dis(rng) < 0.5);
        
        vector<bool> chk(n);
        if (pick) {
            for (int i = 0; i < n-1; ++i) {
                vint s1, s2;
                for (int j=0; j<=i; ++j) s1.pb(source[j]);
                for (int j=i+1; j<n; ++j) s2.pb(source[j]);

                for (int j=0; j<n; ++j) {
                    if (dest[j] == source[i]) chk[j] = true;
                }

                vint d1, d2;
                for (int j=0; j<n; ++j) {
                    if (chk[j]) d1.pb(dest[j]);
                    else d2.pb(dest[j]);
                }

                Seq s1_seq = source_to_seq(s1, rider_type), d1_seq = dest_to_seq(d1, rider_type);
                Seq s2_seq = source_to_seq(s2, rider_type), d2_seq = dest_to_seq(d2, rider_type);
                
                bundle_cuts[rider_num].pb({s1, d1, s1_seq, d1_seq});
                bundle_cuts[rider_num].pb({s2, d2, s2_seq, d2_seq});
            }
        }
        else {
            for (int i = 0; i < n-1; ++i) {
                vint d1, d2;
                for (int j=0; j<=i; ++j) d1.pb(dest[j]);
                for (int j=i+1; j<n; ++j) d2.pb(dest[j]);

                for (int j=0; j<n; ++j) {
                    if (source[j] == dest[i]) chk[j] = true;
                }

                vint s1, s2;
                for (int j=0; j<n; ++j) {
                    if (chk[j]) s1.pb(source[j]);
                    else s2.pb(source[j]);
                }

                Seq s1_seq = source_to_seq(s1, rider_type), d1_seq = dest_to_seq(d1, rider_type);
                Seq s2_seq = source_to_seq(s2, rider_type), d2_seq = dest_to_seq(d2, rider_type);

                bundle_cuts[rider_num].pb({s1, d1, s1_seq, d1_seq});
                bundle_cuts[rider_num].pb({s2, d2, s2_seq, d2_seq});
            }
        }
    }
    
    Seq concat(vint & s1, vint & s2, vint & d1, vint & d2, Seq & s_seq1, Seq & s_seq2, Seq & d_seq1, Seq & d_seq2)
    {
        Seq s = s_seq1;
        s = concat_seq(s, s_seq2, s1.back(), s2.front());
        s = concat_seq(s, d_seq1, s2.back(), d1.front()+K);
        s = concat_seq(s, d_seq2, d1.back()+K, d2.front()+K);
        return s;
    }
    
    Bundle merge(vint & s1, vint & s2, vint & d1, vint & d2, Seq & seq)
    {
        string rider_type = seq.rider_type;
        
        vint rider = rider_info.get_info(rider_type);
        float cost = rider[1] + rider[2] * (float)seq.C / 100;
        
        vint source = s1; for (int x : s2) source.pb(x);
        vint dest = d1; for (int x : d2) dest.pb(x);
        return {rider_type, cost, source, dest};
    }
    
    void mix(SequenceInfo & info1, SequenceInfo & info2, BundleStorage & bundle_storage)
    {
        vint s1, d1; Seq s_seq1, d_seq1; tie(s1, d1, s_seq1, d_seq1) = info1.extract();
        vint s2, d2; Seq s_seq2, d_seq2; tie(s2, d2, s_seq2, d_seq2) = info2.extract();
        
        for (int x : s1) overlap[x] = true;
        for (int x : s2) if (overlap[x]) return;
        for (int x : s1) overlap[x] = false;
        
        std::uniform_int_distribution<std::size_t> dist(0, 3);
        int id = dist(rng);
        
        if (id == 0) {
            // s1 s2 d1 d2
            Seq seq1 = concat(s1, s2, d1, d2, s_seq1, s_seq2, d_seq1, d_seq2);
            if (is_valid(seq1)) bundle_storage.append(merge(s1, s2, d1, d2, seq1));
        }
        else if (id == 1) {
            // s1 s2 d2 d1
            Seq seq2 = concat(s1, s2, d2, d1, s_seq1, s_seq2, d_seq2, d_seq1);
            if (is_valid(seq2)) bundle_storage.append(merge(s1, s2, d2, d1, seq2));
        }
        else if (id == 2) {
            // s2 s1 d1 d2
            Seq seq3 = concat(s2, s1, d1, d2, s_seq2, s_seq1, d_seq1, d_seq2);
            if (is_valid(seq3)) bundle_storage.append(merge(s2, s1, d1, d2, seq3));
        }
        else if (id == 3) {
            // s2 s1 d2 d1
            Seq seq4 = concat(s2, s1, d2, d1, s_seq2, s_seq1, d_seq2, d_seq1);
            if (is_valid(seq4)) bundle_storage.append(merge(s2, s1, d2, d1, seq4));
        }
    }
    
    void search(Solution & sol, BundleStorage & bundle_storage) // 2-opt
    {
        vector<vector<SequenceInfo>> bundle_cuts(3);
        for (int bundle_id = 0; bundle_id < sz(sol.solutions); bundle_id++) {
            std::uniform_real_distribution<> dis(0, 1);
            bool pick = (dis(rng) < 0.3);
            if (!pick) continue;
            
            Bundle bundle = sol.solutions[bundle_id];
            cut(bundle, bundle_cuts);
        }
        
        REP0(rider_num, 3) {
            int n = sz(bundle_cuts[rider_num]);
            vint indices(n); REP0(i, n) indices[i] = i;
            std::shuffle(all(indices), rng);

            for (int ii = 0; ii < min(n, 20); ++ii) {
                int i = indices[ii];
                for (int j = 0; j < n; ++j) {
                    if (i == j) continue;
                    mix(bundle_cuts[rider_num][i], bundle_cuts[rider_num][j], bundle_storage);
                }
            }
        }
    }
};

#endif