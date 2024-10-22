#ifndef UTILS_HPP
#define UTILS_HPP

#include <bits/stdc++.h>
using namespace std;


double get_time(std::chrono::time_point<std::chrono::high_resolution_clock> start) {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    return duration.count() / 1000000.0;
}



void print(map<string, int> & m)
{
    for(auto & x : m)
    {
        cout<<x.first<<" : "<<x.second<<" , ";
    }
    cout<<endl;
}


#endif