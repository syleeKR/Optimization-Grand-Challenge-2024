#ifndef HYPERPARAM_HPP
#define HYPERPARAM_HPP

struct Hyperparameter
{
    // Time for polish
    double time_margin = 0.5; // need fix when n>300 

    // Adaptive Destroyer
    float rho = 0.5; //updating parameter
    float smoothing_ratio = 0.01; // smoothing for insertion operators
    float min_prob = 0.05; // destroy prob min 
    float max_prob = 0.15; // destroy prob max
    int max_destroy = 1000; //just infinite value
    float route_removal_ratio = 0.2; // route removal ratio
    int update_period = 100; // period for updating prob

    int score1 = 20; // operator score when a better overall solution is found
    int score2 = 10; // operator score when a better local solution is found
    int score3 = 2;  // operator score when a worst-cost solution is accepted
    
    //// Shaw removal
    float shaw_noise = 6;
    float shaw_d = 9;
    float shaw_t = 3;
    float shaw_l = 2;
    
    //// Worst removal
    float worst_noise = 3;

    // Initial Semi-Parallel Solution. Add up to 1.
    float alpha1 = 1.0; // cost
    float alpha2 = 0.0; // empty space
    
    int consider_size = 12;
    
    float wt_weight = 0.2;
    float tw_weight = 1.0;
    
    bool use_worst = true;
    bool use_old = false;
    
    Hyperparameter() {}
    
    Hyperparameter(map<string, float> & hparam)
    {
        if (hparam.count("time_margin")) time_margin = hparam["time_margin"];
            
        if (hparam.count("rho")) rho = hparam["rho"];
        if (hparam.count("smoothing_ratio")) smoothing_ratio = hparam["smoothing_ratio"];
        if (hparam.count("min_prob")) min_prob = hparam["min_prob"];
        if (hparam.count("max_prob")) max_prob = hparam["max_prob"];
        if (hparam.count("max_destroy")) max_destroy = hparam["max_destroy"];
        if (hparam.count("route_removal_ratio")) route_removal_ratio = hparam["route_removal_ratio"];
        if (hparam.count("update_period")) update_period = hparam["update_period"];
        
        if (hparam.count("score1")) score1 = hparam["score1"];
        if (hparam.count("score2")) score2 = hparam["score2"];
        if (hparam.count("score3")) score3 = hparam["score3"];
        
        if (hparam.count("shaw_noise")) shaw_noise = hparam["shaw_noise"];
        if (hparam.count("shaw_d")) score3 = hparam["shaw_d"];
        if (hparam.count("shaw_t")) score3 = hparam["shaw_t"];
        if (hparam.count("shaw_l")) score3 = hparam["shaw_l"];
        
        if (hparam.count("worst_noise")) worst_noise = hparam["worst_noise"];
        
        if (hparam.count("alpha1")) alpha1 = hparam["alpha1"];
        if (hparam.count("alpha2")) alpha2 = hparam["alpha2"];
        
        if (hparam.count("consider_size")) consider_size = hparam["consider_size"];
        
        if (hparam.count("wt_weight")) wt_weight = hparam["wt_weight"];
        if (hparam.count("tw_weight")) tw_weight = hparam["tw_weight"];
        
        if (hparam.count("use_worst")) use_worst = hparam["use_worst"];
        if (hparam.count("use_old")) use_old = hparam["use_old"];

    }
};

#endif