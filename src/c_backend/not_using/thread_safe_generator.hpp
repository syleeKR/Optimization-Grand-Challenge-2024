#ifndef THREAD_SAFE_GENERATOR
#define THREAD_SAFE_GENERATOR

#include <iostream>
#include <random>
#include <mutex>

struct ThreadSafeRandomEngine {
    using result_type = std::default_random_engine::result_type;
    
    std::default_random_engine rng;
    std::mutex rng_mutex;

    // Constructor that accepts a seed
    ThreadSafeRandomEngine(int seed = 42) : rng(seed) {}

    // Overload the () operator to act like a function
    result_type operator()() {
        std::lock_guard<std::mutex> lock(rng_mutex);  // Ensure thread-safe access to rng
        return rng();
    }

    // Conform to the UniformRandomBitGenerator concept
    static constexpr result_type min() {
        return std::default_random_engine::min();
    }

    static constexpr result_type max() {
        return std::default_random_engine::max();
    }

    // Optional: A function to reset the seed if needed
    void seed(int new_seed) {
        std::lock_guard<std::mutex> lock(rng_mutex);
        rng.seed(new_seed);
    }
};

#endif