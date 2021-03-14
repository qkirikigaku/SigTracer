#pragma once
#include <random>

using namespace std;

namespace Sampler{
    random_device SeedGen;
    mt19937 mt(SeedGen());
    
    double uniform(double min, double max){
        uniform_real_distribution<double> rand(min, max);
        return rand(mt);
    }

}

