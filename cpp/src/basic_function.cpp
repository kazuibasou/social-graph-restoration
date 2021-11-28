#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <string.h>
#include <random>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <float.h>
#include <math.h>
#include <numeric>
#include "basic_function.h"

int generate_rand(const int N){
	if(N == 0){
		printf("Error: A given integer should be positive to generate random numbers.\n");
		exit(0);
	}
    
	srand((unsigned int)time(NULL));
	std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> randN(0, N-1);

    return randN(mt);
}

int my_round(const double value){

	//return int(std::round(value)); //simple rounding
    return std::nearbyint(value); //round to even
}

int random_key_with_smallest_value(const std::vector<std::pair<int, double>> &pairs){
    if(int(pairs.size()) == 0){
        printf("Error: The size of a given object is zero.\n");
        exit(0);
    }

    double min_value = INFINITY;
    std::vector<int> keys;
    int key;
    double value;

    for(std::pair<int, double> pair:pairs){
        key = pair.first;
        value = pair.second;
        if(value < min_value){
            min_value = double(value);
            keys.clear();
            keys.push_back(key);
        }
        else if(value == min_value){
            keys.push_back(key);
        }
    }

    return random_choice(keys);
}

int min_key_with_smallest_value(const std::vector<std::pair<int, double>> &pairs){
    if(int(pairs.size()) == 0){
        printf("Error: The size of a given object is zero.\n");
        exit(0);
    }
    
    double min_value = INFINITY;
    std::vector<int> keys;
    int key;
    double value;

    for(std::pair<int, double> pair:pairs){
        key = pair.first;
        value = pair.second;
        if(value < min_value){
            min_value = double(value);
            keys.clear();
            keys.push_back(key);
        }
        else if(value == min_value){
            keys.push_back(key);
        }
    }

    return *std::min_element(keys.begin(), keys.end());
}

int random_choice(const std::vector<int> &v){

	int index = generate_rand(int(v.size()));

	return v[index];
}
