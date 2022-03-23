#ifndef _BAR_H 
#define _BAR_H 

#include <vector> 
#include <algorithm>
#include <iostream>
#include <cmath>

int foo (int a);
void Pintar(std::vector<long> v);
std::vector<long> operator-(const std::vector<long>& v1, const std::vector<long>& v2);

bool belongByGens(std::vector<long> x, std::vector<std::vector<long>> gen);
std::vector<std::vector<long>> computeMSG(std::vector<std::vector<long>> generators);

#endif

