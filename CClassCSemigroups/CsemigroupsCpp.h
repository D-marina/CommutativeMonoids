#ifndef _BAR_H 
#define _BAR_H 

#include <vector> 
#include <algorithm>
#include <iostream>
#include <cmath>

int foo (int a);
void Pintar(std::vector<long> v);
std::vector<long> operator-(const std::vector<long>& v1, const std::vector<long>& v2);
long gcd(long a, long b);
long gcdL(std::vector<long> v);


bool belongByGens(std::vector<long> x, std::vector<std::vector<long>> gen);
std::vector<std::vector<long>> computeMSG(std::vector<std::vector<long>> generators);
long belongAxis(std::vector<long> x,std::vector<long> r);
bool axisIsSemigroup(std::vector<std::vector<long>> gen,std::vector<long> r);
bool axisAreSemigroup(std::vector<std::vector<long>> gen,std::vector<std::vector<long>> setR);

std::vector<long> MultiplicityAxis(std::vector<std::vector<long>> generators, std::vector<long> ray);
std::vector<std::vector<long>> diamond(std::vector<std::vector<long>> mult);

//std::vector<std::vector<long>> computeGaps(std::vector<std::vector<long>> generators, std::vector<std::vector<long>> rays, std::vector<std::vector<long>> hyperplanes);

//std::vector<std::vector<long>> diamondMultiplicity(std::vector<std::vector<long>> generators, std::vector<std::vector<long>> rays);


#endif

