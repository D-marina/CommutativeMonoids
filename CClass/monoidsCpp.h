#ifndef _BAR_H 
#define _BAR_H 

#include <vector> 
#include <algorithm>
#include <iostream>
#include <cmath>

long FrobeniusNumber(std::vector<long> lgen);

std::vector<std::vector<long> > FSolve(std::vector<long> lgen, long x);

std::vector<long> smgS(std::vector<long> gen);

void Pintar(std::vector<long> v);
void Pintar(std::vector<double> v);

bool Belong(std::vector<long> generators,long x);

bool Belong(std::vector<long> generators,long x, long fNumber);

long ComputeD(std::vector<long> generators);

long ComputeNs(std::vector<long> a);

long gcd(long a, long b);

long gcdL(std::vector<long> v);

long maximum(std::vector<long> v);
double maximum(std::vector<double> v);

#endif

