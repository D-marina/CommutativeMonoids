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

long Lambda1(std::vector<long>);
long Lambda1(std::vector<long> lgen, long Ns);

long Lambda2(std::vector<long>);
long Lambda2(std::vector<long> lgen, long Ns);

long ComputeN0(std::vector<long> lgen);
long ComputeN0(std::vector<long> lgen, long Ns);

std::vector<std::vector<long>> f1(long e, long n);

std::vector<long> W(std::vector<long> smg, long n);

std::vector<long> L(std::vector<long> lgen, long x);

std::vector<long> Nu(std::vector<long> smg, long n);

std::vector<long> Delta(std::vector<long> laux);

std::vector<long> ComputeDeltaNu(std::vector<long> lgen, long n);

std::vector<long> OrdenaSet(std::vector<long> l1);

#endif

