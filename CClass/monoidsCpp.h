#ifndef _BAR_H 
#define _BAR_H 

#include <vector> 
#include <algorithm>
#include <iostream>

long FrobeniusNumber(std::vector<long> lgen);

std::vector<std::vector<long> > FSolve(std::vector<long> lgen, long x);

std::vector<long> smgS(std::vector<long> gen);

void Pintar(std::vector<long> v);

bool Belong(std::vector<long> generators,long x);

#endif

