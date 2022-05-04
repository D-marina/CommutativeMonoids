#ifndef _BAR_H 
#define _BAR_H 

#include <vector> 
#include <algorithm>
#include <iostream>
#include <cmath>

int foo (int a);
void foo2(std::vector<std::vector<double>> gen);
void Pintar(std::vector<long> v);
std::vector<long> operator-(const std::vector<long>& v1, const std::vector<long>& v2);
long gcd(long a, long b);
long gcdL(std::vector<long> v);
long prodEsc(const std::vector<long>& v1, const std::vector<long>& v2);
double prodEsc(const std::vector<long>& v1, const std::vector<double>& v2);

bool belongByGens(std::vector<long> x, std::vector<std::vector<long>> gen);
std::vector<std::vector<long>> computeMSG(std::vector<std::vector<long>> generators);
long belongAxis(std::vector<long> x,std::vector<long> r);
bool axisIsSemigroup(std::vector<std::vector<long>> gen,std::vector<long> r);
bool axisAreSemigroup(std::vector<std::vector<long>> gen,std::vector<std::vector<long>> setR);


std::vector<long> multiplicityAxis(std::vector<std::vector<long>> generators, std::vector<long> ray);
std::vector<std::vector<long>> multiplicityAllAxes(std::vector<std::vector<long>> generators, std::vector<std::vector<long>> rays);
std::vector<std::vector<long>> diamond(std::vector<std::vector<long>> mult);

bool pointBelongsDiamond(std::vector<long> pt, std::vector<std::vector<double>> eq);
std::vector<std::vector<long>> filterPoints(std::vector<std::vector<long>> points, std::vector<std::vector<double>> eq);
std::vector<std::vector<long>> eqRay(std::vector<long> ray, std::vector<std::vector<long>> hyperplanes);

std::vector<std::vector<long>> deleteRowZero(std::vector<std::vector<long>> m);
std::vector<std::vector<long>> affineTerm(std::vector<std::vector<long>> eq, std::vector<std::vector<long>> d);

bool studyRays(std::vector<std::vector<long>> rays, std::vector<std::vector<long>> hyperplanes, std::vector<std::vector<long>> integerDiamond, std::vector<std::vector<long>> generators);

bool existGenerator(std::vector<std::vector<long>> equationsRay, std::vector<long> affine, std::vector<std::vector<long>> generators);

//std::vector<std::vector<long>> computeGaps(std::vector<std::vector<long>> generators, std::vector<std::vector<long>> rays, std::vector<std::vector<long>> hyperplanes);

//std::vector<std::vector<long>> diamondMultiplicity(std::vector<std::vector<long>> generators, std::vector<std::vector<long>> rays);


#endif

