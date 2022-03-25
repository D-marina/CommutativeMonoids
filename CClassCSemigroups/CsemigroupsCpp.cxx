#include "CsemigroupsCpp.h" 

using namespace std;

int foo(int a)
{
    return a+1;
}

void Pintar(vector<long> v)
{
	cout<<"*";
	for(unsigned ii=0;ii<v.size();ii++)
		cout<<","<<v[ii];
	cout<<"*";
	cout<<endl;
}

vector<long> operator-(const vector<long>& v1, const vector<long>& v2)
{
    unsigned n;
    n = v1.size();
    vector<long> aux(n);
    for(unsigned ii=0;ii<n;ii++)
    {
        aux[ii] = v1[ii]-v2[ii];
    }
    return aux;
}

long gcd(long a, long b)
{
	if(a == 0)
		return b;
	return gcd(b%a,a);
}


///


long gcdL(vector<long> v) 
{
	long n, result;
	n = v.size();
	result = v[0];
	for(unsigned i = 1; i < n; i++) 
        	result = gcd(v[i], result); 
	if(result<0)
		return -result;
	return result; 
}

/////////////////////////////////////////////////////////////////////


bool belongByGens(vector<long> x, vector<vector<long>> gen)
{
    /*
    This function checks if x can be generated using the
    elementos of gen.
    */
    unsigned n,m;
    bool allZero;
    n = gen.size();
    m = x.size();
    if(n == 0)
        return false;
    allZero = true;
    for(unsigned ii=0;ii<m;ii++)
    {
        if(x[ii]!=0)
            allZero = false;
        if(x[ii]<0)
            return false;
    }
    if(allZero)
        return true;
    if(belongByGens(x-gen[0], gen))
    {
        return true;
    }
    else
    {
        if(gen.size() == 1)
        {
            return false;
        }
        else
        {
            vector<vector<long>> aux(n-1);
            for(unsigned ii=0;ii<n-1;ii++)
                aux[ii] = gen[ii+1];
            return belongByGens(x,aux);
        }
    }
    return false;
}

vector<std::vector<long>> computeMSG(vector<vector<long>> generators)
{
 /*
 This function returns a minimal sistem of generators from the set generators.
 */
    unsigned n;
    n = generators.size();
    vector<vector<long>> minimales;
    for(unsigned ii=0; ii<n; ii++)
    {
        vector<vector<long>> aux(n-1); // Creo un vector para meter todos los generadores menos el ii-ésimo.
        unsigned kk =0; // Contador para saber la posición que toca rellenar.
        for(unsigned jj=0; jj<n; jj++)
        {
            if(ii != jj)
            {
                    aux[kk] = generators[jj];
                    kk++;
            }
        }
        if( !belongByGens(generators[ii],aux))
        {
            minimales.push_back(generators[ii]);
        }
    }
    return minimales;
}

long belongAxis(vector<long> x,vector<long> r)
{
/*
   INPUT:
    - x: Value for checking if it is in the ray.
    - r: Minimal value in the ray.
  OUTPUT:
    - 0: If not belongs to the ray.
    - The coefficient if it belongs
*/    
    unsigned n;
    double coef;
    coef = 0;
    n = x.size();
    
    for(unsigned ii=0;ii<n;ii++)
    {
        if(x[ii] != 0 && r[ii] != 0)
        {
            coef = double(x[ii])/r[ii];
            break;
        }
    }
    
    if(coef == 0)
    {
        return 0;
    }
    vector<long> aux2(n);
    
    for(unsigned ii=0;ii<n;ii++)
    {
        aux2[ii] = x[ii]/coef;
    }
    
    if(aux2 == r)
    {
       if(int(coef) == coef)
       {
            return int(coef);
       }
       return coef;
    }
    else
    {
        return 0;
    }
}

bool axisIsSemigroup(vector<vector<long>> gen,vector<long> r)
{
/*
    # Vemos para un eje si se forma un semigrupo.
    # INPUT:
    #   - gen: Set of generators.
    #   - r: Minimal value in the ray.
    # OUTPUT:
    #   - True/False.   
*/
    
    vector<long> aux;
    unsigned n;
    long coef;
    n = gen.size();
    for(unsigned ii=0;ii<n;ii++)
    {
        coef = belongAxis(gen[ii],r);
        if(coef != 0)
        {
            aux.push_back(coef);
        }
    }
    if(gcdL(aux) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool axisAreSemigroup(vector<vector<long>> gen,vector<vector<long>> setR)
{
    /*
    # Vemos si los ejes forman un semigrupo.
    # INPUT:
    #   - gen: Set of generators.
    #   - setR: Set of rays.
    # OUTPUT:
    #   - True/False.
    */
    
    unsigned n;
    n = setR.size();
    for(unsigned ii=0;ii<n;ii++)
    {
            if(!axisIsSemigroup(gen,setR[ii]))
            {
                return false;
            }
    }
    return true;
}