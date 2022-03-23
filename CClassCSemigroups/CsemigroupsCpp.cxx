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

long belongAxis(std::vector<long> x,std::vector<long> r)
{
/*
   INPUT:
    - x: Value for checking if it is in the ray.
    - r: Minimal value in the ray.
  OUTPUT:
    - 0: If not belongs to the ray.
    - The coefficient if it belongs
*/
    
    /*
        coef = 0
    for i in range(len(x)):
        if x[i] != 0 and r[i] != 0:
            coef = x[i]/r[i]
            break
    if coef == 0:
        return 0
    aux2 = [j/coef for j in x]
    if aux2 == r:
        if(int(coef)==coef):
            return int(coef)
        return coef
    else:
        return 0
    */
    
    unsigned n;
    int coef;
    coef = 0;
    n = x.size();
    
    for(unsigned ii=0;ii<n;ii++)
    {
        if(x[ii] != 0 && r[ii] != 0)
        {
            cout<<"x="<<double(x[ii])<<", r="<<r[ii]<<endl;
            coef = double(x[ii])/r[ii];
            break;
        }
    }
    
    cout <<"coef = "<<coef<<endl;
    
    if(coef == 0)
    {
        return 0;
    }
    vector<long> aux2(n);
    for(unsigned ii=0;ii<n;ii++)
    {
           aux2.push_back(x[ii]/coef);
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