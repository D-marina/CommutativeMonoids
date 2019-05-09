#include "monoidsCpp.h" 

using namespace std;

long FrobeniusNumber(vector<long> generators){
	unsigned int eDimension=generators.size();
	vector<long> lGen(generators);
	long h,t,a;
	h=t=a=generators[0];
	vector<long> S(a*lGen[lGen.size()-1]);
	vector<long> Q(a);
	vector<long> P(a);
	vector<long> aux(a);
		
	for(unsigned int i=0;i<Q.size();i++)
		Q[i]=0;
	for(unsigned int i=0;i<S.size();i++)
		S[i]=a*lGen[eDimension-1];
	S[a-1]=0;
	for(unsigned int i=0;i<P.size();i++)
		P[i]=i;
	P[a-1]=lGen.size();

	while(h!=0){
		long QHaux=Q[h-1];
		long v=h;
		Q[h-1]=0;
		if(h==t)
			h=0;
		else
			h=QHaux;
		for(int j=1;j<P[v-1]+1;j++){
			long e=(lGen[j-1]+v) % a;
			long w=lGen[j-1]+S[v-1];
			if(w<S[e-1] && e!=0){
				S[e-1]=w;
				P[e-1]=j;
				if(Q[e-1]==0){
					if(h==0){
						h=e;
						Q[e-1]=e;
						t=e;
					}
					else
					{
						Q[t-1]=e;
						Q[e-1]=e;
						t=e;
					}
				}
			}
		}
	}
	for(unsigned int i=0;i<a;i++)
		aux[i]=S[i];
	long fNumber=aux[0];
	for(unsigned int i=0;i<aux.size();i++)
		if(aux[i]>fNumber)
			fNumber=aux[i];
	return fNumber-a;
}

int isZero(vector<long> v){
	for(unsigned int i=0;i<v.size();i++)
		if(v[i]!=0)
			return 0;
	return 1;
}


///


vector<vector<long> > FSolve(vector<long> lgen, long x)
{
    int posAModificar=0;
    bool sumando = true;
    int xaux = x;
    int dim = lgen.size();

    vector<long> tuplaActual(dim);    
    vector<vector<long> > soluciones;

    for(unsigned int i=0;i<tuplaActual.size();++i){
	tuplaActual[i]=0;
    }

    int k=0;

    while( sumando || !isZero(tuplaActual) ) 
    {
	k++;
        if(sumando)
        {
            if(xaux>=lgen[0])
            {
                tuplaActual[posAModificar]=tuplaActual[posAModificar]+1;
                xaux=xaux-lgen[posAModificar];
                if(xaux==0)
                {
                    soluciones.push_back(tuplaActual);
                    sumando = false;
                    continue;
                }
            	}
            	if(xaux!=0 && xaux<lgen[0])
            	{            
                	sumando = false;
                	continue;
            	}
        }
        else
        {
            if(( posAModificar == dim-1 ) && ( tuplaActual[posAModificar] > 0 ))
            {
                xaux = tuplaActual[posAModificar]*lgen[posAModificar]+xaux;
                tuplaActual[posAModificar] = 0;
                posAModificar = posAModificar-1;
                continue;
            }
            if(( posAModificar < dim-1 ) && ( tuplaActual[posAModificar] > 0 ))
            {
                xaux = xaux+lgen[posAModificar];
                tuplaActual[posAModificar] = tuplaActual[posAModificar]-1;
                posAModificar = posAModificar+1;
                sumando = true;
                continue;
            }
            if (( posAModificar < dim-1 ) && ( tuplaActual[posAModificar] == 0 ))
            {
		if(tuplaActual.size()==0)
			posAModificar=0;
		for(unsigned int i=0;i<tuplaActual.size();i++)
		{	
			if(tuplaActual[i]!=0)
			{
				posAModificar = i;
			}
		}
                xaux = lgen[posAModificar]+xaux;
                tuplaActual[posAModificar] = tuplaActual[posAModificar]-1;
                sumando = true;
                posAModificar = posAModificar+1;
                continue;
            }
        }
    }
    return soluciones;
}


///


vector<long> smgS(vector<long> gen)
{
	vector<long> smgS, generators;
	generators = gen;
	sort(generators.begin(),generators.end());
	long aux;
	aux = generators[0];
	for(unsigned i=1;i<generators.size();i++)
	{
		if(aux == generators[i])
		{
			generators.erase(generators.begin()+i);
			i=1;
		}
		aux=generators[i];
	}
	return generators;
}
