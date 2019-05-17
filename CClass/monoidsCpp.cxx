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


vector<long> smgS(vector<long> generators)
{
	vector<long> smgS, sgS;
	sgS = generators;
	sort(sgS.begin(),sgS.end());
	long aux;
	aux = sgS[0];
	// Delete duplicate
	for(unsigned i=1;i<sgS.size();i++)
	{
		if(aux == sgS[i])
		{
			sgS.erase(sgS.begin()+i);
			i=1;
		}
		aux=sgS[i];
	}
	
	if(sgS.size()==1)
		return sgS;
	smgS.push_back(sgS[0]);
	unsigned i;
	for(i=1;i<sgS.size();i++)
	{
		if(sgS[i]%smgS[0] != 0)
		{
			smgS.push_back(sgS[i]);
			break;
		}
	}

	for(unsigned j=i+1;j<sgS.size();j++)
	{
		if(FSolve(smgS, sgS[j]).size() == 0)
		{
			smgS.push_back(sgS[j]);
		}
	}

	return smgS;
}


///


bool Belong(vector<long> generators,long x)
{
	long multiplicity, fNumber;
	multiplicity =generators[0];
	fNumber = FrobeniusNumber(generators);
	if(x == 0)
		return true;
	if((0 < x) && (x < multiplicity))
		return false;
	for(unsigned i=0;i<generators.size();i++)
	{
		if(x==generators[i])
			return true;
	}
	if((fNumber != 0) && (x>fNumber))
		return true;
	vector<vector<long> > expression;
	expression = FSolve(generators,x);
	if(expression.size() > 0)
		return true;
	else
		return false;
}


///

bool Belong(vector<long> generators,long x, long fNumber)
{
	long multiplicity;
	multiplicity =generators[0];
	if(x == 0)
		return true;
	if((0 < x) && (x < multiplicity))
		return false;
	for(unsigned i=0;i<generators.size();i++)
	{
		if(x==generators[i])
			return true;
	}
	if((fNumber != 0) && (x>fNumber))
		return true;
	vector<vector<long> > expression;
	expression = FSolve(generators,x);
	if(expression.size() > 0)
		return true;
	else
		return false;
}

///


long ComputeD(vector<long> generators)
{
	vector<long> dif;
	long dimension;
	dimension = generators.size();
	for(long i=0;i<dimension-1;i++)
	{
		dif.push_back(generators[i+1]-generators[i]);
	}
	return gcdL(dif);
}


///


long ComputeNs(vector<long> a)
{
	long p, d;
	vector<double> Si;
	vector<double> Siprime;
	p = a.size();
	d = ComputeD(a);
	double max1, max2, max;

	for(int i=1;i<p-1;i++)
	{
		vector<long> aux;
		aux.push_back(a[i]-a[0]);
		aux.push_back(a[0]-a[p-1]);
		aux.push_back(a[p-1]-a[i]);
		
		
		double aux11, aux12, aux21,aux22;
		aux11 = (double)(-a[1]*(a[0]*d*gcdL(aux)+(p-2)*(a[0]-a[i])*(a[0]-a[p-1])));
		aux12 = (double)(((a[0]-a[1])*gcdL(aux)));
		aux21 = (double)(a[p-2]*((p-2)*(a[0]-a[p-1])*(a[p-1]-a[i])-a[p-1]*d*gcdL(aux)));
		aux22 = (double)((a[p-2]-a[p-1])*gcdL(aux));

		Si.push_back(aux11/aux12);
		Siprime.push_back(aux21/aux22);
	}

	max1 = maximum(Si);
	max2 = maximum(Siprime);
	max = max1;

	if(max1 < max2)
		max = max2;
	return (long)ceil(max);
}


///


long ComputeNs(vector<long> a, long d)
{
	long p;
	vector<double> Si;
	vector<double> Siprime;
	p = a.size();
	double max1, max2, max;

	for(int i=1;i<p-1;i++)
	{
		vector<long> aux;
		aux.push_back(a[i]-a[0]);
		aux.push_back(a[0]-a[p-1]);
		aux.push_back(a[p-1]-a[i]);
		
		
		double aux11, aux12, aux21,aux22;
		aux11 = (double)(-a[1]*(a[0]*d*gcdL(aux)+(p-2)*(a[0]-a[i])*(a[0]-a[p-1])));
		aux12 = (double)(((a[0]-a[1])*gcdL(aux)));
		aux21 = (double)(a[p-2]*((p-2)*(a[0]-a[p-1])*(a[p-1]-a[i])-a[p-1]*d*gcdL(aux)));
		aux22 = (double)((a[p-2]-a[p-1])*gcdL(aux));

		Si.push_back(aux11/aux12);
		Siprime.push_back(aux21/aux22);
	}

	max1 = maximum(Si);
	max2 = maximum(Siprime);
	max = max1;

	if(max1 < max2)
		max = max2;
	return (long)ceil(max);
}


///


long Lambda1(vector<long> lgen)
{
	long Ns,dimension;
	double c11, c12, c41,c42, c1, c4;
	Ns = ComputeNs(lgen);
	dimension = lgen.size();
	c11 = (double) (lgen[dimension-1]-lgen[dimension-2])*Ns;
	c12 = (double) lgen[dimension-2];
	c1 = c11/c12;
	c41 = (double)(-Ns)*(lgen[dimension-2]*lgen[dimension-1]*lgen[0]-lgen[dimension-2]*lgen[dimension-1]*lgen[1] + lgen[dimension-2]*lgen[0]*lgen[1]-lgen[dimension-1]*lgen[0]*lgen[1]);
	c42 = (double) lgen[dimension-2]*lgen[dimension-1]*lgen[1];
	c4= c41/c42;	

	return max(ceil(c1),ceil(c4));	
}

///


long Lambda1(vector<long> lgen, long Ns)
{
	long dimension;
	double c11, c12, c41,c42, c1, c4;
	dimension = lgen.size();
	c11 = (double) (lgen[dimension-1]-lgen[dimension-2])*Ns;
	c12 = (double) lgen[dimension-2];
	c1 = c11/c12;
	c41 = (double)(-Ns)*(lgen[dimension-2]*lgen[dimension-1]*lgen[0]-lgen[dimension-2]*lgen[dimension-1]*lgen[1] + lgen[dimension-2]*lgen[0]*lgen[1]-lgen[dimension-1]*lgen[0]*lgen[1]);
	c42 = (double) lgen[dimension-2]*lgen[dimension-1]*lgen[1];
	c4= c41/c42;
	return max(ceil(c1),ceil(c4));	
}


///


long Lambda2(vector<long> lgen)
{
	long Ns,dimension;
	double c21, c22, c31,c32, c2, c3;
	Ns = ComputeNs(lgen);
	dimension = lgen.size();
	c21 = (double)(lgen[0]-lgen[1])*Ns;
	c22 = (double) lgen[1];
	c2 = c21/c22;
	c31 = (double) (Ns*(lgen[dimension-2]*lgen[dimension-1]*lgen[0]-lgen[dimension-2]*lgen[dimension-1]*lgen[1]+lgen[dimension-2]*lgen[0]*lgen[1]-lgen[dimension-1]*lgen[0]*lgen[1]));
	c32 = (double) lgen[dimension-2]*lgen[0]*lgen[1];
	c3 = c31/c32;
	return (long)(-min(c3,c2));
}


///


long Lambda2(vector<long> lgen, long Ns)
{
	long dimension;
	double c21, c22, c31,c32, c2, c3;
	dimension = lgen.size();
	c21 = (double)(lgen[0]-lgen[1])*Ns;
	c22 = (double) lgen[1];
	c2 = c21/c22;
	c31 = (double) (Ns*(lgen[dimension-2]*lgen[dimension-1]*lgen[0]-lgen[dimension-2]*lgen[dimension-1]*lgen[1]+lgen[dimension-2]*lgen[0]*lgen[1]-lgen[dimension-1]*lgen[0]*lgen[1]));
	c32 = (double) lgen[dimension-2]*lgen[0]*lgen[1];
	c3 = c31/c32;
	return (long)(-min(c3,c2));
}


///


long ComputeN0(std::vector<long> lgen)
{
	long dimension, Ns, a1, ap, lambda1, lambda2, N0;
	Ns = ComputeNs(lgen);
	dimension = lgen.size();
	a1 = lgen[0];
	ap = lgen[dimension-1];
	lambda1 = Lambda1(lgen, Ns);
	lambda2 = Lambda2(lgen, Ns);	
	N0 = max(Ns/a1, (ap-a1+lambda1+lambda2)/(ap-a1));
	return N0;
}


///


long ComputeN0(std::vector<long> lgen, long Ns)
{
	long dimension, a1, ap, lambda1, lambda2, N0;
	dimension = lgen.size();
	a1 = lgen[0];
	ap = lgen[dimension-1];
	lambda1 = Lambda1(lgen, Ns);
	lambda2 = Lambda2(lgen, Ns);
	double aux11, aux12, aux21, aux22;
	aux11 = (double) Ns;
	aux12 = (double) a1;
	aux21 = (double) (ap-a1+lambda1+lambda2);
	aux22 = (double) (ap-a1);
	N0 = max(ceil(aux11/aux12), ceil(aux21/aux22));
	return N0;
}


vector<long> ComputeDeltaNu(vector<long> lgen, long n)
{
	long Ns, N0, dimension, l1, l2, x1, x2;
	Ns = ComputeNs(lgen);
	N0 = ComputeN0(lgen,Ns);
	dimension = lgen.size();
	if(N0>n)
	{
		Pintar(Nu(lgen,n));
		Pintar(Delta(Nu(lgen,n)));
		return OrdenaSet(Delta(Nu(lgen,n)));
	}

	l1 = Lambda1(lgen,Ns);
	l2 = Lambda2(lgen,Ns);
	x1 = lgen[0]*n+l1;
	x2 = lgen[dimension-1]*n-l2;

	vector<long> longitudes1, longitudes2;
	long cotaB1, cotaB3;
	double num1, den1;
	num1 = (double) x1;
	den1 = (double) lgen[dimension-1];
	cotaB1 = ceil(num1/den1);
	cotaB3 = x2/lgen[0];

	vector<vector<long> > v;
	// Calculamos las longitudes del trozo 1
	for(long i=n*lgen[0]; i<x1+1; i++)
	{
		v = FSolve(lgen,i);
		vector<long> w;
		for(long j=0;j<(long) v.size();j++)
		{
			vector<long> factorizacion;
			factorizacion = v[j];
			long suma;
			suma = 0;
			for(long k=0;k<(long)factorizacion.size();k++)
			{
				suma +=factorizacion[k];
			}
			w.push_back(suma);
		}
		w = OrdenaSet(w);
		for(long j=0;j<(long) w.size();j++)
		{
			if(w[j]==n)
			{
				for(long k=0;k<(long) w.size();k++)
				{
					longitudes1.push_back(w[k]);
				}
				break;
			}
		}
	}
	longitudes1 = OrdenaSet(longitudes1);
	for(long i=0;i<(long)longitudes1.size();i++)
	{
		if(longitudes1[i]>cotaB1)
		{
			longitudes1.erase(longitudes1.begin()+i);
			i--;
		}
	}
	// Calculamos las longitudes para el trozo 2
	for(long i=x2-1; i<n*lgen[dimension-1]; i++)
	{
		v = FSolve(lgen,i);
		vector<long> w;
		for(long j=0;j<(long) v.size();j++)
		{
			vector<long> factorizacion;
			factorizacion = v[j];
			long suma;
			suma = 0;
			for(long k=0;k<(long)factorizacion.size();k++)
			{
				suma +=factorizacion[k];
			}
			w.push_back(suma);
		}
		w = OrdenaSet(w);
		for(long j=0;j<(long) w.size();j++)
		{
			if(w[j]==n)
			{
				for(long k=0;k<(long) w.size();k++)
				{
					longitudes2.push_back(w[k]);
				}
				break;
			}
		}
	}
	longitudes2 = OrdenaSet(longitudes2);
	for(long i=0;i<(long)longitudes2.size();i++)
	{
		if(longitudes2[i]<cotaB3)
		{
			longitudes2.erase(longitudes2.begin()+i);
			i--;
		}
	}
	longitudes1 = OrdenaSet(longitudes1);
	longitudes2 = OrdenaSet(longitudes2);
	vector<long> dif1, dif3;
	dif1 = Delta(longitudes1);
	dif3 = Delta(longitudes2);
	for(long i=0;i<(long)dif3.size();i++)
	{
		dif1.push_back(dif3[i]);
	}
	return OrdenaSet(dif1);
}

/// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// ///


void Pintar(vector<long> v)
{
	cout<<"*";
	for(unsigned ii=0;ii<v.size();ii++)
		cout<<","<<v[ii];
	cout<<"*";
	cout<<endl;
}

void Pintar(vector<double> v)
{
	cout<<"*";
	for(unsigned ii=0;ii<v.size();ii++)
		cout<<","<<v[ii];
	cout<<"*";
	cout<<endl;
}

///


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


///


long maximum(vector<long> v)
{
	long n;
	long aux;
	
	n = v.size();
	if( n == 0)
		return 0;
	aux = v[0];
	for(long i=1;i<n;i++)
	{
		if(aux < v[i])
		{
			aux = v[i];
		}
	}
	return aux;
}


///


double maximum(vector<double> v)
{
	long n;
	double aux;

	n = v.size();
	if( n == 0)
		return 0;
	aux = v[0];
	for(long i=1;i<n;i++)
	{
		if(aux < v[i])
		{
			aux = v[i];
		}
	}
	return aux;
}


///


vector<vector<long> > f1(long e, long n)
{
	vector<long> aux;
	for(long i=0;i<e;i++)
	{
		aux.push_back(1);
	}
	return FSolve(aux,n);
}


///


vector<long> W(vector<long> smg, long n)
{
	long dim;
	dim = smg.size();
	vector<vector<long> > xx;
	xx = f1(dim, n);
	long aux, suma;
	vector<long> laux, x;
	aux = xx.size();

	for(long i=0; i<aux; i++)
	{
		x = xx[i];
		suma = 0;
		for(long j=0; j<dim; j++)
		{
			suma += x[j]*smg[j]; 
		}
		laux.push_back(suma);
	}
	sort(laux.begin(),laux.end());

	// Delete duplicate
	aux = laux[0];
	for(unsigned i=1;i<laux.size();i++)
	{
		if(aux == laux[i])
		{
			laux.erase(laux.begin()+i);
			i=1;
		}
		aux=laux[i];
	}
	return laux;
}


///


vector<long> L(vector<long> lgen, long x)
{
	long suma;
	vector<vector<long> > l1;
	l1 = FSolve(lgen,x);
	vector<long> l2;
	for(unsigned i=0;i<l1.size();i++)
	{
		suma = 0;
		for(unsigned j=0;j<l1[i].size();j++)
		{
			suma += l1[i][j];
		}
		l2.push_back(suma);
	}
	sort(l2.begin(),l2.end());
	long aux;
	aux = l2[0];
	for(unsigned i=1;i<l2.size();i++)
	{
		if(aux == l2[i])
		{
			l2.erase(l2.begin()+i);
			i=1;
		}
		aux=l2[i];
	}
	return l2;
}


///


vector<long> Nu(vector<long> smg, long n)
{
	vector<long> waux, longaux, laux;
	waux = W(smg,n);
	for(long i=0;i<(long)waux.size();i++)
	{
		laux = L(smg,waux[i]);
		for(long j=0;j<(long)laux.size();j++)
		{
			longaux.push_back(laux[j]);
		}
	}
	
	long aux;
	sort(longaux.begin(),longaux.end());
	// Delete duplicate
	aux = longaux[0];
	for(unsigned i=1;i<longaux.size();i++)
	{
		if(aux == longaux[i])
		{
			longaux.erase(longaux.begin()+i);
			i=1;
		}
		aux=longaux[i];
	}
	return longaux;
}


///


vector<long> Delta(vector<long> laux)
{
	vector<long> l1;
	for(long i=1;i<(long)laux.size();i++)
	{
		l1.push_back(laux[i]-laux[i-1]);
	}
	long aux;
	sort(l1.begin(),l1.end());
	// Delete duplicate
	aux = l1[0];
	for(unsigned i=1;i<l1.size();i++)
	{
		if(aux == l1[i])
		{
			l1.erase(l1.begin()+i);
			i=1;
		}
		aux=l1[i];
	}
	return l1;
}


///


vector<long> OrdenaSet(vector<long> l1)
{
	long aux;
	sort(l1.begin(),l1.end());
	// Delete duplicate
	aux = l1[0];

	for(unsigned i=1;i<l1.size();i++)
	{
		if(aux == l1[i])
		{
			l1.erase(l1.begin()+i);
			i=0;
		}
		aux=l1[i];
	}
	return l1;
}
