#include "CsemigroupsCpp.h" 

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


int foo(int a)
{
    return a+1;
}

void foo2(std::vector<std::vector<double>> gen)
{
       cout<<"Funciona"<<endl;
}



void Pintar(vector<long> v)
{
	cout<<"*";
	for(unsigned ii=0;ii<v.size();ii++)
		cout<<","<<v[ii];
	cout<<"*";
	cout<<endl;
}

void Pintar(vector<vector<long>> m)
{
    unsigned n;
    cout<<"**";
    for(unsigned ii=0;ii<m.size();ii++)
    {
        cout<<"++";
        n = m[ii].size();
        for(unsigned jj=0;jj<n;jj++)
        {
            cout<<","<<m[ii][jj];
        }
        cout<<"++"<<endl;
    }
	cout<<"**";
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

vector<long> operator+(const vector<long>& v1, const vector<long>& v2)
{
    unsigned n;
    n = v1.size();
    vector<long> aux(n);
    for(unsigned ii=0;ii<n;ii++)
    {
        aux[ii] = v1[ii]+v2[ii];
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

long prodEsc(const std::vector<long>& v1, const std::vector<long>& v2)
{
    long suma;
    int n;
    suma = 0;
    n = v1.size();
    for(unsigned ii=0;ii<n;ii++)
    {
           suma += v1[ii]*v2[ii];
    }
    return suma;
}

double prodEsc(const std::vector<long>& v1, const std::vector<double>& v2)
{
    long suma;
    int n;
    suma = 0;
    n = v1.size();
    for(unsigned ii=0;ii<n;ii++)
    {
           suma += v1[ii]*v2[ii];
    }
    return suma;
}

vector<vector<long>> deleteDuplicates(vector<vector<long>> m)
{
    vector<vector<long>> aux;
    unsigned mSize, auxSize;
    bool isIn;
    mSize = m.size();
    auxSize = 1;
    aux.push_back(m[0]);
    for(unsigned ii=1;ii<mSize;ii++)
    {
        isIn = false;
        for(unsigned jj=0;jj<auxSize;jj++)
        {
            if(m[ii] == aux[jj])
            {
                isIn = true;
            }
        }
        if(!isIn)
        {
            aux.push_back(m[ii]);
            auxSize++;
        }
    }
    return aux;
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
            if(x[ii]%r[ii] != 0)
            {
                return 0;
            }
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

///

vector<long> multiplicityAxis(vector<vector<long>> generators, vector<long> ray)
{
    /*
    # Calculamos la multiplicidad en cada rayo.
    # INPUT:
    #   - ray: rayo del cono.
    #   - smg: sistema minimal de generadores.
    # OUTPUT:
    #   - Elemento minimal en el rayo.
    */
    vector<long> multiplicity1;
    vector<vector<long>> multiplicity2;
    long n;
    n = generators.size();
    for(unsigned ii=0; ii<n;ii++)
    {
        long aux;
        aux = belongAxis(generators[ii],ray);
        if(aux != 0)
        {
            multiplicity1.push_back(aux);
            multiplicity2.push_back(generators[ii]);
        }
    }
    
    int minElementIndex = min_element(multiplicity1.begin(),multiplicity1.end()) - multiplicity1.begin(); // GetMinimumElement
    
    
    return multiplicity2[minElementIndex];
}


vector<vector<long>> multiplicityAllAxes(vector<vector<long>> generators, vector<vector<long>> rays)
{
    vector<vector<long>> multiplicities;
    unsigned nRays;
    nRays = rays.size();
    for(unsigned ii=0;ii<nRays;ii++)
    {
        multiplicities.push_back(multiplicityAxis(generators,rays[ii]));
    }
    return multiplicities;
}


vector<vector<long>> diamond(vector<vector<long>> mult)
{
    /*
    # En primer lugar tenemos que calcular el diamante.
    # INPUT:
    #   - a: Minimal elements of the ray.
    # OUTPUT:
    #   - Points of the diamond
    */
    
    unsigned n = mult.size();
    unsigned m = mult[0].size();
    vector<long> vZero(m, 0);
    vector<vector<long>> aux(n);
    for(unsigned ii=0;ii<n;ii++)
    {
            aux[ii]=mult[ii];
    }
    aux.push_back(vZero);
    
    for(unsigned ii=0;ii<n;ii++)
    {
        for(unsigned jj=ii+1;jj<n;jj++)
        {
            aux.push_back(mult[ii]+mult[jj]);
        }
    }
    
    return aux;
}

bool pointBelongsDiamond(vector<long> pt, vector<vector<double>> eq)
{

    unsigned dim,nEq;
    dim = pt.size();
    nEq = eq.size();
    double sum;
    for(unsigned ii=0;ii<nEq;ii++)
    {
        sum = 0;
        for(unsigned jj=0;jj<dim;jj++)
        {
            sum += pt[jj]*eq[ii][jj];
        }
        sum += eq[ii][dim];
        if(sum > 0.0000000001) //Avoid "false zeros"
        {
            return(false);
        }
    }
       return true;
}

// This function check which points of points verify equations eq.
vector<vector<long>> filterPoints(vector<vector<long>> points, vector<vector<double>> eq)
{
    vector<vector<long>> filter;
    int n;
    n = points.size();
    for(unsigned ii=0;ii<n;ii++)
    {
        if(pointBelongsDiamond(points[ii], eq))
        {
         filter.push_back(points[ii]);   
        }
    }
    return(filter);
}

vector<vector<long>> eqRay(vector<long> ray, std::vector<std::vector<long>> hyperplanes)
{
    vector<vector<long>> eq;
    int n;
    n = hyperplanes.size();
    for(unsigned ii=0;ii<n;ii++)
    {
        if(prodEsc(ray,hyperplanes[ii])==0)
        {
            eq.push_back(hyperplanes[ii]);
        }
    }
    return(eq);
}

vector<vector<long>> deleteRowZero(vector<vector<long>> m)
{
    vector<vector<long>> aux;
    int n1,n2;
    bool allzero;
    n1 = m.size();
    for(unsigned ii=0;ii<n1;ii++)
    {
        allzero = true;
        n2 = m[ii].size();
        for(unsigned jj=0;jj<n2;jj++)
        {
            if(m[ii][jj] != 0)
            {
                allzero = false;
            }
        }
        if(!allzero)
        {
            aux.push_back(m[ii]);
        }
    }
    return(aux);
}

vector<vector<long>> affineTerm(vector<vector<long>> eq, vector<vector<long>> d)
{

    int neq, dsize;
    neq = eq.size();
    //dim = eq[0].size();
    dsize = d.size();
    vector<vector<long>> afin;
    for(unsigned ii=0;ii<dsize;ii++)
    {
        vector<long> aux;
        for(unsigned jj=0;jj<neq;jj++)
        {
            aux.push_back(prodEsc(d[ii],eq[jj]));
        }
        afin.push_back(aux);
    }
    return(deleteRowZero(afin));
}

bool existGenerator(vector<vector<long>> equationsRay, vector<long> affine, vector<vector<long>> generators)
{
    /*
        # Calculamos si en cada recta afín paralela al rayo hay un generador.
        # INPUT:
        #   - eqray: ecuaciones de un rayo.
        #   - afinset: valor afín de la recta.
        #   - smg: generadores del semigrupo.
        # OUTPUT:
        #   True/False si en esa recta hay un generador.
    def ExistGenerator(eqray, afin,smg):
        aux = AffineTerm(eqray,smg)
        if afin in aux:
            return True
        else:
            return False
    */
    vector<vector<long>> aux;
    aux = affineTerm(equationsRay,generators);
    int n;
    n = aux.size();
    for(unsigned ii=0;ii<n;ii++)
    {
        if(affine == aux[ii])
            return(true);
    }        
    return(false);
}

bool studyRays(vector<vector<long>> rays, vector<vector<long>> hyperplanes, vector<vector<long>> integerDiamond, vector<vector<long>> generators)
{
    /*
        for ray in rayos:
            # Calculamos las ecuaciones del rayo.
            eqrayo = EqRay(ray,hp)
            # Calculamos los términos afines de los puntos enteros del diamantes.
            afinesDiamante = AffineTerm(eqrayo,diamanteEntero)
            # Nos quedamos con los generadores afines.
            genAfines = AffineSemigroup(afinesDiamante, "generators").getMSG()
            # Comprobamos que por cada afín pasa un generador
            for afin in genAfines:
                if not ExistGenerator(eqrayo, afin,smg):
                    return False
        return True
    */
    int nRays;
    nRays = rays.size();
    for(unsigned ii=0;ii<nRays;ii++)
    {
        Pintar(rays[ii]);
        vector<vector<long>> equationsRay;
        equationsRay = eqRay(rays[ii],hyperplanes);
        vector<vector<long>> afDiamond, genAf;
        afDiamond = affineTerm(equationsRay,integerDiamond);
        genAf = computeMSG(afDiamond);
        int numGen;
        numGen = genAf.size();
        for(unsigned jj=0; jj<numGen; jj++)
        {
            if(!existGenerator(equationsRay,genAf[jj],generators))
            {
                return(false);
            }
        }        
    }
    return(true);
}

vector<long> conductorAxis(vector<vector<long>> generators, vector<long> ray)
{
    /*
    
    # Calculamos el conductor.
    # INPUT:
        #   - gen: Set of generators.
    #   - r: Minimal value in the ray.
    # OUTPUT:
    #   - conductor
    */
    vector<long> aux, aux2, conductor;
    unsigned nGen, nRay;
    long candidate, frobNumber;
    nGen = generators.size();
    nRay = ray.size();
    
    for(unsigned ii=0;ii<nGen;ii++)
    {
        candidate = belongAxis(generators[ii],ray);
        if(candidate != 0)
        {
            aux.push_back(candidate);
        } 
    }
    
    frobNumber = FrobeniusNumber(aux);
    
    for(unsigned ii=0;ii<nRay;ii++)
    {
        conductor.push_back(ray[ii]*(frobNumber+1));
    }
    return conductor;
}

vector<vector<long>> computeXDiamond(vector<vector<long>> generators, vector<vector<long>> rays, vector<vector<long>> equations,vector<vector<long>> diamondMultiplicity)
{
    unsigned nRays, nAffineDiamond;
    long maxMinNorm;
    vector<long> conductor;
    vector<vector<long>> eqRays, affineDiamond;
    nRays = rays.size();
    for(unsigned ii=0;ii<nRays;ii++)
    {
        conductor = conductorAxis(generators,rays[ii]);
        eqRays = eqRay(rays[ii],equations);
        affineDiamond = deleteDuplicates(affineTerm(eqRays, diamondMultiplicity));
        maxMinNorm = 0;
        nAffineDiamond = affineDiamond.size();
        for(unsigned jj=0;jj<nAffineDiamond;jj++)
        {
            
        }
    }
    return eqRays;
    /*
    for ray in rayos:
        # Calculamos el conductor.
        conductor = ConductorAxis(gen,ray)
        # Calculamos las ecuaciones del rayo.
        eq = EqRay(ray,hp)
        # Calculamos los términos afines.
        afinesDiamante = DeleteDuplicates(AffineTerm(eq,diamanteM))
        # Calculamos la máxima norma mínima en cada recta afin.
        maxMinNorm = 0
        for afin in afinesDiamante:
            # En primer lugar miramos la norma minima en el diamante.
            minimumNormDiamond = MinimumPointAffineDiamond(diamanteM,afin,eq)
            # Ahora buscamos el elemento mínimo en el diamante en dicha recta.
            minimumNormSG = ComputeMinimumNormInSemigroupLine(ray,minimumNormDiamond,gen)
            if maxMinNorm < minimumNormSG:
                maxMinNorm = minimumNormSG
        n = NormOne(conductor) + maxMinNorm
        x.append(ComputeBoundRay(n,ray,gen))
    # Calculamos el diamante de esos valores X
    diamanteX = Diamond(x)
    */
}


//vector<vector<long>> computeGaps(vector<vector<long>> generators, vector<vector<long>> rays,vector<vector<long>> hyperplanes)
//{
       /*
       Compute the set of gaps of a CSemigroup
       INPUT: Generators, rays and hyperplanes
       OUTPUT: Set of gaps
       */
    
    // Compute the diamonds of multiplicities.
    
    
        /*
    # Calculamos del diamante de multiplicidades.
    diamanteM = DiamondMultiplicity(rayos,gen)
    x = []
    # Para cada rayo del cono.
    for ray in rayos:
        # Calculamos el conductor.
        conductor = ConductorAxis(gen,ray)
        # Calculamos las ecuaciones del rayo.
        eq = EqRay(ray,hp)
        # Calculamos los términos afines.
        afinesDiamante = DeleteDuplicates(AffineTerm(eq,diamanteM))
        # Calculamos la máxima norma mínima en cada recta afin.
        maxMinNorm = 0
        for afin in afinesDiamante:
            # En primer lugar miramos la norma minima en el diamante.
            minimumNormDiamond = MinimumPointAffineDiamond(diamanteM,afin,eq)
            # Ahora buscamos el elemento mínimo en el diamante en dicha recta.
            minimumNormSG = ComputeMinimumNormInSemigroupLine(ray,minimumNormDiamond,gen)
            if maxMinNorm < minimumNormSG:
                maxMinNorm = minimumNormSG
        n = NormOne(conductor) + maxMinNorm
        x.append(ComputeBoundRay(n,ray,gen))
    # Calculamos el diamante de esos valores X
    diamanteX = Diamond(x)
    # Creamos un cubo que contenga al diamante.
    cotaX = Cube(diamanteX)
    # Calculamos las ecuaciones del diamante.
    hullX = ConvexHull(diamanteX)
    eqDiamanteX = [list(x) for x in hullX.equations]
    # Calculamos el diamante entero X.
    diamanteEnteroX = IntegerDiamond(eqDiamanteX,cotaX)
    # Inicialmente el conjunto de huecos es vacío.
    gapset = []
    # Definimos el semigrupo afín.
    afseg = AffineSemigroup(gen, "generators")
    for x in diamanteEnteroX:
        if not afseg.belongs(x):
            gapset.append(x)
    return gapset
        */
//}