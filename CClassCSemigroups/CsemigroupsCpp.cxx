#include "CsemigroupsCpp.h" 

using namespace std;

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

///

vector<long> MultiplicityAxis(vector<vector<long>> generators, vector<long> ray)
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
    
    cout<<minElementIndex<<endl;
    
    return multiplicity2[minElementIndex];
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
    /*
    dim = len(pt)
    for x in eq:
        sum = 0
        for i in range(dim):
            sum = sum + pt[i]*x[i]
        sum = round(sum+x[-1],10)
        if sum > 0:
            return False
    return True
    */
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
        cout<<sum<< endl;
        if(sum > 0.0000000001) //Avoid "false zeros"
        {
            cout<<" * "<<sum<<" * "<< endl;
            return(false);
        }
    }
       return true;
}
    

//vector<vector<long>> diamondMultiplicity(vector<vector<long>> generators, vector<vector<long>> rays)
//{
    /*
    # Calculamos el diamante entero de las multiplicidades.
    # INPUT:
    #   - cone: rayos del cono.
    #   - smg: sistema minimal de generadores.
    # OUTPUT:
    #   - diamante entero de las multiplicidades
    */
    
//    vector<long> extremes;
//    int n;
//    n = cone.size();
//    for(int ii=0; ii<n; i++)
//    {
//        
//    }
    
    /*
        extremos = []
    for ray in cone:
        extremos.append(MultiplicityAxis(ray,smg))
    diamante = Diamond(extremos)
    hull = ConvexHull(diamante)
    eq = [list(x) for x in hull.equations]
    cotas = Cube(diamante)
    return IntegerDiamond(eq,cotas)
    */
//}


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