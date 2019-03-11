def FrobeniusNumber(generators,eDimension=0):
    if eDimension == 0:
        eDimension = len(generators)
    lGen =[generators[i] for i in range(eDimension)]
    h=t=a=generators[0]
    b=lGen
    Q = [0 for i in range(a)]
    S = [a*lGen[-1] for i in range(a*lGen[-1])]
    S[a-1]=0
    P = [i for i in range(a)]
    P[a-1] = len(lGen)
    while(h != 0):
        QHaux = Q[h-1]
        v = h
        Q[h-1] = 0
        if(h == t):
            h = 0
        else:
            h = QHaux
        for j in range(1,P[v-1]+1):
            e = (b[j-1]+v) % a
            w = b[j-1] + S[v-1]
            if(w < S[e-1] and e != 0 ):
                S[e-1] = w
                P[e-1] = j
                if(Q[e-1] == 0):
                    if(h == 0 ):
                        h = e
                        Q[e-1] = e
                        t = e
                    else:
                        Q[t-1] = e
                        Q[e-1] = e
                        t = e
    aux = [S[i] for i in range(0,a)]
    fNumber = max(aux)-a
    return(fNumber)
