def pBounds(p,i,k,e):
    # For fixed p,i,e, is this function monotonic in k?
    ordersRight = [0]*(i*k)
    orders = [0]*(i*k)
    for c in range(i*k-1,0,-1):
        oTop=0
        oRight=0
        if c+e<i*k:
            oTop=ordersRight[c+e]
        if p*c<i*k:
            oRight=ordersRight[p*c]
        ordersRight[c]=max(1+oTop,oRight + i - (p*c)//k)
        orders[c] = ordersRight[c]+i-1-c//k
    return orders

def matrix_precision(M):
    return min([m.precision_absolute() for m in M.coefficients()])

def prismatic_envelope_f(p,E,k,prec,zprec,fprec):
    # Returns B,grN,phi,delta, where B is a polynomial ring in f0,...,fprec-1,
    # grN is a function which measures the Nygaard filtration of an element of
    # B, phi is the graded lift of Frobenius, and delta is the graded
    # delta-ring structure.
    
    # F-weight z is 1, Nygaard of z is 0.
    # F-weight of f_j is k*p^j, Nygaard weight of f_j is p^j.
    
    B=PolynomialRing(A,'f',fprec)
    B.inject_variables()
    variable_names=B.variable_names()
    # A list of the fi for ease of reference.
    fvars=[]
    for j in variable_names:
        fvars.append(B(j))
    fvars_phitilde=[]
    for j in range(len(fvars)-1):
        fvars_phitilde.append(fvars[j]^p+p*fvars[j+1])
    # At the last stop, the delta term is beyond the f-precision
    # so it is not added.
    fvars_phitilde.append(fvars[len(fvars)-1]^p)
    
    def weight(funct):
        # Returns a tuple (F-weight,N-weight), where F-weight is the F-weight up to zprec
        # and N-weight is a *lower bound* on the Nygaard weight of f.
        # For example, this function will not correctly compute the Nygaard weight of large powers of p.
        F_weight=zprec
        # The highest possible Nygaard weight is that of z^{k-1}f_{fprec-1}, which has Nygaard weight p^{fprec-1}.
        nygaard_weight=p^(fprec-1)
        p_powers=[p^j for j in range(fprec)]
        
        def monomial_weight(m):
            # Returns a tuple (F-weight,N-weight), where F-weight is the exact F-weight
            # and N-weight is the exact Nygaard weight of a monomial in the f_i.
            monomial_nygaard_weight=sum([a*b for a,b in zip(p_powers,list(m.degrees()))])
            return (k*monomial_nygaard_weight,monomial_nygaard_weight)
        
        for m in funct.monomials():
            f_weight,w=monomial_weight(m)
            # c is a power series in z.
            c=A(funct.monomial_coefficient(m))
            # The z-adic valuation.
            v=c.valuation()
            if v>0:
                vk=floor(log(v,k))
            else:
                vk=0
            nygaard_weight=min(nygaard_weight,vk+w)
            F_weight=min(F_weight,f_weight+v)
        
        return (F_weight,nygaard_weight)
    
    def phitilde(funct):
        g=0
        for m in funct.monomials():
            c=funct.monomial_coefficient(m)
            g+=c.V(p)*(m(fvars_phitilde))
        return g
                
    def deltatilde(funct):
        return (B(phitilde(funct)-funct^p))*(1/p)
    
    def coefficient_divide_f(funct,gunct,zprec,weight):
        # Takes a polynomial gunct in the fj and divides all the
        # coefficients by funct.
        out=0
        for m in gunct.monomials():
            (f_weight,nygaard_weight)=weight(m)
            new_coefficient=lazy_division(funct,gunct.monomial_coefficient(m),zprec-f_weight)
            out+=new_coefficient*m
        return out
    
    def trim(funct):
        # Here, we use the fact that we have a fixed z-precision and the f_j secretly
        # have z-weight kp^j. So, we can trim lots of information.
        g=0
        for m in funct.monomials():
            w=weight(m)[0]
            new_coefficient=funct.monomial_coefficient(m).add_bigoh(max(zprec-w,0))
            g+=new_coefficient*m
        return B(g)
    
    rels=[]
    if fprec==1:
        rels.append(fvars[0]^p)
    else:
        rels.append(trim(phitilde(B(E))*f1 + deltatilde(B(E))*f0^p))
        for j in range(1,fprec-1):
            rels.append(trim(deltatilde(rels[j-1])))
        rels.append(fvars[fprec-1]^p)
        
    # new_rels is a lookup table for how to rewrite f_j^p in terms of other factors.
    new_rels=[]
    for j in range(len(rels)):
        r=B(rels[j])
        u=r.monomial_coefficient(fvars[j]^p)
        new_rels.append(trim(fvars[j]^p)-(1/u)*r)
 
    def reduce(funct):
        # Run through the monomials of f and rewrite them once using new_rels.
        reduced=0
        for m in funct.monomials():
            new_monomial=1
            degs=m.degrees()
            for exponent in range(len(degs)):
                r=degs[exponent]//p
                new_monomial=new_monomial*new_rels[exponent]^r*fvars[exponent]^(degs[exponent]-p*r)
            reduced+=funct.monomial_coefficient(m)*new_monomial
        return trim(B(reduced))
                
    # Recursively rewrite the f_j^p -> ... dictionary, obtaining reduced forms for those relations.
    #note that the relations we obtain from this are still independent of k.
    for j in range(len(new_rels)):
        prev_value=new_rels[j]
        new_rels[j]=reduce(B(prev_value))
        while prev_value != new_rels[j]:
            prev_value=new_rels[j]
            new_rels[j]=reduce(B(prev_value))
    
    # Contains for each f_j, the value of phi^{{p^j}} = (f_j^p + p f_{j+1})/(E(z)^{p^{j+1}}).
    fvars_phi_divided=[]
    for j in range(len(fvars)-1):
        phitilde_reduced = new_rels[j]+p*fvars[j+1]
        fvars_phi_divided.append(coefficient_divide_f(A(E)^(p^(j+1)),phitilde_reduced,zprec,weight))
    fvars_phi_divided.append(B(0))
    
    def reduce_zk(funct):
        reduced=0
        def coefficient_reduce(c):
            # Takes a coefficient (a power series in z)
            # and returns an element of B obtained by using z^k=f_0
            out=B(c[0])
            for j in range(1,min(zprec,c.precision_absolute())):
                r=j//k
                out+=fvars[0]^r*z^(j-k*r)*c[j]
            return B(out)
        for m in funct.monomials():
             reduced+=coefficient_reduce(funct.monomial_coefficient(m))*m
        return trim(B(reduced))
    
    def recursive_reduce(funct):
        prev_value = funct
        funct = reduce_zk(B(funct))
        funct = reduce(B(funct))
        while prev_value!=funct:
            prev_value = funct
            funct = reduce_zk(B(funct))
            funct = reduce(B(funct))
        return funct

    # Only works for *homogeneous* terms in the f_j, and they get divided by the obvious Nygaard filtration.
    def phi_divided(funct):
        result=0
        for m in funct.monomials():
            c=funct.monomial_coefficient(m)
            result+=B(c.V(p)*m(fvars_phi_divided))
        return result

    return B,fvars,weight,phitilde,phi_divided,deltatilde,reduce,recursive_reduce
    
def nablaP_matrix_OK(p,i,k,E,prec,zprec,gprec):
    # Returns nablaP_OK.
    # Note that this does not technically depend on k, so one can run this once
    # for a large k and slice down to compute for several k at once.
    C,gvarsC,weightC,phiC,deltaC,phiCtilde,deltaCtilde,reduceC,recreduceC,g0_divideC,coefficient_divideC=prismatic_envelope_g(p,E,k,prec,zprec,gprec)
    
    def initialize_u():
        # Takes almost no time.
        if len(gvarsC)==1:
            divided_g0=gvarsC[0]^p
        else:
            divided_g0=(gvarsC[0]^p+p*gvarsC[1])
        divided_g0=recreduceC(divided_g0)
        divided_g0=coefficient_divideC(A(E)^p,divided_g0,zprec,weightC)

        divided_E=C(E(z)-E(z-gvarsC[0]))
        divided_E=g0_divideC(divided_E)
        divided_E=recreduceC(divided_E)

        u=C(1-phiC(divided_E)*divided_g0)
        u=recreduceC(u)
        return u
    
    u=initialize_u()

    def recursive_phiC_product(u):
        v=u
        while True:
            funct=recreduceC(u*phiC(v))
            if funct==v:
                return v
            else:
                v=funct
                                                
    v=recursive_phiC_product(u)
    
    def initialize_bk_factor():
        # Standard powering algorithm leads to substantial speedups.
        result = C(1)
        base = v
        exp = i

        while(exp >= 2):
            if(exp % 2 == 1):
                 result = recreduceC(result*base)
            base = recreduceC(base*base)
            exp = exp // 2
        result = recreduceC(result*base)
        return result
    
    bk_factor=initialize_bk_factor()
    
    # nablaP builder
    # First, nablaP_OK and then a nablaP_OKpik. We get the second using a saturation
    # method using the maps from prismaOK to prism OKpiK.

    # nablaP_OK builder
    # Element of BK twist i is x*s^i, x in Prisma and s is the BK orientation.
    # So, x is a linear combination of powers of z.
    # Have to compute (\eta_L-\eta_R)(z^j*s^i).
    # Split into two things: \eta_L is easy.
    # What matters is coefficient of g_0 in -\eta_R(z^j*s^i) after reduction.
    # Warning: note the sign!
    # eta is multiplicative: so -\eta_R(z^j*s^i)=-\eta_R(z)^j*\eta_R(s)^i.
    # \eta_R(z)=z-g0.
    # \eta_R(s)=v*s,
    # v=BK twist factor=\prod_r\phi^r(u), r>=0.
    # u=\phi(E(z_R))/\phi(E(z_L))=some polynomial in the z_R,z_L and divided F on g0.
    # Interested in coefficient in of g0 in -(z-g0)^j*v^i

    # Figure out r bound.
    # Figure gprec.
    # Can probably get away with floor(log(k*i-1,p)).

    def initialize_nablaP_OK():
        nablaP_OK=Matrix(W,k*i-1,k*i-1)
        old_column_to_process=-bk_factor
        for n in range(1,i*k):
            column_to_process=recreduceC((z-gvarsC[0])*old_column_to_process)
            old_column_to_process=column_to_process
            input_build=[1]+[0 for i in range(len(gvarsC)-1)]
            send_g0_to_one_dict={gvarsC[0]:1}
            for j in range(1,len(gvarsC)):
                send_g0_to_one_dict[gvarsC[j]]=0
            send_g0_to_one=A(column_to_process.coefficient(send_g0_to_one_dict))
            for m in range(0,i*k-1):
                nablaP_OK[m,n-1]=send_g0_to_one[m]
        return nablaP_OK
            
    return initialize_nablaP_OK()

def syntomic_matrices(p,i,k,E,prec,zprec,fprec,gprec,nablaP_OK=False):
    # Returns syn0,syn1,nablaN,nablaP
    # The slowest part of the computation is the computation of the preliminary
    # matrix nablaP_OK, so this can alternatively be passed as an argument
    # in the case it is precomputed.
    B,fvars,weightB,phiBtilde,phi_dividedB,deltaBtilde,reduceB,recreduceB=prismatic_envelope_f(p,E,k,prec,zprec,fprec)

    # Basis z^c\prod_{j=0}^{fprec-1}f_j^{a_j}, 0<=a_j<=p-1, 0<=c<=k-1 with one exception: no 1.

    # can0 builder
    can0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=1
        for j in range(fprec):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(E(z)^(i-(n//k))*z^c*fprod)
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=1
            for j in range(fprec):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            can0[m-1,n-1]=coefficient_to_process[a]
    #print('can0 is {}'.format(can0))

    # phi0 builder
    phi0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=1
        for j in range(fprec):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(z^(p*c)*phi_dividedB(fprod))
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=1
            for j in range(fprec):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            phi0[m-1,n-1]=coefficient_to_process[a]

    # can1 builder
    # Same as above but with basis starting with 1 and going up to z^{k-2}*fprec or z^{k-1}f_{prec-1}
    # Does not make great sense for k=1.
    can1=Matrix(W,k*i-1,k*i-1)
    for n in range(0,k*i-1):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=1
        for j in range(fprec):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(E(z)^(i-(n//k)-1)*z^c*fprod)
        for m in range(0,k*i-1):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=1
            for j in range(fprec):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            can1[m,n]=coefficient_to_process[a]
    #print('can1 is {}'.format(can1))
            
    # nablaP builder
    # First, nablaP_OK and then a nablaP_OKpik. We get the second using a saturation
    # method using the maps from prismaOK to prism OKpiK.
    if nablaP_OK==False:
        nablaP_OK=nablaP_matrix_OK(p,i,k,E,prec,zprec,gprec)
    
    OKtoOKmodpi0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        column_to_process=recreduceB(B(z^n))
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=1
            for j in range(fprec):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            OKtoOKmodpi0[m-1,n-1]=coefficient_to_process[a]
    #print('OKtoOKmodpi0 is {}'.format(OKtoOKmodpi0))

    OKtoOKmodpi1=Matrix(W,k*i-1,k*i-1)
    for n in range(0,k*i-1):
        column_to_process=recreduceB(B(z^n))
        for m in range(0,k*i-1):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=1
            for j in range(fprec):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            OKtoOKmodpi1[m,n]=coefficient_to_process[a]
    #print('OKtoOKmodpi1 is {}'.format(OKtoOKmodpi1))
    
    success=0

    precision_loss1=max([valuation(x,p) for x in OKtoOKmodpi0.elementary_divisors()])
    try:
        OKtoOKmodpi0,OKtoOKmodpi1,nablaP_OK,nablaP=square_complete_bottom(OKtoOKmodpi0,OKtoOKmodpi1,nablaP_OK)
    except ValueError:
        print('Error defining nablaP, increasing the precision by {}.'.format(2*precision_loss1))
        return 0,0,0,0,success,[precision_loss1,0,0]

    syn0=can0-phi0
    precision_loss2=max([valuation(x,p) for x in can1.elementary_divisors()])
    try:
        can0,can1,nablaN,nablaP=square_complete_top(can0,can1,nablaP)
    except ValueError:
        print('Error defining can1, upping the precision by {}.'.format(2*(precision_loss2)))
        return 0,0,0,0,success,[precision_loss1,precision_loss2,0]

    precision_loss3=max([valuation(x,p) for x in nablaN.elementary_divisors()])
    try:
        syn0,syn1,nablaN,nablaP=square_complete_right(syn0,nablaN,nablaP)
    except ValueError:
        print('Error defining syn1, upping the precision by {}.'.format(2*(precision_loss3)))
        print('nablaN precision is {}'.format(matrix_precision(nablaN)))
        print('nablaP precision is {}'.format(matrix_precision(nablaP)))
        print('nablaP_OK precision is {}'.format(matrix_precision(nablaP_OK)))
        print('syn0 precision is {}'.format(matrix_precision(syn0)))
        print('can0 precision is {}'.format(matrix_precision(can0)))
        print('phi0 precision is {}'.format(matrix_precision(phi0)))
        print('can1 precision is {}'.format(matrix_precision(can1)))
        print('OKtoOKmodpi0 precision is {}'.format(matrix_precision(OKtoOKmodpi0)))
        print('OKtoOKmodpi1 precision is {}'.format(matrix_precision(OKtoOKmodpi1)))

        return 0,0,0,0,success,[precision_loss1,precision_loss2,precision_loss3]

    success=1
    return syn0,syn1,nablaN,nablaP,success,[precision_loss1,precision_loss2,precision_loss3]
    
def standard_projection_matrix(a,b):
    if a>b:
        raise TypeError('The input should be a pair (a,b) where a<=b.')
    else:
        M=identity_matrix(b)
        return M[0:a,:]

def saturation(M):
    # Input: an r-by-s matrix M, which can be defined
    # over any ring which supports smith normal form.
    # Output: matrices Mtilde, P, Q, where
    # Mtilde is the saturation of M,
    # P is the rank-by-r projection operator P such that P*M=Mtilde,
    # Q is the inclusion operator such that Q*Mtilde=M.
    R=M.base_ring()
    r=rank(M)
    target_dim=M.dimensions()[0]
    if r==target_dim:
        return M,Matrix(R,diagonal_matrix([1 for i in range(r)])),Matrix(R,diagonal_matrix([1 for i in range(r)]))
    else:
        D,S,T=M.smith_form()
        P=standard_projection_matrix(r,target_dim)
        Dtilde=P*D
        return Matrix(R,Dtilde*(T^(-1))),Matrix(R,P*S),Matrix(R,S^(-1)*P.transpose())
    
def saturate_square(syn0,syn1,etaN,etaP):
    # Takes four matrices such that syn1*etaN=etaP*syn0,
    # and returns four matrices syn0,syn1tilde,etaNtilde,etaPtilde,
    # where etaNtilde and etaPtilde are the saturations of etaN and etaP,
    # and where syn1tilde is the induced map between the saturations
    # so that syn1tilde*etaNtilde=etaPtilde&syn0.
    etaNtilde,etaN_P,etaN_Q=saturation(etaN)
    etaPtilde,etaP_P,etaP_Q=saturation(etaP)
    syn1tilde=etaP_P*syn1*etaN_Q
    return syn0,syn1tilde,etaNtilde,etaPtilde

def homology_total(syn0,syn1,etaN,etaP):
    # Out of laziness, the current implementation assumes
    # the rational homology of the total complex is 0.
    # This means that H^2 can be computed using elementary divisors of d^1
    # and then H^0 and H^1 can be computed using the elementary divisors
    # of the saturation of d^0.
    # Returns a dictionary 'hi':(rank,[list of torsion orders])
    null=Matrix(ZZ,0,0)
    if syn0==null or syn1==null or etaN==null or etaP==null:
        raise NotImplementedError('Funcionality is only implemented if the matrices are nonempty. Probably they need positive rank too.')
        return False
    syn0,syn1,etaN,etaP=saturate_square(syn0,syn1,etaN,etaP)
    d0=block_matrix([[syn0],[etaN]])
    d1=block_matrix([[etaP,-syn1]])
    dict={'h0':[],'h1':[],'h2':[]}
    el_divs=d1.elementary_divisors()
    r=0
    torsion=[]
    for d in el_divs:
        if d==0:
            r+=1
        elif d==1:
            pass
        else:
            torsion.append(d)
    dict['h2']=(r,torsion)
    d0tilde,P,Q=saturation(d0)
    el_divs=d0tilde.elementary_divisors()
    r1=0
    torsion=[]
    for d in el_divs:
        if d==0:
            r1+=1
        elif d==1:
            pass
        else:
            torsion.append(d)
    dict['h1']=(r1,torsion)
    dict['h0']=(0,[])
    return dict

def matrix_exponent(M):
    # Input: a square matrix of full rank.
    # Output: the smallest exponent N such that if y is in
    # the target space, Ny is in the image.
    # Warning: this does not check that M has full rank.
    el_divs=M.elementary_divisors()
    N=1
    for d in el_divs:
        N=lcm(N,d)
    return N

def square_complete_top(can0,can1,nablaP):
    # Completes the top of a partially commutative diagram of
    # matrices of full rank.
    R=nablaP.base_ring()
    return can0,can1,Matrix(R,can1^(-1)*nablaP*can0),nablaP

def square_complete_right(syn0,nablaN,nablaP):
    R=nablaP.base_ring()
    return syn0,Matrix(R,nablaP*syn0*nablaN^(-1)),nablaN,nablaP

def square_complete_bottom(OKtoOKmodpi0,OKtoOKmodpi1,nablaP_OK):
    # Completes the bottom of a partially commutative diagram of
    # matrices of full rank.
    R=nablaP_OK.base_ring()
    return OKtoOKmodpi0,OKtoOKmodpi1,nablaP_OK,Matrix(R,OKtoOKmodpi1*nablaP_OK*OKtoOKmodpi0^(-1))

def syntomic_cohomology(syn0,syn1,nablaN,nablaP):
    # Out of laziness, the current implementation assumes
    # the rational homology of the total complex is 0.
    # This means that H^2 can be computed using elementary divisors of d^1
    # and then H^0 and H^1 can be computed using the elementary divisors
    # of the saturation of d^0.
    # Returns a dictionary 'hi':(rank,[list of torsion orders])
    
    # This is really only here because padic matrices
    # don't have ranks :(
    d0=block_matrix([[syn0],[nablaN]])
    d1=block_matrix([[nablaP,-syn1]])
    dict={'h0':[],'h1':[],'h2':[]}
    el_divs=d1.elementary_divisors()
    r=0
    torsion=[]
    for d in el_divs:
        if d==0:
            r+=1
        elif d==1:
            pass
        else:
            torsion.append(d)
    dict['h2']=(r,torsion)
    el_divs=d0.elementary_divisors()
    r1=0
    torsion=[]
    for d in el_divs:
        if d==0:
            r1+=1
        elif d==1:
            pass
        else:
            torsion.append(d)
    dict['h1']=(r1,torsion)
    dict['h0']=(0,[])
    return dict

def lazy_division(funct,gunct,zprec):
    # Returns quotient + O(z^zprec) such that quotient*funct=gunct mod z^zprec.
    a=funct.valuation()
    b=gunct.valuation()
    quotient=0
    new_gunct=gunct
    for j in range(b-a,zprec):
        try:
            # Possible p-adic precision loss here.
            quotient+=z^j*W(new_gunct[j+a]/funct[a])
        except ValueError:
            raise ValueError('Function g={gstr} not divisible by f={fstr} up to precision={hstr}.'.format(gstr=gunct, fstr=funct, hstr=zprec))
        new_gunct=gunct-funct*quotient
    return quotient

def prismatic_envelope_g(p,E,k,prec,zprec,gprec):
    # ToDo: flatten away d before returning.
    # For working with O_K.
    # Returns C,grN,phi,delta, where C is a polynomial ring in g0,...,gprec-1,
    # grN is a function which measures the Nygaard filtration of an element of
    # C, phi is the graded lift of Frobenius, and delta is the graded
    # delta-ring structure.
    
    # z has F-weight 1 and Nygaard-weight 0.
    # d has F-weight 0 and Nygaard-weight 1.
    # gj has F-weight p^j and Nygaard-weight p^j.
    A0.<d>=PolynomialRing(A)
    A0.inject_variables()
    C=PolynomialRing(A0,'g',gprec)
    C.inject_variables()
    variable_names=C.variable_names()
    # A list of the gi for ease of reference.
    gvars=[]
    for j in variable_names:
        gvars.append(C(j))
    gvars_phitilde=[]
    for j in range(len(gvars)-1):
        gvars_phitilde.append(gvars[j]^p+p*gvars[j+1])
    # At the last stop, the delta term is beyond the g-precision
    # so it is not added.
    gvars_phitilde.append(C(0))
    
    # Gently coerce E into C. Occasionally, z is replaced by d otherwise, which
    # leads to nonsense.
    E_C=C(A(E))
      
    def weight(funct):
        # Returns a tuple (F-weight,N-weight), where F-weight is the exact (up to precision) F-weight
        # and N-weight is a *lower bound* on the Nygaard weight of f.
        # For example, this function will not correctly compute the Nygaard weight E(z), although
        # it will of d.
        nygaard_weight=p^(gprec-1)
        f_weight=zprec
        p_powers=[p^j for j in range(gprec)]
        
        def monomial_weight(m):
            # Returns the Nygaard weight of a monomial in the gj.
            monomial_nygaard_weight=sum([a*b for a,b in zip(p_powers,list(m.degrees()))])
            return (monomial_nygaard_weight,monomial_nygaard_weight)
        
        for m in funct.monomials():
            monomial_f_weight,monomial_nygaard_weight=monomial_weight(m)
            # c is a polynomial in d on power series in z.
            c=funct.monomial_coefficient(m)
            coefficient_nygaard_weight=p^(gprec-1)
            coefficient_f_weight=zprec
            deg=c.degree(d)
           
            for n in c.monomials():
                x=c.monomial_coefficient(n)
                v=x.valuation()
                coefficient_nygaard_weight=min(coefficient_nygaard_weight,n.degree(d))
                coefficient_f_weight=min(coefficient_f_weight,v)
            nygaard_weight=min(nygaard_weight,monomial_nygaard_weight+coefficient_nygaard_weight)
            f_weight=min(f_weight,monomial_f_weight+coefficient_f_weight)
        
        return (f_weight,nygaard_weight)
    
    def phitilde(funct): 
        g=0
        for m in funct.monomials():
            c=funct.monomial_coefficient(m)
            h=0
            # Iterate over powers n of d.
            for n in c.monomials():
                x=c.monomial_coefficient(n)
                h+=x.V(p)*n^p
            g+=C(m(gvars_phitilde)*h)
        return g
    
    def deltatilde(funct):
        return (C(phitilde(funct)-funct^p))*(1/p)

    def trim(funct):
        # Here, we use the fact that we have a fixed z-precision and the g_j secretly
        # have z-weight p^j. So, we can trim lots of information.
        trimmed=0
        for m in funct.monomials():
            w=weight(m)[0]
            new_coefficient=0
            c=funct.monomial_coefficient(m)
            for n in c.monomials():
                to_trim=c.monomial_coefficient(n).add_bigoh(max(zprec-w,0))
                new_coefficient+=to_trim*n
            trimmed+=new_coefficient*m
        return C(trimmed)
    
    rels=[]
    if gprec==1:
        rels.append(gvars[0]^p)
    else:
        if p==2:
            rels.append(trim(phitilde(E_C)*g1+(deltatilde(E_C)+E_C^2)*g0^2-z*E_C*d*g0))
        else:
            new_rel=phitilde(E_C)*g1+deltatilde(E_C)*g0^p
            for j in range(1,p):
                new_rel+=(-1)^(j+1)*(binomial(p,j)//p)*z^j*E_C^(p-j)*d^j*gvars[0]^(p-j)
            rels.append(trim(new_rel))
        for s in range(1,gprec-1):
            rels.append(trim(deltatilde(rels[s-1])))
        rels.append(gvars[gprec-1]^p)
        
    # new_rels is a lookup table for how to rewrite g_j^p in terms of other factors.
    if gprec==1:
        new_rels=rels
    else:
        new_rels=[]
        for j in range(len(rels)):
            r=C(rels[j])
            u=r.monomial_coefficient(gvars[j]^p)
            u=A(u)
            new_rels.append(trim(gvars[j]^p-(1/u)*r))
    
    def coefficient_divide(funct,gunct,zprec,weight):
        # Takes a polynomial gunct in the d,gj and divides all the
        # coefficients by funct.
        out=0
        for m in gunct.monomials():
            for n in gunct.monomial_coefficient(m).monomials():
                (f_weight,nygaard_weight)=weight(n*m)
                out+=lazy_division(funct,(gunct.monomial_coefficient(m)).monomial_coefficient(n),zprec-f_weight)*n*m
        return out
    
    def coefficient_reduce(c):
        # Takes a coefficient (a polynomial in d with coefficients in power series in z)
        # and returns an element of C obtained by using d=E
        new_coefficient=0
        for n in c.monomials():
            new_coefficient+=c.monomial_coefficient(n)*E_C^(n.degree(d))
        return new_coefficient    

    def reduce_d(funct):
        reduced=0
        for m in funct.monomials():
            reduced+=coefficient_reduce(funct.monomial_coefficient(m))*m
        return trim(C(reduced))
    
    #now, reduce new_rels by rewriting all d's by E(z), since we're done with delta_tilde:
    for j in range(len(new_rels)):
        new_rels[j]=reduce_d(new_rels[j])   
    
    def reduce(funct):
        # Run through the monomials of g and rewrite them once using new_rels.
        # assumes they are already reduce_d'd.
        reduced=0
        for m in funct.monomials():
            new_monomial=1
            degs=m.degrees()
            for exponent in range(len(degs)):
                r=degs[exponent]//p
                new_monomial=new_monomial*new_rels[exponent]^r*gvars[exponent]^(degs[exponent]-p*r)
            reduced+=funct.monomial_coefficient(m)*new_monomial
        reduced=trim(C(reduced))                
        return reduced
    
    #now, we reduce the relations in new_rels, using themselves.
    #We don't do this using recursive_reduce below,
    #since every step we can overwrite the previous new_rels to speed up convergence
    for j in range(len(new_rels)):
        prev_value=new_rels[j]
        new_rels[j]=reduce(C(prev_value))
        while prev_value != new_rels[j]:
            prev_value=new_rels[j]
            new_rels[j]=reduce(C(prev_value))
    
    #recursively reduces an arbitrary expression using reduce_d and the reduced new_rels.
    def recursive_reduce(funct):
        funct=reduce_d(funct)
        prevFunct=funct
        funct=reduce(C(funct))
        while prevFunct != funct:
            prevFunct = funct
            funct=reduce(C(funct))
        return funct
    
    def initialize_gvars_phi():
        # phi(g_j) = phi(E(z))^{p^j} * (g_j^p+p g_{j+1})/E(z)^{p^{j+1}}
        gvars_phi=[]
        for j in range(len(gvars)):
            X=phitilde(E_C)^(p^j)*phitilde(gvars[j])
            X=recursive_reduce(X)
            gvars_phi.append(coefficient_divide(A(E)^(p^(j+1)),X,zprec,weight))
        return gvars_phi
    gvars_phi=initialize_gvars_phi()

    def phi(funct):
        g=0
        for m in funct.monomials():
            c=funct.monomial_coefficient(m)
            h=0
            # Iterate over powers n of d.
            for n in c.monomials():
                x=c.monomial_coefficient(n)
                deg=n.degree(d)
                h+=x.V(p)*C(E(z^p))^deg
            g+=C(m(gvars_phi)*h)
        return g
    
    def delta(funct):
        return (C(phi(funct)-funct^p))*(1/p)
    
    def g0_divide(funct):
        # In a polynomial of the form g0*F(g0) returns F(g0).
        divided_funct=0
        for m in funct.monomials():
            degs=list(m.degrees())
            if degs[0]==0:
                raise TypeError('Cannot divide by g0.')
            else:
                degs[0]=degs[0]-1
                divided_funct+=funct.monomial_coefficient(m)*C.monomial(*tuple(degs))
        return divided_funct
        
    return C,gvars,weight,phi,delta,phitilde,deltatilde,reduce,recursive_reduce,g0_divide,coefficient_divide
