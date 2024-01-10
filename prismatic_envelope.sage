from precision import *
from matrix_functions import *
set_verbose(-1)

"""
A module to compute p-adic syntomic cohomology of OK/pi^n.

WARNING: many of the functions in this module assume that the base coefficient ring W
and the base Breuil-Kisin ring A are in the namespace.
"""

####################
# Helper functions #
####################

def eisenstein_normalization(p,E):
    """Renormalizes Eisenstein polynomial.

    Takes an Eisenstein poylnomial E(z)=...+pu and returns (1/u)*E so that the constant term is
    +p.
    """
    return E/(E(0)/p)

def lazy_division(funct,gunct,Fprec):
    """Returns quotient + O(z^Fprec) such that quotient*funct=gunct mod z^Fprec."""
    a=funct.valuation()
    b=gunct.valuation()
    quotient=0
    new_gunct=gunct
    for j in range(b-a,Fprec):
        try:
            # Possible p-adic precision loss here.
            quotient+=z^j*W(new_gunct[j+a]/funct[a])
        except ValueError:
            raise ValueError('Function g={gstr} not divisible by f={fstr} up to precision={hstr}.'.format(gstr=gunct, fstr=funct, hstr=Fprec))
        new_gunct=gunct-funct*quotient
    return quotient


#######################
# Prismatic envelopes #
#######################

def prismatic_envelope_f(p,E,k,prec,Fprec,debug=False):
    """Returns B,grN,phi,delta, where B is a polynomial ring in f0,...,num_f-1,
    grN is a function which measures the Nygaard filtration of an element of
    B, phi is the graded lift of Frobenius, and delta is the graded
    delta-ring structure.

    F-weight z is 1, Nygaard of z is 0.
    F-weight of f_j is k*p^j, Nygaard weight of f_j is p^j.
    """
    if i==1:
        num_f=0
    else:
        num_f=floor(log(i-1,p))+1
    if num_f==1:
        B=PolynomialRing(A,'f0',1)
    else:
        B=PolynomialRing(A,'f',num_f)
    B.inject_variables(verbose=False)
    C=PolynomialRing(B,'d_tilde')
    variable_names=B.variable_names()
    C.inject_variables(verbose=False)
    # A list of the fi for ease of reference.
    fvars=[]
    for j in variable_names:
        fvars.append(B(j))
    # Viewing f_j has having Nygaard weight p^j, we have the graded Frobenius
    # phitilde(x)=x^p+p\delta_i(x) for a homogeneous element of Nygaard weight i on the direct
    # sum of the Nygaard pieces; see Lemma 4.3 of KZPN. [REF]
    #
    # We record the value of phitilde on f_j. As f_{j+1}=\delta_{p^j}(f_j) by definition, we have
    # phitilde(f_j)=f_j^p+p f_{j+1}.
    fvars_phitilde=[]
    for j in range(num_f-1):
        fvars_phitilde.append(fvars[j]^p+p*fvars[j+1])
    # At the last stop, the delta term is beyond the f-precision
    # so it is not added.
    # TODO: isn't fvars[num_f-1]^p *also* beyond the f-precision, so could we not set it to zero?
    if num_f>0:
        fvars_phitilde.append(fvars[num_f-1]^p)

    if debug:
        print("fvars_phitilde:")
        print(fvars_phitilde)
        print('\n')
    
    def weight(funct):
        """The weight of a function.

        Returns a tuple (F-weight,N-weight), where F-weight is the F-weight up to Fprec
        and N-weight is a *lower bound* on the Nygaard weight of f.
        For example, this function will not correctly compute the Nygaard weight of large powers of p.
        """
        F_weight=Fprec
        # The highest possible Nygaard weight is that of z^{k-1}f_{num_f-1}, which has Nygaard weight p^{num_f-1}.
        if num_f==0:
            nygaard_weight=0
        else:
            nygaard_weight=p^(num_f-1)
        p_powers=[p^j for j in range(num_f)]
        
        def monomial_weight(m):
            """The weight of a monomial.

            Returns a tuple (F-weight,N-weight), where F-weight is the exact F-weight
            and N-weight is the exact Nygaard weight of a monomial in the f_i.
            """
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
    
    def coefficient_divide_f(funct,gunct,Fprec,weight):
        """Takes a polynomial gunct in the fj and divides all the coefficients by funct."""
        out=0
        for m in gunct.monomials():
            (f_weight,nygaard_weight)=weight(m)
            new_coefficient=lazy_division(funct,gunct.monomial_coefficient(m),Fprec-f_weight)
            out+=new_coefficient*m
        return out
    
    def trim(funct):
        """Trims function.

        Here, we use the fact that we have a fixed z-precision and the f_j secretly
        have z-weight kp^j. So, we can trim lots of information.
        """
        g=0
        for m in funct.monomials():
            w=weight(m)[0]
            new_coefficient=funct.monomial_coefficient(m).add_bigoh(max(Fprec-w,0))
            g+=new_coefficient*m
        return B(g)
    
    rels=[]
    if num_f==0:
        pass
    elif num_f==1:
        # CHANGED: rels.append(fvars[0]^p)
        # TODO: could just be set to 0, right?
        rels.append(fvars_phitilde[0])
    else:
        # The following is the first relation in Lemma 4.9 [REF] under the assumption that r has
        # rank 1 so that \delta(r)=0.
        rels.append(trim(phitilde(B(E))*f1 + deltatilde(B(E))*f0^p))
        for j in range(1,num_f-1):
            rels.append(trim(deltatilde(rels[j-1])))
        # CHANGED: rels.append(fvars[num_f-1]^p)
        # TODO: could just be set to 0, right?
        rels.append(fvars_phitilde[num_f-1])

    if debug:
        print("rels:")
        print(rels)
        print('\n')
        
    # new_rels is a lookup table for how to rewrite f_j^p in terms of other factors.
    new_rels=[]
    for j in range(len(rels)):
        r=B(rels[j])
        u=r.monomial_coefficient(fvars[j]^p)
        new_rels.append(trim(fvars[j]^p)-(1/u)*r)
 
    def reduce(funct):
        """Run through the monomials of f and rewrite them once using new_rels."""
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
        fvars_phi_divided.append(coefficient_divide_f(A(E)^(p^(j+1)),phitilde_reduced,Fprec,weight))
    if num_f>0:
        fvars_phi_divided.append(B(0))
    if debug:
        print("fvars divided frobeniuses")
        print(fvars_phi_divided)
        print('\n')
    
    def reduce_zk(funct):
        reduced=0
        def coefficient_reduce(c):
            # Takes a coefficient (a power series in z)
            # and returns an element of B obtained by using z^k=f_0
            if num_f>0:
                out=B(c[0])
                for j in range(1,min(Fprec,c.precision_absolute())):
                    r=j//k
                    out+=fvars[0]^r*z^(j-k*r)*c[j]
                return B(out)
            else:
                return c
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
        result=B(0)
        for m in funct.monomials():
            c=funct.monomial_coefficient(m)
            result+=B(c.V(p)*m(fvars_phi_divided))
        return result

    def expand_d_tilde(target_nygaard_weight,funct):
        new_funct=C(0)
        for m in funct.monomials():
            m_coefficient=funct.monomial_coefficient(m)
            for t in m_coefficient.monomials():
                a=min(weight(B(t))[1]+m.degree()-target_nygaard_weight,m.degree())
                new_funct+=m_coefficient.monomial_coefficient(t)*E^a*t*d_tilde^(m.degree()-a)
        return new_funct

    def coefficient_reduce(funct):
        new_funct=C(0)
        for m in funct.monomials():
            new_funct+=C(recursive_reduce(funct.monomial_coefficient(m)))*m
        return new_funct

    def recreduce_nygaard(target_nygaard_weight,funct):
        # Reduces a polynomial in d_tilde^A z^j\prod f_j^{a_j}
        # subject to being in N^{>=target_nygaard_weight}.
        for m in funct.monomials():
            if weight(funct.monomial_coefficient(m))[1]+m.degree()<j:
                raise ValueError("Input funct has Nygaard weight smaller than target_nygaard_weight.")
        prev_value = funct
        funct = coefficient_reduce(funct)
        funct = expand_d_tilde(target_nygaard_weight,funct)
        while prev_value!=funct:
            prev_value = funct
            funct = coefficient_reduce(funct)
            funct = expand_d_tilde(target_nygaard_weight,funct)
        return funct

    return B,C,fvars,weight,phitilde,phi_divided,deltatilde,reduce,recursive_reduce,recreduce_nygaard


    
def prismatic_envelope_g(p,E,k,prec,Fprec,debug=False):
    """For working with O_K.

    Returns C,grN,phi,delta, where C is a polynomial ring in g0,...,num_g-1,
    grN is a function which measures the Nygaard filtration of an element of
    C, phi is the graded lift of Frobenius, and delta is the graded
    delta-ring structure.
    
    z has F-weight 1 and Nygaard-weight 0.
    d has F-weight 0 and Nygaard-weight 1.
    gj has F-weight p^j and Nygaard-weight p^j.
    """
    A0.<d>=PolynomialRing(A)
    A0.inject_variables(verbose=False)
    num_g=floor(log(n*i-1,p))+1
    C=PolynomialRing(A0,'g',num_g)
    C.inject_variables(verbose=False)
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
        """Returns a tuple of weights.

        Returns a tuple (F-weight,N-weight), where F-weight is the exact (up to precision) F-weight
        and N-weight is a *lower bound* on the Nygaard weight of f.
        For example, this function will not correctly compute the Nygaard weight E(z), although
        it will of d.
        """
        nygaard_weight=p^(num_g-1)
        f_weight=Fprec
        p_powers=[p^j for j in range(num_g)]
        
        def monomial_weight(m):
            """Returns the Nygaard weight of a monomial in the gj."""
            monomial_nygaard_weight=sum([a*b for a,b in zip(p_powers,list(m.degrees()))])
            return (monomial_nygaard_weight,monomial_nygaard_weight)
        
        for m in funct.monomials():
            monomial_f_weight,monomial_nygaard_weight=monomial_weight(m)
            # c is a polynomial in d on power series in z.
            c=funct.monomial_coefficient(m)
            coefficient_nygaard_weight=p^(num_g-1)
            coefficient_f_weight=Fprec
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
        """Trim function.

        Here, we use the fact that we have a fixed z-precision and the g_j secretly
        have z-weight p^j. So, we can trim lots of information.
        """
        trimmed=0
        for m in funct.monomials():
            w=weight(m)[0]
            new_coefficient=0
            c=funct.monomial_coefficient(m)
            for n in c.monomials():
                to_trim=c.monomial_coefficient(n).add_bigoh(max(Fprec-w,0))
                new_coefficient+=to_trim*n
            trimmed+=new_coefficient*m
        return C(trimmed)
    
    rels=[]
    if num_g==1:
        rels.append(gvars[0]^p)
    else:
        # This time, the relations from Lemma 4.9 of KZPN [REF] are more complicated.
        if p==2:
            rels.append(trim(phitilde(E_C)*g1+(deltatilde(E_C)+E_C^2)*g0^2-z*E_C*d*g0))
        else:
            new_rel=phitilde(E_C)*g1+deltatilde(E_C)*g0^p
            for j in range(1,p):
                new_rel+=(-1)^(j+1)*(binomial(p,j)//p)*z^j*E_C^(p-j)*d^j*gvars[0]^(p-j)
            rels.append(trim(new_rel))
        for s in range(1,num_g-1):
            rels.append(trim(deltatilde(rels[s-1])))
        rels.append(gvars[num_g-1]^p)

    if debug:
        print("rels in prismatic_envelope_g:")
        print(rels)
        print('\n')
        
    # new_rels is a lookup table for how to rewrite g_j^p in terms of other factors.
    if num_g==1:
        new_rels=rels
    else:
        new_rels=[]
        for j in range(len(rels)):
            r=C(rels[j])
            u=r.monomial_coefficient(gvars[j]^p)
            u=A(u)
            new_rels.append(trim(gvars[j]^p-(1/u)*r))

    if debug:
        print("new_rels before reduction in prismatic_envelope_g:")
        print(new_rels)
        print('\n')
    
    def coefficient_divide(funct,gunct,Fprec,weight):
        """Takes a polynomial gunct in the d,gj and divides all the coefficients by funct."""
        out=0
        for m in gunct.monomials():
            for n in gunct.monomial_coefficient(m).monomials():
                (f_weight,nygaard_weight)=weight(n*m)
                out+=lazy_division(funct,(gunct.monomial_coefficient(m)).monomial_coefficient(n),Fprec-f_weight)*n*m
        return out
    
    def coefficient_reduce(c):
        """Reduces coefficient

        Takes a coefficient (a polynomial in d with coefficients in power series in z)
        and returns an element of C obtained by using d=E
        """
        new_coefficient=0
        for n in c.monomials():
            new_coefficient+=c.monomial_coefficient(n)*E_C^(n.degree(d))
        return new_coefficient    

    def reduce_d(funct):
        reduced=0
        for m in funct.monomials():
            reduced+=coefficient_reduce(funct.monomial_coefficient(m))*m
        return trim(C(reduced))
    
    # Now, reduce new_rels by rewriting all d's by E(z), since we're done with delta_tilde:
    for j in range(len(new_rels)):
        new_rels[j]=reduce_d(new_rels[j])   
    
    def reduce(funct):
        """Run through the monomials of g and rewrite them once using new_rels.

        assumes they are already reduce_d'd.
        """
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
    
    # Now, we reduce the relations in new_rels, using themselves.
    # We don't do this using recursive_reduce below,
    # since in every step we can overwrite the previous new_rels to speed up convergence
    for j in range(len(new_rels)):
        prev_value=new_rels[j]
        new_rels[j]=reduce(C(prev_value))
        while prev_value != new_rels[j]:
            prev_value=new_rels[j]
            new_rels[j]=reduce(C(prev_value))

    if debug:
        print("new_rels in prismatic_envelope_g:")
        print(new_rels)
        print('\n')
    
    def recursive_reduce(funct):
        """Recursively reduces an arbitrary expression using reduce_d and the reduced new_rels."""
        funct=reduce_d(funct)
        prevFunct=funct
        funct=reduce(C(funct))
        while prevFunct != funct:
            prevFunct = funct
            funct=reduce(C(funct))
        return funct
    
    def initialize_gvars_phi():
        """phi(g_j) = phi(E(z))^{p^j} * (g_j^p+p g_{j+1})/E(z)^{p^{j+1}}"""
        gvars_phi=[]
        for j in range(len(gvars)):
            X=phitilde(E_C)^(p^j)*phitilde(gvars[j])
            X=recursive_reduce(X)
            gvars_phi.append(coefficient_divide(A(E)^(p^(j+1)),X,Fprec,weight))
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
        """In a polynomial of the form g0*F(g0) returns F(g0)."""
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


###############################################
# Build nabla on prismatic cohomology for OK. #
###############################################

def nablaP_matrix_OK(p,i,k,E,prec,Fprec,debug=False):
    """Returns nablaP_OK.

    Note that this does not technically depend on k, so one can run this once
    for a large k and slice down to compute for several k at once.
    """
    C,gvarsC,weightC,phiC,deltaC,phiCtilde,deltaCtilde,reduceC,recreduceC,g0_divideC,coefficient_divideC=prismatic_envelope_g(p,E,k,prec,Fprec,debug=debug)
    
    def initialize_u():
        # Takes almost no time.
        if len(gvarsC)==1:
            divided_g0=gvarsC[0]^p
        else:
            divided_g0=(gvarsC[0]^p+p*gvarsC[1])
        divided_g0=recreduceC(divided_g0)
        divided_g0=coefficient_divideC(A(E)^p,divided_g0,Fprec,weightC)

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
    if debug:
        print("bk factor is")
        print(bk_factor)
    
    # nablaP builder
    # First, nablaP_OK and then a nablaP_OKpik. We get the second using a saturation
    # method using the maps from prismaOK to prism OKpiK.

    # nablaP_OK builder
    # Element of BK twist i is x*s^i, x in Prism and s is the BK orientation.
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
    # Figure num_g.
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
            
    nablaP_OK=initialize_nablaP_OK()
    if debug:
        print("nablaP_OK is")
        print(nablaP_OK)
    return nablaP_OK


#########################################################
# The main function which collates everything together. #
#########################################################

def syntomic_matrices(p,i,k,E,prec,Fprec,nablaP_OK=False,debug=False):
    """Returns syn0,syn1,nablaN,nablaP.

    The slowest part of the computation is the computation of the preliminary
    matrix nablaP_OK, so this can alternatively be passed as an argument
    in the case it is precomputed.
    """
    B,C,fvars,weightB,phiBtilde,phi_dividedB,deltaBtilde,reduceB,recreduceB,recreduceN=prismatic_envelope_f(p,E,k,prec,Fprec,debug=debug)
    if i==1:
        num_f=0
    else:
        num_f=floor(log(i-1,p))+1

    # Basis z^c\prod_{j=0}^{num_f-1}f_j^{a_j}, 0<=a_j<=p-1, 0<=c<=k-1 with one exception: no 1.

    # can0 builder
    can0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(E(z)^(i-(n//k))*z^c*fprod)
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            can0[m-1,n-1]=coefficient_to_process[a]

    if debug:
        print("can0:")
        print(can0)
        print('\n')

    # phi0 builder
    phi0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(z^(p*c)*phi_dividedB(fprod))
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            phi0[m-1,n-1]=coefficient_to_process[a]

    # can1 builder
    # Same as above but with basis starting with 1 and going up to z^{k-2}*num_f or z^{k-1}f_{prec-1}
    # Does not make great sense for k=1.
    can1=Matrix(W,k*i-1,k*i-1)
    for n in range(0,k*i-1):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(E(z)^(i-(n//k)-1)*z^c*fprod)
        for m in range(0,k*i-1):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            can1[m,n]=coefficient_to_process[a]

    if debug:
        print('can1:')
        print(can1)
        print('\n')

    # nablaP builder
    # First, nablaP_OK and then a nablaP_OKpik. We get the second using a saturation
    # method using the maps from prismaOK to prism OKpiK.
    if nablaP_OK==False:
        nablaP_OK=nablaP_matrix_OK(p,i,k,E,prec,Fprec,debug=debug)
   
    OKtoOKmodpi0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        column_to_process=recreduceB(B(z^n))
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            OKtoOKmodpi0[m-1,n-1]=coefficient_to_process[a]

    OKtoOKmodpi1=Matrix(W,k*i-1,k*i-1)
    for n in range(0,k*i-1):
        column_to_process=recreduceB(B(z^n))
        for m in range(0,k*i-1):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            OKtoOKmodpi1[m,n]=coefficient_to_process[a]
    
    # Given nablaP for OK, compute it for OK/pi^k.
    OKtoOKmodpi0,OKtoOKmodpi1,nablaP_OK,nablaP=square_complete_bottom(OKtoOKmodpi0,OKtoOKmodpi1,nablaP_OK)

    # The relative syntomic matrix.
    syn0=can0-phi0

    # Compute nablaN for OK/pi^k.
    can0,can1,nablaN,nablaP=square_complete_top(can0,can1,nablaP)

    # Compute syn1 by completing the square.
    syn0,syn1,nablaN,nablaP=square_complete_right(syn0,nablaN,nablaP)

    return syn0,syn1,nablaN,nablaP


def v1_matrices(p,i,k,E,prec,Fprec,debug=False):
    # WARNING: padic precision calcluations for v1 have not been checked.
    if i-p+1<=0:
        raise NotImplementedError
    B,C,fvars,weightB,phiBtilde,phi_dividedB,deltaBtilde,reduceB,recreduceB,recreduceN=prismatic_envelope_f(p,E,k,prec,Fprec,debug=debug)
    if i==1:
        num_f=0
    else:
        num_f=floor(log(i-1,p))+1

    # v1_on_P0 builder
    v1P0=Matrix(W,k*i-1,k*(i-p+1)-1)
    for n in range(1,(i-p+1)*k):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        # We use c+1 to denote that we've multiplied by 'z^p'.
        column_to_process=recreduceB(z^(c+p)*fprod)
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            v1P0[m-1,n-1]=coefficient_to_process[a]

    # v1_on_N0 builder
    v1N0=Matrix(W,k*i-1,k*(i-p+1)-1)
    for n in range(1,(i-p+1)*k):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        # We use c+1 to denote that we've multiplied by 'z^p'.
        reduced_form=recreduceN(i,C(z^c*fprod)*d_tilde^(i-p+1-(n//k)+p))
        column_to_process=B(0)
        for cffcnt in reduced_form.coefficients():
            column_to_process+=B(cffcnt)
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            v1N0[m-1,n-1]=coefficient_to_process[a]

    # v1_on_P1 builder
    v1P1=Matrix(W,k*i-1,k*(i-p+1)-1)
    for n in range(0,(i-p+1)*k-1):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        # We use c+1 to denote that we've multiplied by 'z^p'.
        column_to_process=recreduceB(z^(c+p)*fprod)
        for m in range(0,i*k-1):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            v1P1[m,n]=coefficient_to_process[a]

    # v1_on_N1 builder
    v1N1=Matrix(W,k*i-1,k*(i-p+1)-1)
    for n in range(0,(i-p+1)*k-1):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        print("fprod and c are {} and {}".format(fprod,c))
        # We use c+1 to denote that we've multiplied by 'z^p'.
        print("input is {}".format(C(z^c*fprod)*d_tilde^(i-p-(n//k)+p)))
        reduced_form=recreduceN(i-1,C(z^c*fprod)*d_tilde^(i-p-(n//k)+p))
        print("reduced_form is {}".format(reduced_form))
        column_to_process=B(0)
        for cffcnt in reduced_form.coefficients():
            column_to_process+=B(cffcnt)
        print(column_to_process)
        print('\n')
        for m in range(0,i*k-1):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            print("gprod and a are {} and {}".format(gprod,a))
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            print("coefficient to process is {}".format(coefficient_to_process))
            v1N1[m,n]=coefficient_to_process[a]

        return v1N0,v1P0,v1N1,v1P1

###########
# Old v1. #
###########

def syntomic_matrices_v1(p,i,k,E,prec,Fprec,debug=False):
    # WARNING: this is not correct.
    if i-p+1<=0:
        raise NotImplementedError
    B,C,fvars,weightB,phiBtilde,phi_dividedB,deltaBtilde,reduceB,recreduceB,recreduceN=prismatic_envelope_f(p,E,k,prec,Fprec,debug=debug)
    if i==1:
        num_f=0
    else:
        num_f=floor(log(i-1,p))+1

    # First compute in weight i.

    # Basis z^c\prod_{j=0}^{num_f-1}f_j^{a_j}, 0<=a_j<=p-1, 0<=c<=k-1 with one exception: no 1.

    # can0 builder
    b_can0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(E(z)^(i-(n//k))*z^c*fprod)
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            b_can0[m-1,n-1]=coefficient_to_process[a]

    if debug:
        print("b_can0:")
        print(b_can0)
        print('\n')

    # phi0 builder
    b_phi0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(z^(p*c)*phi_dividedB(fprod))
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            b_phi0[m-1,n-1]=coefficient_to_process[a]

    # can1 builder
    # Same as above but with basis starting with 1 and going up to z^{k-2}*num_f or z^{k-1}f_{prec-1}
    # Does not make great sense for k=1.
    b_can1=Matrix(W,k*i-1,k*i-1)
    for n in range(0,k*i-1):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(E(z)^(i-(n//k)-1)*z^c*fprod)
        for m in range(0,k*i-1):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            b_can1[m,n]=coefficient_to_process[a]

    if debug:
        print('b_can1:')
        print(b_can1)
        print('\n')

    b_nablaP_OK=nablaP_matrix_OK(p,i,k,E,prec,Fprec,debug=debug)
   
    b_OKtoOKmodpi0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        column_to_process=recreduceB(B(z^n))
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            b_OKtoOKmodpi0[m-1,n-1]=coefficient_to_process[a]

    b_OKtoOKmodpi1=Matrix(W,k*i-1,k*i-1)
    for n in range(0,k*i-1):
        column_to_process=recreduceB(B(z^n))
        for m in range(0,k*i-1):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            b_OKtoOKmodpi1[m,n]=coefficient_to_process[a]
    
    # Given nablaP for OK, compute it for OK/pi^k.
    b_OKtoOKmodpi0,b_OKtoOKmodpi1,b_nablaP_OK,b_nablaP=square_complete_bottom(b_OKtoOKmodpi0,b_OKtoOKmodpi1,b_nablaP_OK)

    # The relative syntomic matrix.
    b_syn0=b_can0-b_phi0

    # Compute nablaN for OK/pi^k.
    b_can0,b_can1,b_nablaN,b_nablaP=square_complete_top(b_can0,b_can1,b_nablaP)

    # Compute syn1 by completing the square.
    b_syn0,b_syn1,b_nablaN,b_nablaP=square_complete_right(b_syn0,b_nablaN,b_nablaP)

    # Then, compute in weight i-p+1.
    i=i-p+1

    # Basis z^c\prod_{j=0}^{num_f-1}f_j^{a_j}, 0<=a_j<=p-1, 0<=c<=k-1 with one exception: no 1.

    # can0 builder
    a_can0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(E(z)^(i-(n//k))*z^c*fprod)
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            a_can0[m-1,n-1]=coefficient_to_process[a]

    if debug:
        print("a_can0:")
        print(a_can0)
        print('\n')

    # phi0 builder
    a_phi0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(z^(p*c)*phi_dividedB(fprod))
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            a_phi0[m-1,n-1]=coefficient_to_process[a]

    # can1 builder
    # Same as above but with basis starting with 1 and going up to z^{k-2}*num_f or z^{k-1}f_{prec-1}
    # Does not make great sense for k=1.
    a_can1=Matrix(W,k*i-1,k*i-1)
    for n in range(0,k*i-1):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        column_to_process=recreduceB(E(z)^(i-(n//k)-1)*z^c*fprod)
        for m in range(0,k*i-1):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            a_can1[m,n]=coefficient_to_process[a]

    if debug:
        print('a_can1:')
        print(a_can1)
        print('\n')

    a_nablaP_OK=nablaP_matrix_OK(p,i,k,E,prec,Fprec,debug=debug)
   
    a_OKtoOKmodpi0=Matrix(W,k*i-1,k*i-1)
    for n in range(1,i*k):
        column_to_process=recreduceB(B(z^n))
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            a_OKtoOKmodpi0[m-1,n-1]=coefficient_to_process[a]

    a_OKtoOKmodpi1=Matrix(W,k*i-1,k*i-1)
    for n in range(0,k*i-1):
        column_to_process=recreduceB(B(z^n))
        for m in range(0,k*i-1):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            a_OKtoOKmodpi1[m,n]=coefficient_to_process[a]
    
    # Given nablaP for OK, compute it for OK/pi^k.
    a_OKtoOKmodpi0,a_OKtoOKmodpi1,a_nablaP_OK,a_nablaP=square_complete_bottom(a_OKtoOKmodpi0,a_OKtoOKmodpi1,a_nablaP_OK)

    # The relative syntomic matrix.
    a_syn0=a_can0-a_phi0

    # Compute nablaN for OK/pi^k.
    a_can0,a_can1,a_nablaN,a_nablaP=square_complete_top(a_can0,a_can1,a_nablaP)

    # Compute syn1 by completing the square.
    a_syn0,a_syn1,a_nablaN,a_nablaP=square_complete_right(a_syn0,a_nablaN,a_nablaP)

    # Now, reset i.
    i=i+p-1

    # v1_on_P0 builder
    v1P0=Matrix(W,k*i-1,k*(i-p+1)-1)
    for n in range(1,(i-p+1)*k):
        n=ZZ(n)
        c=n.mod(k)
        d=W(n//k)
        fprod=B(1)
        for j in range(num_f):
            fprod=fprod*fvars[j]^(d[j])
        # We use c+1 to denote that we've multiplied by 'z^p'.
        column_to_process=recreduceB(z^(c+p)*fprod)
        for m in range(1,i*k):
            m=ZZ(m)
            a=m.mod(k)
            b=W(m//k)
            gprod=B(1)
            for j in range(num_f):
                gprod=gprod*fvars[j]^(b[j])
            coefficient_to_process=column_to_process.monomial_coefficient(gprod)
            v1P0[m-1,n-1]=coefficient_to_process[a]

    print(a_can0)
    print(b_can0)
    print(v1P0)

    a_can0,b_can0,v1N0,v1P0=square_complete_top(a_can0,b_can0,v1P0)
    a_nablaP,b_nablaP,v1P1,v1P0=square_complete_bottom(a_nablaP,b_nablaP,v1P0)
    a_nablaN,b_nablaN,v1N1,v1N0=square_complete_bottom(a_nablaN,b_nablaN,v1N0)

    return a_syn0,a_syn1,a_nablaN,a_nablaP,b_syn0,b_syn1,b_nablaN,b_nablaP,v1N0,v1P0,v1N1,v1P1


           
###########################################################################
# Compute the cohomology of the syntomic complex via elementary divisors. #
###########################################################################

def syntomic_cohomology(syn0,syn1,nablaN,nablaP):
    """Returns the cohomology of the total complex together with the precision of the input matrices.

    Out of laziness, the current implementation assumes
    the rational homology of the total complex is 0.
    This means that H^2 can be computed using elementary divisors of d^1
    and then H^0 and H^1 can be computed using the elementary divisors
    of the saturation of d^0.
    Returns a dictionary 'hi':(rank,[list of torsion orders]) and the precision of the matrices.
    """
    d0=block_matrix([[syn0],[nablaN]])
    d1=block_matrix([[nablaP,-syn1]])
    # This is really only here because padic matrices
    # don't have ranks :(
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
    return dict,min(matrix_precision(syn0),
            matrix_precision(syn1),
            matrix_precision(nablaN),
            matrix_precision(nablaP))
