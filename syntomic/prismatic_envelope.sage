class PrismaticEnvelopeF():
    def __init__(self,p,n,eisenstein, prec_F, prec_p, debug=False):
        self.debug = debug
        self.p = p
        self.n = n
        self.eisenstein = eisenstein
        self.prec_F = prec_F
        self.prec_p = prec_p
        self._init_ring()
        self._init_lambda()
        self._init_relations()
        self._reduce_relations()
        self._init_nygaard_ring()
        self._init_phi_divided()
        self._compute_basis()
    def _init_ring(self):
        if(self.debug):
            print('initializing ring')
        self.basering = self.eisenstein.parent()
        self.z = self.basering.gen()
        if (self.prec_F-1)//self.n == 0:
            num_f = 0
        else:
            num_f=floor(log((self.prec_F-1)//self.n,self.p))+1
        if num_f==1:
            self.ring=PolynomialRing(self.basering,'f0',1)
        else:
            self.ring=PolynomialRing(self.basering,'f',num_f)
        self.f = self.ring.gens()
    def weight_F(self,m):
        result = 0
        for i in range(len(self.f)):
            result += self.n * self.p^i * m.degrees()[i]
        return result
    def weight_N(self,m):
        result = 0
        for i in range(len(self.f)):
            result += self.p^i * m.degrees()[i]
        return result
    def phi_base(self,val):
        #this should ideally be in the "delta-ring base class" at some point
        return val.V(self.p)
    def _divide_coefficients_base(self,val):
        modp = val.map_coefficients(lambda x: x % self.p)
        assert(modp==0)
        result = val.map_coefficients(lambda x: x // self.p)
        return result
    def delta_base(self,val):
        assert(val.parent() == self.basering)
        result = self._divide_coefficients_base(self.phi_base(val)-val^self.p)
        return result
    def trim(self,val):
        result = self.ring(0)
        for m in val.monomials():
            w=self.weight_F(m)
            if w>=self.prec_F:
                continue
            new_c = val.monomial_coefficient(m).add_bigoh(self.prec_F-w)
            result += self.ring.from_base_ring(new_c)*m
        return result
    def _init_relations(self):
        if(self.debug):
            print('initializing relations')
        if len(self.f) == 0:
            self._relations = []
            return
        if len(self.f) == 1:
            self._relations = [self.ring(0)]
            return
        self._relations=[]
        d = self.ring.from_base_ring(self.eisenstein)
        dpow = d
        for j in range(len(self.f)-1):
            dpow = dpow^self.p
            val = (-self.p + dpow * self._lambda[j])*self.f[j+1]
            self._relations.append(val)
        self._relations.append(self.ring(0))
    def _reduce_coefficient(self,val):
        result = self.ring(0)
        poly = val.polynomial()
        for m in poly.monomials():
            c = poly.monomial_coefficient(m)
            deg = m.degree()
            r = deg//self.n
            if(r == 0):
                result += self.ring.from_base_ring(c*m)
            elif(len(self.f)>0):
                result += self.ring.from_base_ring(c*self.z^(deg-self.n*r))*self.f[0]^r
            #else, z^n=0 and we add nothing
        return result
    def _reduce_single_step(self,val):
        result=self.ring(0)
        for m in val.monomials():
            new_summand=self._reduce_coefficient(val.monomial_coefficient(m))
            degs=m.degrees()
            for exponent in range(len(degs)):
                r=degs[exponent]//self.p
                new_summand*=self._relations[exponent]^r*self.f[exponent]^(degs[exponent]-self.p*r)
            result += new_summand
        return self.trim(result)
    # old version. Much slower, since a lot of terms in f0 experience Nygaard explosion using unoptimized higher relations in the process.
    # new version below cycles the relations to alleviate that.
    # the even older version which we used to produce the tables actually does not use the relation z^n=f_0 at all while reducing relations.
    # this makes creating the prismatic envelope fast, but leaves lots of Nygaard explosions happening in the computation of the syntomic matrices (I think)
    #def _reduce_relations(self):
    #    for j in range(len(self._relations)):
    #        prev_val = self._relations[j]
    #        self._relations[j] = self._reduce_single_step(prev_val)
    #        while(prev_val != self._relations[j]):
    #            prev_val = self._relations[j]
    #            self._relations[j] = self._reduce_single_step(prev_val)


    def _reduce_relations_other(self):
        done = False
        while not done:
            done = True
            for j in range(len(self._relations)):
                prev_val = self._relations[j]
                coeff = self._relations[j].monomial_coefficient(self.f[j]^self.p)
                if coeff != 0:
                    self._relations[j] = (1/(1-coeff)) * (self._relations[j] - coeff * self.f[j]^self.p)
                else:
                    self._relations[j] = self._reduce_single_step(prev_val)
                if(prev_val != self._relations[j]):
                    if(self.debug):
                        print("f{} not done".format(j))
                    done=False
                else:
                    if(self.debug):
                        print("f{} done".format(j))
    def _reduce_relations(self):
        if(self.debug):
            print('reducing relations')
        done = False
        while not done:
            done = True
            for j in range(len(self._relations)):
                prev_val = self._relations[j]
                self._relations[j] = self._reduce_single_step(prev_val)
                if(prev_val != self._relations[j]):
                    if(self.debug):
                        print("f{} not done".format(j))
                    done=False
                else:
                    if(self.debug):
                        print("f{} done".format(j))
    def reduce(self,val):
        prev_val=val
        result=self._reduce_single_step(prev_val)
        while(prev_val != result):
            prev_val = result
            result = self._reduce_single_step(prev_val)
        return result
    def _init_nygaard_ring(self):
        if(self.debug):
            print('initializing nygaard ring')
        self.nygaard_ring = PolynomialRing(self.ring, 'dtilde')
        self.dtilde = self.nygaard_ring.gen()
    def _expand_dtilde(self,target_nygaard, val):
        result = self.nygaard_ring(0)
        for m in val.monomials():
            coefficient=val.monomial_coefficient(m)
            for t in coefficient.monomials():
                a=min(self.weight_N(t)+m.degree()-target_nygaard,m.degree())
                if(a<0):
                    raise ValueError('Input value has Nygaard weight smaller than target_nygaard')
                d = self.nygaard_ring.from_base_ring(self.ring.from_base_ring(self.eisenstein))
                result+=coefficient.monomial_coefficient(t)*d^a*t*self.dtilde^(m.degree()-a)
        return result

    def _nygaard_reduce_f_monomials(self,val):
        val = self.nygaard_ring(val)
        result=self.nygaard_ring(0)
        for m in val.monomials():
            result+=self.nygaard_ring.from_base_ring(self.reduce(val.monomial_coefficient(m)))*m
        return result
        
    def nygaard_reduce(self,target_nygaard, val):
        assert val.parent()==self.nygaard_ring, "in nygaard_reduce: type {}, should be {}".format(val.parent(), self.nygaard_ring)
        prev_value = val
        result = self._nygaard_reduce_f_monomials(prev_value)
        result = self._expand_dtilde(target_nygaard,result)
        while prev_value != result:
            prev_value = result
            result = self._nygaard_reduce_f_monomials(prev_value)
            result = self._expand_dtilde(target_nygaard,result)
        return result
    def can(self,val):
        assert val.parent()==self.nygaard_ring, "in can: type {}, should be {}".format(val.parent(), self.nygaard_ring)
        reduced = self.nygaard_reduce(0, val)
        result = reduced.monomial_coefficient(self.nygaard_ring(1))
        return result
    def _init_lambda(self):
        if(self.debug):
            print('initializing lambda')
        d = self.eisenstein
        self._lambda=[-1/self.delta_base(d)]
        d_pow = d
        for j in range(1,len(self.f)):
            d_pow = d_pow^self.p
            next_lambda = self._lambda[j-1]^self.p
            next_lambda *= 1/(1 - self.delta_base(d_pow*self._lambda[j-1]))
            self._lambda.append(next_lambda)
    def _init_phi_divided(self):
        if(self.debug):
            print('initializing phi_divided')

        self._phi_divided_f=[]
        for j in range(len(self.f)-1):
            val = self._lambda[j] * self.f[j+1]
            self._phi_divided_f.append(val)
        self._phi_divided_f.append(self.ring(0))
    def phi_divided(self,i,val):
        val = self.nygaard_reduce(i,val)
        result=self.ring(0)
        for dpow in val.monomials():
            weight_dtilde = dpow.degree()
            c = val.monomial_coefficient(dpow)
            for mon in c.monomials():
                weight_f = self.weight_N(mon)
                c2 = c.monomial_coefficient(mon)
                if(len(self.f)>0):
                    result += self.ring.from_base_ring(self.phi_base(c2) * self.phi_base(self.eisenstein)^(weight_dtilde + weight_f - i)) * mon(self._phi_divided_f)
                else:    
                    result += self.ring.from_base_ring(self.phi_base(c2) * self.phi_base(self.eisenstein)^(weight_dtilde - i))
        return result
    def _compute_basis(self):
        if(self.debug):
            print('computing basis')

        self._basis = []
        for j in range(self.prec_F):
            e = j
            monomial = self.ring(1)
            monomial *= self.z^(e % self.n)
            e = e//self.n
            i = 0
            while(e > 0):
                monomial *= self.f[i]^(e % self.p)
                e = e//self.p
                i += 1
            self._basis.append(monomial)
    def basis(self):
        return self._basis
    def nygaard_basis(self,i):
        result = []
        for j in range(self.prec_F):
            monomial = self.nygaard_ring(self._basis[j])
            a = max(i - j // self.n, 0) 
            monomial *= self.dtilde^a
            result.append(monomial)
        return result
    def element_to_vector(self,val):
        assert(val.parent() == self.ring)
        val = self.reduce(val)
        result = self.prec_F * [0]
        for m in val.monomials():
            deg_f = self.weight_F(m)
            c = val.monomial_coefficient(m).polynomial()
            for zpow in c.monomials():
                deg_z = zpow.degree()
                coeff = c.monomial_coefficient(zpow)
                result[deg_f+deg_z] += coeff
        return vector(self.basering.base_ring(), result)

    def nygaard_element_to_vector(self,i,val):
        assert(val.parent() == self.nygaard_ring)
        val = self.nygaard_reduce(i,val)
        result = self.prec_F * [0]
        for dpow in val.monomials():
            c1 = val.monomial_coefficient(dpow)
            for fmon in c1.monomials():
                deg_f = self.weight_F(fmon)
                c2 = val.monomial_coefficient(fmon).polynomial()
                for zpow in c2.monomials():
                    deg_z = zpow.degree()
                    coeff = c2.monomial_coefficient(zpow)
                    result[deg_f+deg_z] += coeff
        return vector(self.basering.base_ring(), result)

    def vector_to_element(self,v):
        result = self.ring(0)
        for i in range(v.length()):
            result += v[i] * self._basis[i]
        return result

class PrismaticEnvelopeG():
    def __init__(self,p,eisenstein, prec_F, prec_p, debug=False):
        self.debug=debug
        self.p = p
        self.eisenstein = eisenstein
        self.prec_F = prec_F
        self.prec_p = prec_p
        self._init_ring()
        self._init_nygaard_ring()
        self._init_lambda()
        self._init_R_terms()
        self._init_relations()
        self._reduce_relations_other()
        self._init_bk_unit()
    def _init_ring(self):
        if(self.debug):
            print('initializing ring')
        self.basering = self.eisenstein.parent()
        self.z = self.basering.gen()
        if(self.prec_F-1)==0:
            num_g=0
        else:
            num_g=floor(log(self.prec_F-1,self.p))+1
        if(num_g)==1:
            self.ring=PolynomialRing(self.basering,'g0',1)
        else:
            self.ring=PolynomialRing(self.basering,'g',num_g)
        self.g=self.ring.gens()
    def _init_nygaard_ring(self):
        if(self.debug):
            print('initializing nygaard ring')
        self.nygaard_ring = PolynomialRing(self.ring, 'dtilde')
        self.dtilde = self.nygaard_ring.gen()
    def weight_F(self,m):
        result=0
        for i in range(len(self.g)):
            result += self.p^i * m.degrees()[i]
        return result
    def weight_N(self,m):
        result=0
        for i in range(len(self.g)):
            result += self.p^i * m.degrees()[i]
        return result
    def trim(self,val):
        if(val.parent() == self.ring):
            return self._trim(val)
        if(val.parent() == self.nygaard_ring):
            return self._trim_nygaard(val)
    def _trim(self,val):
        result = self.ring(0)
        for m in val.monomials():
            w=self.weight_F(m)
            if w>=self.prec_F:
                continue
            new_c = val.monomial_coefficient(m).add_bigoh(self.prec_F-w)
            result += self.ring.from_base_ring(new_c)*m
        return result
    def _trim_nygaard(self,val):
        result = self.nygaard_ring(0)
        for m in val.monomials():
            c=val.monomial_coefficient(m)
            result += self._trim(c)*m
        return result
    def phi_base(self,val):
        #this should ideally be in the "delta-ring base class" at some point
        return val.V(self.p)
    def _divide_coefficients_base(self,val):
        assert(val.parent() == self.basering)
        modp = val.map_coefficients(lambda x: x % self.p)
        assert(modp==0),"_divide_coefficients_base {}".format(val)
        result = val.map_coefficients(lambda x: x // self.p)
        return result
    def delta_base(self,val):
        assert(val.parent() == self.basering)
        try:
            result = self._divide_coefficients_base(self.phi_base(val)-val^self.p)
        except:
            print("Error in delta_base, val={}".format(val))
        return result
    def phi(self,val):
        result = self.ring(0)
        for m in val.monomials():
            c = val.monomial_coefficient(m)
            result += self.phi_base(c) * m(self._phi_g)
        if(hasattr(self,'_relations')):
            result = self.reduce(result)
        else:    
            result = self.trim(result)
        return self._trim(result)
    def _divide_coefficients(self,val):
        assert(val.parent() == self.ring)
        modp = val.map_coefficients(lambda x: x.map_coefficients(lambda x: x % self.p))
        assert(modp==0),"_divide_coefficients {}".format(val)
        result = val.map_coefficients(lambda x: x.map_coefficients(lambda x: x // self.p))
        return result
    def delta(self,val):
        assert(val.parent() == self.ring)
        pdelta = (self.phi(val) - val^self.p)
        pdelta = self.trim(pdelta)
        try:
            result = self._divide_coefficients(pdelta)
        except:
            print("Error in delta, val={}".format(val))
            print("_phi_g={}".format(self._phi_g))
            print("_R_term={}".format(self._R_term))
        return result
    def _init_R_terms(self):
        if(self.debug):
            print('initializing R terms')
        if len(self.g) == 0:
            self._delta_g = []
            self._phi_g = []
            self._phi_divided_g = []
            return
        if len(self.g) == 1:
            self._delta_g = [self.ring(0)]
            self._phi_g = [self.ring(0)]
            self._phi_divided_g = [self.ring(0)]
            return
        self._delta_g = []
        self._phi_g = len(self.g) * [None]   #this is so that self.phi works in the recursion, but will complain if we use it in the wrong way.
        d = self.eisenstein
        first_delta = self.ring(0)
        for j in range(1,self.p):
            first_delta += -(binomial(self.p,j)//self.p) * self.z^j * (-self.g[0])^(self.p-j)
        if self.p==2:
            first_delta += -self.g[0]^2
        self._delta_g.append(first_delta)
        self._phi_g[0] = self.g[0]^self.p + self.p*first_delta
        self._R_term = [self.delta(self.g[0]) * (1 / self.delta_base(d))]
        self._phi_divided_g = [self._lambda[0]*self.g[1] + self._R_term[0]]
        d_pow = d
        for j in range(1,len(self.g)-1):
            d_pow = d_pow^self.p
            delta_d_pow = self.delta_base(d_pow)
            deltaR = self.delta(self._R_term[j-1])
            w1 = self.ring(0)
            for i in range(1,self.p):
                w1 += -(binomial(self.p,i) // self.p) * (self.g[j] * self._lambda[j-1])^i * self._R_term[j-1]^(self.p-i)
            self._R_term.append((1/self._lambda_denominator[j-1]) * (deltaR + w1))
            #print('j={}\ndpow={}\ndeltadpow={}'.format(j,d_pow,delta_d_pow))

            delta_g = (1+delta_d_pow*self._lambda[j])*self.g[j+1] + delta_d_pow * self._R_term[j]
            self._delta_g.append(delta_g)
            self._phi_g[j] = self.g[j]^self.p + self.p*delta_g
            self._phi_divided_g.append(self._lambda[j]*self.g[j+1] + self._R_term[j])
        self._R_term.append(self.ring(0))
        self._phi_g[len(self.g)-1] = self.ring(0)
        self._phi_divided_g.append(self.ring(0))
    def _init_relations(self):
        if(self.debug):
            print('initializing relations')
        if len(self.g) == 0:
            self._relations = []
            return
        if len(self.g) == 1:
            self._relations = [self.ring(0)]
            return
        d = self.eisenstein
        self._relations=[]
        dpow = d
        for j in range(len(self.g)-1):
            dpow = dpow^self.p
            factor = -self.p + dpow * self._lambda[j]
            new_rel = factor * self.g[j+1] + dpow * self._R_term[j]
            self._relations.append(new_rel)
        self._relations.append(self.ring(0))
    def _reduce_single_step(self,val):
        result=self.ring(0)
        for m in val.monomials():
            new_summand=val.monomial_coefficient(m)
            degs=m.degrees()
            for exponent in range(len(degs)):
                r=degs[exponent]//self.p
                new_summand*=self._relations[exponent]^r*self.g[exponent]^(degs[exponent]-self.p*r)
            result += new_summand
        return self.trim(result)
    def _reduce_relations(self):
        if(self.debug):
            print('reducing relations')
        for j in range(len(self._relations)):
            prev_val = self._relations[j]
            self._relations[j] = self._reduce_single_step(prev_val)
            while(prev_val != self._relations[j]):
                coeff = self._relations[j].monomial_coefficient(self.g[j]^self.p)
                if coeff != 0:
                    self._relations[j] = (1/(1-coeff)) * (self._relations[j] - coeff * self.g[j]^self.p)
                else:
                    prev_val = self._relations[j]
                    self._relations[j] = self._reduce_single_step(prev_val)
    def _reduce_relations_other(self):
        if(self.debug):
            print('reducing relations')
        done = False
        while not done:
            done = True
            for j in range(len(self._relations)):
                prev_val = self._relations[j]
                coeff = self._relations[j].monomial_coefficient(self.g[j]^self.p)
                if coeff != 0:
                    self._relations[j] = (1/(1-coeff)) * (self._relations[j] - coeff * self.g[j]^self.p)
                else:
                    self._relations[j] = self._reduce_single_step(prev_val)
                if(prev_val != self._relations[j]):
                    if(self.debug):
                        print("g{} not done".format(j))
                    done=False
                else:
                    if(self.debug):
                        print("g{} done".format(j))
    def reduce(self,val):
        prev_val=val
        result=self._reduce_single_step(prev_val)
        while(prev_val != result):
            prev_val = result
            result = self._reduce_single_step(prev_val)
        return result
    def _init_lambda(self):
        if(self.debug):
            print('initializing lambda')
        d = self.eisenstein
        self._lambda=[-1/self.delta_base(d)]
        self._lambda_denominator=[]
        d_pow = d
        for j in range(1,len(self.g)):
            d_pow = d_pow^self.p
            self._lambda_denominator.append(1 - self.delta_base(d_pow*self._lambda[j-1]))
            self._lambda.append(self._lambda[j-1]^self.p / self._lambda_denominator[j-1])
    def _init_bk_unit(self):
        if(len(self.g)==0):
            self.bk_unit = self.ring(1)
            return
        if(self.debug):
            print('initializing bk unit')

        phiu = 1 - self._phi_divided_g[0]
        previous = phiu
        phiv = self.reduce(phiu * self.phi(phiu))
        while(previous != phiv):
            previous = phiv
            phiv = self.reduce(phiu * self.phi(phiv))
        self.bk_unit = phiv
    def basis(self):
        result = []
        zpow = self.ring(1)
        for i in range(self.prec_F):
            result.append(zpow)
            zpow*=self.z
        return result
#    def right_unit(self,val):
#        assert val.parent() == self.basering
#        pol = val.polynomial()
#        result = self.reduce(pol(self.z - self.g[0]))
#        return result

        
        


