from syntomic.sage_import import sage_import

sage_import('syntomic/prismatic_envelope', fromlist=['PrismaticEnvelopeF', 'PrismaticEnvelopeG'])
sage_import('syntomic/matrices', fromlist=['divide_left', 'divide_right'])
sage_import('syntomic/homology', fromlist=['Homology'])



#To compute the syntomic complex Z_p(i),
#we need to know nablaN(i), nablaP(i), syn0(i), syn1(i) to Fprec=i*n
#To compute the THH representatives, we need Fprec=(i+1)*n
#To compute v1*k: Z_p(i-k(p-1)) -> Z_p(i), we need to know all up to Fprec=i*n



class SyntomicComplex():
    def __init__(self,prismaticEnvelopeF, prismaticEnvelopeG,debug=False):
        self.debug = debug
        self.prismaticEnvelopeF = prismaticEnvelopeF
        self.prismaticEnvelopeG = prismaticEnvelopeG
        self.basering = prismaticEnvelopeF.basering.base_ring()
    def _compute_can0(self,i,prec_F):
        if(self.debug):
            print('computing can0')
        if(prec_F > self.prismaticEnvelopeF.prec_F):
            raise ValueError('Not enough prec_F')
        self._can0=Matrix(self.basering, prec_F-1, prec_F-1)
        basis = self.prismaticEnvelopeF.nygaard_basis(i)
        for j in range(1,prec_F):
            can_img = self.prismaticEnvelopeF.can(basis[j])
            column = self.prismaticEnvelopeF.element_to_vector(can_img)
            for k in range(1,prec_F):
                self._can0[k-1,j-1]=column[k]
    def _compute_can1(self,i,prec_F):            
        if(self.debug):
            print('computing can1')
        if(prec_F > self.prismaticEnvelopeF.prec_F):
            raise ValueError('Not enough prec_F')
        self._can1=Matrix(self.basering, prec_F-1, prec_F-1)
        basis = self.prismaticEnvelopeF.nygaard_basis(i-1)
        for j in range(0,prec_F-1):
            can_img = self.prismaticEnvelopeF.can(basis[j])
            column = self.prismaticEnvelopeF.element_to_vector(can_img)
            for k in range(0,prec_F-1):
                self._can1[k,j]=column[k]
    def _compute_phi0(self,i,prec_F):        
        if(self.debug):
            print('computing phi0')
        if(prec_F > self.prismaticEnvelopeF.prec_F):
            raise ValueError('Not enough prec_F')
        self._phi0=Matrix(self.basering, prec_F-1, prec_F-1)
        basis = self.prismaticEnvelopeF.nygaard_basis(i)
        for j in range(1,prec_F):
            phi_img = self.prismaticEnvelopeF.phi_divided(i, basis[j])
            column = self.prismaticEnvelopeF.element_to_vector(phi_img)
            for k in range(1,prec_F):
                self._phi0[k-1,j-1]=column[k]
    def _compute_nablaPOK(self,i,prec_F):
        if(self.debug):
            print('computing nablaPOK')
        if(prec_F > self.prismaticEnvelopeG.prec_F):
            raise ValueError('Not enough prec_F')
        self._nablaPOK=Matrix(self.basering, prec_F-1, prec_F-1)
        z = self.prismaticEnvelopeG.z

        bk_unit_exponent = i
        bk_unit_pow = self.prismaticEnvelopeG.ring(1)
        bk_unit_base = self.prismaticEnvelopeG.bk_unit
        while(bk_unit_exponent > 0):
            if bk_unit_exponent % 2 != 0:
                bk_unit_pow = self.prismaticEnvelopeG.reduce(bk_unit_pow * bk_unit_base)
            bk_unit_exponent = bk_unit_exponent // 2    
            bk_unit_base = self.prismaticEnvelopeG.reduce(bk_unit_base * bk_unit_base)
        previous_term = -bk_unit_pow
        zR = self.prismaticEnvelopeG.z - self.prismaticEnvelopeG.g[0]
        for j in range(1,prec_F):
            previous_term = self.prismaticEnvelopeG.reduce(zR * previous_term)
            column = previous_term.monomial_coefficient(self.prismaticEnvelopeG.g[0])
            for k in range(1,prec_F):
                self._nablaPOK[k-1,j-1]=column[k-1]
    def _compute_red0(self,prec_F):
        if(self.debug):
            print('computing red0')
        if(prec_F > self.prismaticEnvelopeF.prec_F):
            raise ValueError('Not enough prec_F')
        self._red0=Matrix(self.basering, prec_F-1, prec_F-1)
        z = self.prismaticEnvelopeF.ring.from_base_ring(self.prismaticEnvelopeF.z)
        basis = [z^i for i in range(0,prec_F)]
        for j in range(1,prec_F):
            can_img = self.prismaticEnvelopeF.reduce(basis[j])
            column = self.prismaticEnvelopeF.element_to_vector(can_img)
            for k in range(1,prec_F):
                self._red0[k-1,j-1]=column[k]
    def _compute_red1(self,prec_F):
        if(self.debug):
            print('computing red1')
        if(prec_F > self.prismaticEnvelopeF.prec_F):
            raise ValueError('Not enough prec_F')
        self._red1=Matrix(self.basering, prec_F-1, prec_F-1)
        z = self.prismaticEnvelopeF.ring.from_base_ring(self.prismaticEnvelopeF.z)
        basis = [z^i for i in range(0,prec_F)]
        for j in range(0,prec_F-1):
            can_img = self.prismaticEnvelopeF.reduce(basis[j])
            column = self.prismaticEnvelopeF.element_to_vector(can_img)
            for k in range(0,prec_F-1):
                self._red1[k,j]=column[k]
    def _compute_nablaP(self,prec_F):
        if(self.debug):
            print('computing nablaP')
        #nablaP * red0 = red1 * nablaPOK
        self._nablaP = divide_right(self._red1 * self._nablaPOK, self._red0)
    def _compute_nablaN(self,prec_F):
        if(self.debug):
            print('computing nablaN')
        #nablaP * can0 = can1 * nablaN
        self._nablaN = divide_left(self._can1, self._nablaP * self._can0)
    def _compute_phi1(self,prec_F):
        if(self.debug):
            print('computing phi1')
        #nablaP * phi0 = phi1 * nablaN
        self._phi1 = divide_right(self._nablaP*self._phi0, self._nablaN)
    def _compute_syn0(self,prec_F):
        if(self.debug):
            print('computing syn0')
        self._syn0 = self._can0 - self._phi0
    def _compute_syn1(self,prec_F):
        if(self.debug):
            print('computing syn1')
        self._syn1 = self._can1 - self._phi1
    def _compute_matrices(self,i,prec_F):
        if(self.debug):
            print('computing matrices')
        self._compute_can0(i,prec_F)
        self._compute_can1(i,prec_F)
        self._compute_phi0(i,prec_F)
        self._compute_red0(prec_F)
        self._compute_red1(prec_F)
        self._compute_nablaPOK(i,prec_F)
        self._compute_nablaP(prec_F)
        self._compute_nablaN(prec_F)
        self._compute_phi1(prec_F)
        self._compute_syn0(prec_F)
        self._compute_syn1(prec_F)
        self.max_prec = prec_F
    def print_matrices(self,prec):
        print("can0:\n{}\n\n".format(self._can0.apply_map(lambda x: x.add_bigoh(prec))))
        print("can1:\n{}\n\n".format(self._can1.apply_map(lambda x: x.add_bigoh(prec))))
        print("phi0:\n{}\n\n".format(self._phi0.apply_map(lambda x: x.add_bigoh(prec))))
        print("red0:\n{}\n\n".format(self._red0.apply_map(lambda x: x.add_bigoh(prec))))
        print("red1:\n{}\n\n".format(self._red1.apply_map(lambda x: x.add_bigoh(prec))))
        print("nablaPOK:\n{}\n\n".format(self._nablaPOK.apply_map(lambda x: x.add_bigoh(prec))))
        print("nablaP:\n{}\n\n".format(self._nablaP.apply_map(lambda x: x.add_bigoh(prec))))
        print("nablaN:\n{}\n\n".format(self._nablaN.apply_map(lambda x: x.add_bigoh(prec))))
        print("phi1:\n{}\n\n".format(self._phi1.apply_map(lambda x: x.add_bigoh(prec))))
        print("syn0:\n{}\n\n".format(self._syn0.apply_map(lambda x: x.add_bigoh(prec))))
        print("syn1:\n{}\n\n".format(self._syn1.apply_map(lambda x: x.add_bigoh(prec))))
    def _compute_chain_complex(self,prec_F):
        d0=block_matrix([[self._syn0],[self._nablaN]])
        d1=block_matrix([[self._nablaP, -self._syn1]])
        self.complex=ChainComplex({0:d0,1:d1})
    def _compute_chain_complex_mod_p(self,prec_F):
        d0=matrix(GF(p),block_matrix([[self._syn0],[self._nablaN]]))
        d1=matrix(GF(p),block_matrix([[self._nablaP, -self._syn1]]))
        self.complex_mod_p=ChainComplex({0:d0,1:d1})
    def _compute_homology(self):
        self.homology=Homology(self.complex)
    def _compute_homology_mod_p(self):
        self.homology_mod_p=Homology(self.complex_mod_p)
