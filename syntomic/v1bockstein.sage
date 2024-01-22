from syntomic.sage_import import sage_import

sage_import('syntomic/prismatic_envelope', fromlist=['PrismaticEnvelopeF', 'PrismaticEnvelopeG'])
sage_import('syntomic/syntomic_complex', fromlist=['SyntomicComplex'])
sage_import('syntomic/spectralsequence', fromlist=['FilteredComplex'])

class BocksteinV1():
    def __init__(self,p,n,eisenstein,prec_p,max_weight,debug=False):
        self.basering = eisenstein.parent().base_ring()
        self.debug=debug
        self.p=p
        self.n=n
        self.eisenstein=eisenstein
        self.prec_p=prec_p
        self.max_weight=max_weight
        self.prec_F = (max_weight+1)*n
        self._compute_envelopeF()
        self._compute_envelopeG()
        self._compute_syntomic_complexes()
        self._compute_complexes_mod_p()
        self._compute_v1_maps()
        self._initialize_spectral_sequences()
    def _compute_envelopeF(self):
        if(self.debug):
            print("computing envelopeF")
        self.envelopeF =  PrismaticEnvelopeF(self.p,self.n,self.eisenstein,self.prec_F,self.prec_p, self.debug)
    def _compute_envelopeG(self):
        if(self.debug):
            print("computing envelopeF")
        self.envelopeG =  PrismaticEnvelopeG(self.p,self.eisenstein,self.prec_F,self.prec_p, self.debug)
    def _compute_syntomic_complexes(self):
        if(self.debug):
            print("computing syntomic_complexes")
        self.syntomic_complexes = {}
        for i in range(0,max_weight+1):
            if(self.debug):
                print("computing syntomic complex in weight {}".format(i))
            syn = SyntomicComplex(self.prec_F,i, self.envelopeF, self.envelopeG, self.debug)
            self.syntomic_complexes[i]=syn
    def _compute_complexes_mod_p(self):
        if(self.debug):
            print("computing complexes mod p")
        for syn in self.syntomic_complexes.values():
            syn._compute_chain_complex_mod_p()
    def _compute_v1_maps(self):
        if(self.debug):
            print("computing v1 maps")
        self.v1maps = {}
        for i in range(0,self.max_weight-p+2):
            v1N0 = self._compute_v1N0(i)
            v1P0 = self._compute_v1P0(i)
            v1N1 = self._compute_v1N1(i)
            v1P1 = self._compute_v1P1(i)
            maps = {}
            maps[0] = v1N0
            maps[1] = block_matrix([[v1P0,0],[0,v1N1]])
            maps[2] = v1P1
            self.v1maps[i] = maps
    def _compute_v1P0(self,i):
        if(self.debug):
            print("computing v1P0 in weight {}".format(i))
        #weight 0 is special, todo: revisit that
        v1P0=Matrix(self.basering.residue_field(), self.prec_F-1, self.prec_F-1)
        basis = self.envelopeF.basis()
        for j in range(1,self.prec_F):
            image = self.envelopeF.reduce(self.eisenstein^self.p * basis[j])
            column = self.envelopeF.element_to_vector(image)
            for k in range(1,self.prec_F):
                v1P0[k-1,j-1]=column[k]
        return v1P0

    def _compute_v1N0(self,i):
        if(self.debug):
            print("computing v1N0 in weight {}".format(i))
        v1N0=Matrix(self.basering.residue_field(), self.prec_F-1, self.prec_F-1)
        basis = self.envelopeF.nygaard_basis(i)
        dtilde = self.envelopeF.dtilde
        for j in range(1,self.prec_F):
            image = self.envelopeF.nygaard_reduce(i+self.p-1, dtilde^p * basis[j])
            column = self.envelopeF.nygaard_element_to_vector(i+self.p-1, image)
            for k in range(1,self.prec_F):
                v1N0[k-1,j-1] = column[k]
        return v1N0

    def _compute_v1P1(self,i):
        if(self.debug):
            print("computing v1P1 in weight {}".format(i))
        v1P1=Matrix(self.basering.residue_field(), self.prec_F-1, self.prec_F-1)
        basis = self.envelopeF.basis()
        for j in range(0,self.prec_F-1):
            image = self.envelopeF.reduce(self.eisenstein^self.p * basis[j])
            column = self.envelopeF.element_to_vector(image)
            for k in range(0,self.prec_F-1):
                v1P1[k,j]=column[k]
        return v1P1

    def _compute_v1N1(self,i):
        if(self.debug):
            print("computing v1N1 in weight {}".format(i))
        v1N1=Matrix(self.basering.residue_field(), self.prec_F-1, self.prec_F-1)
        basis = self.envelopeF.nygaard_basis(i-1)
        dtilde = self.envelopeF.dtilde
        for j in range(0,self.prec_F-1):
            image = self.envelopeF.nygaard_reduce(i-1+self.p-1, dtilde^p * basis[j])
            column = self.envelopeF.nygaard_element_to_vector(i-1+self.p-1, image)
            for k in range(0,self.prec_F-1):
                v1N1[k,j] = column[k]
        return v1N1
    def _initialize_spectral_sequences(self):
        if(self.debug):
            print("initializing spectral sequences")
        self.spectral_sequences={}
        for i in range(self.p-1):
            self._initialize_spectral_sequence_weight_residue(i)
    def _initialize_spectral_sequence_weight_residue(self,weight):
        if(self.debug):
            print("initializing spectral sequence in weight {} mod {}".format(weight,p-1))
        complexes = {}
        maps = {}
        for i in range(weight, self.max_weight+1, p-1):
            index = (i-weight)//(p-1)
            complexes[-index] = self.syntomic_complexes[i].complex_mod_p
        for i in range(weight, self.max_weight-p+2, p-1):
            index = (i-weight)//(p-1)
            maps[-index] = self.v1maps[i]
        self.spectral_sequences[weight] = FilteredComplex(self.basering.residue_field(), complexes, maps)
        
    

