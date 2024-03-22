from syntomic.sage_import import sage_import

sage_import('syntomic/prismatic_envelope', fromlist=['PrismaticEnvelopeF', 'PrismaticEnvelopeG'])
sage_import('syntomic/syntomic_complex', fromlist=['SyntomicComplex'])
sage_import('syntomic/spectralsequence', fromlist=['FilteredComplex'])
sage_import('syntomic/name', fromlist=['name_vector'])

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
        self._syntomic_complexes = {}
        for i in range(0,self.max_weight+1):
            if(self.debug):
                print("computing syntomic complex in weight {}".format(i))
            syn = SyntomicComplex(self.prec_F,i, self.envelopeF, self.envelopeG, self.debug)
            self._syntomic_complexes[i]=syn
    def _compute_complexes_mod_p(self):
        if(self.debug):
            print("computing complexes mod p")
        for syn in self._syntomic_complexes.values():
            syn._compute_chain_complex_mod_p()
    def _compute_v1_maps(self):
        if(self.debug):
            print("computing v1 maps")
        self.v1maps = {}
        for i in range(0,self.max_weight-self.p+2):
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
        if(i==0):
            v1P0=Matrix(self.basering.residue_field(), nrows=self.prec_F-1,ncols=1)
            image = self.envelopeF.reduce(self.envelopeF.ring.from_base_ring(self.eisenstein^self.p))
            column = self.envelopeF.element_to_vector(image)
            for k in range(1,self.prec_F):
                v1P0[k-1,0]=column[k]
            return v1P0
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
        if(i==0):
            v1N0=Matrix(self.basering.residue_field(), nrows=self.prec_F-1,ncols=1)
            dtilde = self.envelopeF.dtilde
            image = self.envelopeF.nygaard_reduce(self.p-1,dtilde^self.p)
            column = self.envelopeF.nygaard_element_to_vector(self.p-1,image)
            for k in range(1,self.prec_F):
                v1N0[k-1,0]=column[k]
            return v1N0
        v1N0=Matrix(self.basering.residue_field(), self.prec_F-1, self.prec_F-1)
        basis = self.envelopeF.nygaard_basis(i)
        dtilde = self.envelopeF.dtilde
        for j in range(1,self.prec_F):
            image = self.envelopeF.nygaard_reduce(i+self.p-1, dtilde^self.p * basis[j])
            column = self.envelopeF.nygaard_element_to_vector(i+self.p-1, image)
            for k in range(1,self.prec_F):
                v1N0[k-1,j-1] = column[k]
        return v1N0

    def _compute_v1P1(self,i):
        if(self.debug):
            print("computing v1P1 in weight {}".format(i))
        if(i==0):
            return Matrix(self.basering.residue_field(),nrows=self.prec_F-1,ncols=0)
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
        if(i==0):
            return Matrix(self.basering.residue_field(),nrows=self.prec_F-1,ncols=0)
        v1N1=Matrix(self.basering.residue_field(), self.prec_F-1, self.prec_F-1)
        basis = self.envelopeF.nygaard_basis(i-1)
        dtilde = self.envelopeF.dtilde
        for j in range(0,self.prec_F-1):
            image = self.envelopeF.nygaard_reduce(i-1+self.p-1, dtilde^self.p * basis[j])
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
            print("initializing spectral sequence in weight {} mod {}".format(weight,self.p-1))
        complexes = {}
        maps = {}
        for i in range(weight, self.max_weight+1, self.p-1):
            index = (i-weight)//(self.p-1)
            if(i==0):
                complexes[-index] = ChainComplex(base_ring=self.basering.residue_field(), data={0: Matrix(self.basering.residue_field(), 1,1)})
                continue
            complexes[-index] = self._syntomic_complexes[i].complex_mod_p
        max_index = index
        for i in range(weight, self.max_weight-self.p+2, self.p-1):
            index = (i-weight)//(self.p-1)
            maps[-index] = self.v1maps[i]
        self.spectral_sequences[weight] = FilteredComplex(self.basering.residue_field(), complexes, maps, -max_index)
    def page(self,r,i,n):
        #reindexing by weight
        imod = i % (self.p-1)
        index = (i-imod) // (self.p-1)
        return self.spectral_sequences[imod].page(r, -index,n)
    def differential(self,r,i,n):
        imod = i % (self.p-1)
        index = (i-imod) // (self.p-1)
        return self.spectral_sequences[imod].differential(r, -index,n)
    def complex(self,i):
        imod = i % (self.p-1)
        index = (i-imod) // (self.p-1)
        return self.spectral_sequences[imod]._get_complex(-index)
    def name_element_e1_cycle(self, i, n, vect):
        #print(vect)
        cycle = self.page(1,i,n).representatives * vect
        #print(cycle)
        bottom_rank = self.complex(i).free_module_rank(n)
        top_rank = self.complex(i-self.p+1).free_module_rank(n+1)

        top_vector=cycle[:top_rank]
        bottom_vector=cycle[top_rank:]

        #print(bottom_rank)
        #print(top_rank)

        if(bottom_rank > 0):
            if i==0 and n==0:
                bottom_string = str(bottom_vector[0])
            elif i==0 and n==1:
                bottom_string = name_vector(['∂'], bottom_vector)
            else:
                bottom_string = self._syntomic_complexes[i].name_element(n,bottom_vector)
        else:
            bottom_string = ''
        if(top_rank > 0):
            if i-self.p+1==0 and n==0:
                top_string = str(top_vector[0])
            elif i-self.p+1==0 and n==1:
                top_string = name_element(['∂'], top_vector)
            else:
                top_string = self._syntomic_complexes[i-self.p+1].name_element(n+1,top_vector)
        else:
            top_string = ''
        if bottom_string == '' and top_string == '':
            return '0'
        if bottom_string == '':
            return "({})*σ".format(top_string)
        if top_string == '':
            return bottom_string
        return "{} + ({})*σ".format(bottom_string,top_string)
    def name_element_er_e1(self, r, i, n, vector):
        e1rep = self.page(r,i,n).representatives_e1 * vector
        names = [ "e[{},{},{}]".format(i,n,k) for k in range(len(e1rep))]
        return name_vector(names, e1rep)
    def e1_page_decoration(self, i, n): 
        result = {}
        result['i']=int(i)
        result['n']=int(n)
        result['basis']=[]
        result['zero']=False
        e1rank = len(self.page(1,i,n).orders)
        if(e1rank == 0):
            result['zero']=True
            return result
        for k in range(e1rank):
            basis_elt = {}
            basis_elt['name'] = "e[{},{},{}]".format(i,n,k)
            vect = vector(self.basering.residue_field(), e1rank)
            vect[k] = 1
            basis_elt['cycle'] = self.name_element_e1_cycle(i,n,vect)
            result['basis'].append(basis_elt)
        return result
    def er_page_decoration(self,r,i,n):
        result={}
        result['r']=int(r)
        result['i']=int(i)
        result['n']=int(n)
        result['basis']=[]
        result['zero']=False
        errank = len(self.page(r,i,n).orders)
        if(errank == 0):
            result['zero']=True
            return result
        for k in range(errank):
            basis_elt={}
            vect = vector(self.basering.residue_field(), errank)
            vect[k] = 1
            basis_elt['cycle']= self.name_element_er_e1(r,i,n,vect)
            result['basis'].append(basis_elt)
        return result
    def diff_r_decoration(self,r,i,n):
        result={}
        source_rank = len(self.page(r,i,n).orders)
        target_rank = len(self.page(r,i-r*(self.p-1),n+1).orders)
        result['length'] = int(r)
        result['source_weight']=int(i)
        result['target_weight']=int(i-r*(self.p-1))
        result['source_degree']=int(n)
        result['target_degree']=int(n+1)
        result['zero']=True
        result['basis_values']=[]
        for j in range(source_rank):
            diff_line = {}
            vect = vector(self.basering.residue_field(), source_rank)
            vect[j] = 1
            diff_line['from'] = self.name_element_er_e1(r,i,n,vect)
            hit_vect = self.differential(r,i,n) * vect
            if(hit_vect != 0):
                result['zero'] = False
            diff_line['to'] = self.name_element_er_e1(r,i-r*(self.p-1),n+1,hit_vect)
            result['basis_values'].append(diff_line)
        return result
    def all_decorations(self):
        result={}
        result['pages']={} #indexed by bidegree
        result['differentials']={}
        for i in range(self.max_weight+1):
            for n in range(-1,3):
                e1page = self.e1_page_decoration(i,n)
                if(e1page['zero']):
                    continue
                index="{},{}".format(i,n)
                result['pages'][index]={}
                result['pages'][index]['e1']=e1page
                result['pages'][index]['er']=[]
                for r in range(2,self.max_weight // (self.p-1)):
                    erpage = self.er_page_decoration(r,i,n)
                    result['pages'][index]['er'].append(erpage)
                result['differentials'][index]=[]
                for r in range(1, self.max_weight // (self.p-1)):
                    diff = self.diff_r_decoration(r,i,n)
                    if diff['zero']:
                        continue
                    result['differentials'][index].append(diff)
        return result


