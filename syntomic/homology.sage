class Homology():
    def __init__(self,comp):
        self.complex = comp
        self._compute_homology()
    def _compute_homology(self):
        nz=self.complex.nonzero_degrees()
        if(len(nz)==0):
            self.orders={}
            self.representatives={}
            return
        self._complex_smith_form()
        deg_min=min(nz)
        deg_max=max(nz)
        self.orders={}
        self.representatives={}
        for i in range(deg_min,deg_max+1):
            self.orders[i] = self._homology_smith_form(i)
            self.representatives[i] = self._compute_homology_representatives(i)
    def _integral_inverse(self,f):
        _,u,v=f.smith_form()
        #ufv=id, so f=u^{-1}v^{-1}, so f^{-1}=vu
        return v*u
    def _complex_smith_form(self):
        # Returns a triple (D,f,g) where each differential is in ``Smith form''
        # and f is a map of chain complexes C -> D.
        # Assumes that the degree of the differential is +1.
        nz=self.complex.nonzero_degrees()
        deg_min=min(nz)
        deg_max=max(nz)
        # The number of differentials
        lt=deg_max-deg_min
        self._f_dict={}
        self._g_dict={}
        for i in range(deg_min,deg_max+1):
            self._f_dict[i]=identity_matrix(self.complex.free_module_rank(i))
            self._g_dict[i]=identity_matrix(self.complex.free_module_rank(i))
        d_dict={}
        d_ranks={}
        j=0
        S,U,V=self.complex.differential()[deg_max-(j+1)].smith_form()
        d_dict[deg_max-(j+1)]=S
        rank=0
        while rank < S.ncols() and rank < S.nrows() and S[rank,rank] != 0:
            rank += 1
        d_ranks[deg_max-(j+1)]=rank
        self._f_dict[deg_max-j]=U*self._f_dict[deg_max-j]
        self._f_dict[deg_max-(j+1)]=self._integral_inverse(V)*self._f_dict[deg_max-(j+1)]    
           #Inverting here coerces to p-adic FIELD, and then smith_form later on produces garbage.
           #So we use another smith call in _integral_inverse to determine the inverse correctly (this should also produce precision bounds correctly)
        self._g_dict[deg_max-j]=self._g_dict[deg_max-j]*self._integral_inverse(U)
        self._g_dict[deg_max-(j+1)]=self._g_dict[deg_max-(j+1)]*V    
        for j in range(1,lt):
            new_d=self._f_dict[deg_max-j]*self.complex.differential()[deg_max-(j+1)]
            col_offset=0
            row_offset=d_ranks[deg_max-j]
            col_num=self.complex.differential()[deg_max-(j+1)].ncols()
            row_num=d_dict[deg_max-j].ncols() - d_ranks[deg_max-j]
            new_d_sub=new_d.submatrix(row_offset,col_offset,row_num,col_num)
            S_sub,U_sub,V=new_d_sub.smith_form()
            S=block_matrix([[Matrix(row_offset,col_num)],[S_sub]])
            U=block_matrix([[identity_matrix(row_offset),0],[0,U_sub]])
            self._f_dict[deg_max-j]=U*self._f_dict[deg_max-j]
            self._f_dict[deg_max-(j+1)]=V^(-1)*self._f_dict[deg_max-(j+1)]
            self._g_dict[deg_max-j]=self._g_dict[deg_max-j]*U^(-1)
            self._g_dict[deg_max-(j+1)]=self._g_dict[deg_max-(j+1)]*V
            d_dict[deg_max-(j+1)]=S
            rank=0
            while rank < S_sub.ncols() and rank < S_sub.nrows() and S_sub[rank,rank] != 0:
                rank += 1
            d_ranks[deg_max-(j+1)]=rank
        self._complex_smith = ChainComplex(d_dict)
        self._d_ranks=d_ranks
    def _homology_smith_form(self,i):
        nz=self._complex_smith.nonzero_degrees()
        deg_min=min(nz)
        deg_max=max(nz)
        if i<deg_min or i>deg_max:
            return []
        elif i==deg_max:
            num_gens=self._complex_smith.differential()[i-1].nrows()
            offset=0
            l=[]
            for j in range(num_gens):
                if j<self._complex_smith.differential()[i-1].ncols() and j+offset<self._complex_smith.differential()[i-1].nrows():
                    div=self._complex_smith.differential()[i-1][j+offset,j]
                    if not div.is_unit():
                        l.append(ZZ(div))
                else:
                    l.append(0)
            return l
        elif i>deg_min:
            num_gens=self._complex_smith.differential()[i].ncols()-self._d_ranks[i]
            l=[]
            offset= self._d_ranks[i]
            for j in range(num_gens):
                if j<self._complex_smith.differential()[i-1].ncols() and j+offset<self._complex_smith.differential()[i-1].nrows():
                    div=self._complex_smith.differential()[i-1][j+offset,j]
                    if not div.is_unit():
                        l.append(ZZ(div))
                else:
                    l.append(0)
            return l
        elif i==deg_min:
            num_gens=self._complex_smith.differential()[i].ncols()-self._d_ranks[i]
            l=[0]*num_gens
            return l
    def _compute_homology_representatives(self,i):
        nz=self._complex_smith.nonzero_degrees()
        deg_min=min(nz)
        deg_max=max(nz)
        if i<deg_min or i>deg_max:
            return []
        h = self.orders[i]
        result = []
        for j in range(self._g_dict[i].ncols() - len(h), self._g_dict[i].ncols()):
            column = self._g_dict[i].column(j)
            result.append(column)
        return result
            
        return self._g_dict[i].submatrix(0,self._g_dict[i].ncols() - len(h))
    def compute_cycle(self,i,v):
        mat = self.compute_homology_representatives(i)
        return mat*v
    def normalize_cycle(self,i,v):
        if self.complex.differential()[i]*v != 0:
            print("Not a cycle!")
        h = self.orders[i]    
        mat = self._f_dict[i].submatrix(self._f_dict[i].nrows() - len(h), 0)
        return self.compute_cycle(i, mat*v)


class MorphismHomology():
    def __init__(self, hA, hB, F):
        self.hA = hA
        self.hB = hB
        self.F = F
        self.compute_morphisms()
    def compute_morphisms(self):
        nz=self.F.keys()
        deg_min=min(nz)
        deg_max=max(nz)
        self.hF = {}
        for i in range(deg_min, deg_max+1):
            new_F=self.hB._f_dict[i]*self.F[i]*self.hA._g_dict[i]
            x=len(self.hA.orders[i])
            y=len(self.hB.orders[i])
            self.hF[i] = new_F.submatrix(new_F.nrows()-y,new_F.ncols()-x)
