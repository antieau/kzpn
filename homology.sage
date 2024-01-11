class Homology():
    def __init__(self,C):
        self.complex = C
        self.compute_homology()
    def compute_homology(self):
        self.complex_smith_form()
        nz=self.complex.nonzero_degrees()
        deg_min=min(nz)
        deg_max=max(nz)
        self.homology_orders={}
        self.homology_representatives={}
        for i in range(deg_min,deg_max+1):
            self.homology_orders[i] = self.homology_smith_form(i)
            self.homology_representatives[i] = self.compute_homology_representatives(i)
        
    def complex_smith_form(self):
        # Returns a triple (D,f,g) where each differential is in ``Smith form''
        # and f is a map of chain complexes C -> D.
        # Assumes that the degree of the differential is +1.
        nz=self.complex.nonzero_degrees()
        deg_min=min(nz)
        deg_max=max(nz)
        # The number of differentials
        lt=deg_max-deg_min
        self.f_dict={}
        self.g_dict={}
        for i in range(deg_min,deg_max+1):
            self.f_dict[i]=identity_matrix(self.complex.free_module_rank(i))
            self.g_dict[i]=identity_matrix(self.complex.free_module_rank(i))
        D_dict={}
        j=0
        S,U,V=self.complex.differential()[deg_max-(j+1)].smith_form()
        D_dict[deg_max-(j+1)]=S
        self.f_dict[deg_max-j]=U*self.f_dict[deg_max-j]
        self.f_dict[deg_max-(j+1)]=V^(-1)*self.f_dict[deg_max-(j+1)]
        self.g_dict[deg_max-j]=self.g_dict[deg_max-j]*U^(-1)
        self.g_dict[deg_max-(j+1)]=self.g_dict[deg_max-(j+1)]*V    
        for j in range(1,lt):
            new_d=self.f_dict[deg_max-j]*self.complex.differential()[deg_max-(j+1)]
            col_offset=0
            row_offset=D_dict[deg_max-j].transpose().rank()
            col_num=self.complex.differential()[deg_max-(j+1)].dimensions()[1]
            row_num=D_dict[deg_max-j].transpose().nullity()
            new_d_sub=new_d.submatrix(row_offset,col_offset,row_num,col_num)
            S_sub,U_sub,V=new_d_sub.smith_form()
            S=block_matrix([[Matrix(row_offset,col_num)],[S_sub]])
            U=block_matrix([[identity_matrix(row_offset),0],[0,U_sub]])
            self.f_dict[deg_max-j]=U*self.f_dict[deg_max-j]
            self.f_dict[deg_max-(j+1)]=V^(-1)*self.f_dict[deg_max-(j+1)]
            self.g_dict[deg_max-j]=self.g_dict[deg_max-j]*U^(-1)
            self.g_dict[deg_max-(j+1)]=self.g_dict[deg_max-(j+1)]*V
            D_dict[deg_max-(j+1)]=S
        self.complex_smith = ChainComplex(D_dict)
    def homology_smith_form(self,i):
        nz=self.complex_smith.nonzero_degrees()
        deg_min=min(nz)
        deg_max=max(nz)
        if i<deg_min or i>deg_max:
            return []
        elif i==deg_max:
            num_gens=self.complex_smith.differential()[i-1].nrows()
            offset=0
            l=[]
            for j in range(num_gens):
                if j<self.complex_smith.differential()[i-1].ncols() and j+offset<self.complex_smith.differential()[i-1].nrows():
                    div=self.complex_smith.differential()[i-1][j+offset,j]
                    if not div.is_unit():
                        l.append(div)
                else:
                    l.append(0)
            return l
        elif i>deg_min:
            num_gens=self.complex_smith.differential()[i].ncols()-self.complex_smith.differential()[i].rank()
            l=[]
            offset=self.complex_smith.differential()[i].rank()
            for j in range(num_gens):
                if j<self.complex_smith.differential()[i-1].ncols() and j+offset<self.complex_smith.differential()[i-1].nrows():
                    div=self.complex_smith.differential()[i-1][j+offset,j]
                    if not div.is_unit():
                        l.append(div)
                else:
                    l.append(0)
            return l
        elif i==deg_min:
            num_gens=self.complex_smith.differential()[i].ncols()-self.complex_smith.differential()[i].rank()
            l=[0]*num_gens
            return l
    def compute_homology_representatives(self,i):
        nz=self.complex_smith.nonzero_degrees()
        deg_min=min(nz)
        deg_max=max(nz)
        if i<deg_min or i>deg_max:
            return []
        h = self.homology_orders[i]
        return self.g_dict[i].submatrix(0,self.g_dict[i].ncols() - len(h))

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
            new_F=self.hB.f_dict[i]*self.F[i]*self.hA.g_dict[i]
            x=len(self.hA.homology_orders[i])
            y=len(self.hB.homology_orders[i])
            self.hF[i] = new_F.submatrix(new_F.nrows()-y,new_F.ncols()-x)