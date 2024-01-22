from syntomic.sage_import import sage_import

sage_import('syntomic/matrices', fromlist=['integral_inverse'])

class Homology():
    def __init__(self,comp):
        self._base_ring = comp.base_ring()
        self.complex = comp
        self._compute_homology()
    def _compute_homology(self):
        self._complex_smith_form()
        self._compute_orders_and_maps()
    def _complex_smith_form(self):
        # Computes a triple (D,f,g) where each differential is in ``Smith form''
        # and f is a map of chain complexes C -> D.
        # Assumes that the degree of the differential is +1.
        nz=self.complex.nonzero_degrees()
        self._f_dict={}
        self._g_dict={}
        for i in nz:
            self._f_dict[i]=identity_matrix(self.complex.free_module_rank(i))
            self._g_dict[i]=identity_matrix(self.complex.free_module_rank(i))
        d_dict={}
        d_ranks={}
        for deg in reversed(nz):
            if deg+1 not in nz:
                d_ranks[deg]=0
                d_dict[deg] = matrix(self._base_ring, nrows=0, ncols=self.complex.free_module_rank(deg))
                continue
            transformed_diff = self._f_dict[deg+1] * self.complex.differential()[deg]
            diff_sub = transformed_diff.submatrix(d_ranks[deg+1],0)
            diag,sub_transform1,transform2 = diff_sub.smith_form()
            new_diff = block_matrix([[Matrix(self._base_ring, d_ranks[deg+1], diag.ncols())],[diag]])
            transform1 = block_matrix([[identity_matrix(self._base_ring, d_ranks[deg+1]),0],[0,sub_transform1]])
            self._f_dict[deg+1] = transform1 * self._f_dict[deg+1]
            self._f_dict[deg] = integral_inverse(transform2) * self._f_dict[deg]
            self._g_dict[deg+1] = self._g_dict[deg+1] * integral_inverse(transform1)
            self._g_dict[deg] = self._g_dict[deg] * transform2
            d_dict[deg] = new_diff
            rank = 0
            while rank < min(diag.ncols(), diag.nrows()) and diag[rank,rank] != 0:
                rank += 1
            d_ranks[deg] = rank
        self._complex_smith = ChainComplex(base_ring=self._base_ring, data=d_dict)
        self._d_ranks = d_ranks

    def _compute_orders_and_maps(self):
        self.orders = {}
        self.representatives = {}
        self.projections = {}
        for deg in self._complex_smith.nonzero_degrees():
            kernel_start_index = self._d_ranks[deg]
            kernel_size = self._complex_smith.free_module_rank(deg) - kernel_start_index
            orders = []
            for j in range(kernel_size):
                if j<self._complex_smith.differential()[deg-1].ncols():
                    entry = self._complex_smith.differential()[deg-1][j+self._d_ranks[deg],j]
                    if entry.is_unit():
                        kernel_start_index += 1
                        kernel_size -= 1
                    else:
                        orders.append(entry)
                else:
                    orders.append(0)
            #now, cycle representatives consist of basis vectors from kernel_start_index to kernel_start_index + kernel_size - 1
            self.orders[deg] = orders
            self.representatives[deg] = self._g_dict[deg].submatrix(0,kernel_start_index)
            self.projections[deg] = self._f_dict[deg].submatrix(kernel_start_index,0)
    def order(self,deg):
        if deg in self.orders:
            return self.orders[deg]
        else:
            return []
    def representative(self,deg):
        if deg in self.representatives:
            return self.representatives[deg]
        else:
            return matrix(self._base_ring, ncols=0,nrows=0)
    def projection(self,deg):
        if deg in self.projections:
            return self.projections[deg]
        else:
            return matrix(self._base_ring, ncols=0,nrows=0)

def morphism_homology(source_homology, target_homology, morphism, deg):
    #print(target_homology.projection(deg))
    #print(morphism[deg])
    #print(morphism[deg].dimensions())
    #print(source_homology.representative(deg))
    #print(source_homology.representative(deg).dimensions())
    return target_homology.projection(deg) * morphism[deg] * source_homology.representative(deg)

def morphisms_homology(source_homology, target_homology, morphism):
    result = {}
    for deg in morphism.keys():
        result[deg] = morphism_homology(source_homology, target_homology, morphism, deg)
    return result

#todo: revisit sparseness here (for example, MorphismHomology.morphism_homology should be defined whereever one of source_homology and target_homology is, not necessarily both
#class MorphismHomology():
#    def __init__(self, source_homology, target_homology, morphism):
#        self.source_homology = source_homology
#        self.target_homology = target_homology
#        self.morphism = morphism
#        self.compute_morphism()
#    def compute_morphism(self):
#        nz=self.morphism.keys()
#        deg_min=min(nz)
#        deg_max=max(nz)
#        self.morphism_homology = {}
#        for i in range(deg_min, deg_max+1):
#            new_F=self.target_homology._f_dict[i]*self.morphism[i]*self.source_homology._g_dict[i]
#            x=len(self.source_homology.orders[i])
#            y=len(self.target_homology.orders[i])
#            self.morphism_homology[i] = new_F.submatrix(new_F.nrows()-y,new_F.ncols()-x)
