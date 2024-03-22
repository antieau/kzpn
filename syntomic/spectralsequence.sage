from syntomic.sage_import import sage_import

sage_import('syntomic/matrices', fromlist=['integral_inverse'])
sage_import('syntomic/homology', fromlist=['Homology','morphisms_homology'])

class FilteredComplex():
    def __init__(self, base_ring, complexes, maps, min_filt):
        self._base_ring = base_ring
        self.complexes = complexes
        self.maps = maps
        self.min_filt = min_filt
        max_filt = max(self.complexes.keys())
        self.last_page = max_filt - min_filt
        self._sanitize_maps() 
        self._strata={}
        self._pages={}
        self._differentials={}
        self._compute_pages()
        self._compute_differentials()
    def _zero_map(self,s,deg):
        nrows = self._get_complex(s-1).free_module_rank(deg)
        ncols = self._get_complex(s).free_module_rank(deg)
        self.maps[s][deg] = matrix(self._base_ring, nrows=nrows, ncols=ncols)
    def _sanitize_maps(self):
        for s in self.complexes:
            if not s in self.maps:
                self.maps[s] = {}
            if not s+1 in self.maps:
                self.maps[s+1] = {}
            for deg in self.complexes[s].nonzero_degrees():
                if not deg in self.maps[s]:
                    self._zero_map(s,deg)
                if not deg in self.maps[s+1]:
                    self._zero_map(s+1,deg)
    def _get_complex(self,filtration):
        filtration = max(filtration, self.min_filt)
        if filtration in self.complexes:
            result = self.complexes[filtration]
            return result
        else:
            result = ChainComplex(base_ring=self._base_ring)
            return result
    def _get_single_map(self,source_filtration):
        if source_filtration in self.maps:
            result = self.maps[source_filtration]
            return result
        else:
            result = {}
            source = self._get_complex(source_filtration)
            target = self._get_complex(source_filtration-1)
            for n in set(source.nonzero_degrees() + target.nonzero_degrees()):
                if n in source.nonzero_degrees():
                    source_rk = source.differential()[n].ncols()
                else:
                    source_rk = 0
                if n in target.nonzero_degrees():
                    target_rk = target.differential()[n].ncols()
                else:
                    target_rk = 0
                result[n] = matrix(self._base_ring, nrows=target_rk, ncols=source_rk)
            return result    
    def _get_map(self,source_filtration,target_filtration):
        if source_filtration<target_filtration:
            raise ValueError("_get_map: source_filtration needs to be bigger than target_filtration")
        source_filtration = max(source_filtration, self.min_filt)
        target_filtration = max(target_filtration, self.min_filt)
        if source_filtration==target_filtration:
            result = {}
            source = self._get_complex(source_filtration)
            for n in source.nonzero_degrees():
                result[n] = identity_matrix(self._base_ring,source.differential()[n].ncols())
            return result
        morphism=self._get_single_map(target_filtration+1)
        for i in range(target_filtration+2, source_filtration+1):
            morphism = compose_morphism(morphism,self._get_single_map(i))
        return morphism
    def _get_stratum(self,r,s):    
        bottom = s
        top = s+r
        bottom = max(bottom, self.min_filt)
        top = max(top, self.min_filt)
        target = self._get_complex(bottom)
        source = self._get_complex(top)
        morphism = self._get_map(top,bottom)
        return cofiber(source, target, morphism)
    def _compute_page(self,r,s):
        source = self._get_stratum(r,s)
        target = self._get_stratum(r,s-r+1)

        map_top = self._get_map(s+r,s+1)
        morphism = morphism_cofiber(map_top,self._get_map(s,s-r+1))
        source_homology = Homology(source)
        target_homology = Homology(target)

        morphism_h = morphisms_homology(source_homology,target_homology,morphism)

        #map to stratum(1,s):
        morphism_e1 = morphism_cofiber(map_top, self._get_map(s,s))

        page = {}
        #compute image
        for n,f in morphism_h.items():
            if n in source_homology.orders:
                source_orders = source_homology.orders[n]
            else:
                source_orders=[]
            if n in target_homology.orders:
                target_orders = target_homology.orders[n]
            else:
                target_orders=[]
            image = Image(self._base_ring, source_orders, target_orders, f)
            orders = image.orders
            projection = image.projection * source_homology.projection(n)
            representatives = source_homology.representative(n) * image.representatives

            if n in morphism_e1:
                morphism_e1n=morphism_e1[n]
            else:
                morphism_e1n=matrix(self._base_ring, 0, 0)
            if r>1:
                representatives_e1 = self.page(1,s,n).projection * morphism_e1n * representatives
            else:
                representatives_e1 = identity_matrix(self._base_ring, len(image.orders))
            page[n] = PageData(image.orders, representatives_e1, representatives, projection)
            #page[n]['representatives_e1']
            #to compute, compute map f from stratum(r,s) -> stratum(1,s)
            #then compute induced map on page by following
            #   e(1,s)[n].projection*f*e(r,s)[n].representatives
        return page
    def _compute_differential(self,r,s):
        #Want differential from F^s / F^{s+r} (stratum r,s) to F^{s+r}/F^{s+2r} (stratum r, s+r)
        source = self._get_stratum(r,s)
        target = self._get_stratum(r,s+r)
        page_source = self._compute_page(r,s)
        page_target = self._compute_page(r,s+r)
        #source_homology = Homology(source)
        #target_homology = Homology(target)
        differential = {}
        nz = set([n-1 for n in target.nonzero_degrees()] + list(source.nonzero_degrees()))
        for deg in nz:
            #differential is induced by map of complexes C^n(F^s / F^{s+r}) -> C^{n+1}(F^{s+r}/F^{s+2r}), (-1)^n*block_matrix([0,0],[id,0])
            cols0 = self._get_complex(s+r).free_module_rank(deg+1)
            rows0 = self._get_complex(s+2*r).free_module_rank(deg+2)
            cols1 = self._get_complex(s).free_module_rank(deg)
            rows1 = self._get_complex(s+r).free_module_rank(deg+1)
            assert(rows1 == cols0)
            delta = block_matrix([[matrix(self._base_ring,rows0,cols0),matrix(self._base_ring,rows0,cols1)],[identity_matrix(self._base_ring, rows1), matrix(self._base_ring, rows1, cols1)]])
            if(deg % 2 == 1):
                delta = -delta
            #delta_on_homology = target_homology.projection(deg+1) * delta * source_homology.representative(deg)
            #then we need induced map on image.
            
            if deg in page_source:
                rep_source = page_source[deg].representatives
            else:
                rep_source = matrix(self._base_ring, nrows=delta.ncols(),ncols=0)
            if deg+1 in page_target:
                proj_target = page_target[deg+1].projection
            else:
                proj_target = matrix(self._base_ring, nrows=0, ncols=delta.nrows())
            differential[deg] = proj_target*delta*rep_source
        return differential
    def _compute_pages(self):
        max_filt = max(self.complexes.keys())
        min_filt = self.min_filt
        for r in range(1,max_filt-min_filt+1):
            for s in range(min_filt, max_filt+1):
                self._pages[r,s] = self._compute_page(r,s)
    def _compute_differentials(self):
        max_filt = max(self.complexes.keys())
        min_filt = self.min_filt
        for r in range(1,max_filt-min_filt):
            for s in range(min_filt, max_filt+1):
                self._differentials[r,s] = self._compute_differential(r,s)
    def page(self,r,s,n=None):
        if(r < 1):
            raise ValueError("Page index must be >= 1")
        r = min(r, self.last_page)
        if (r,s) in self._pages:
            page = self._pages[r,s]
        else:
            page = {}
        if n is None:
            return page
        if n in page:
            return page[n]
        stratum = self._get_stratum(r,s)
        rank = stratum.free_module_rank(n) #this could be computed a bit cheaper but I don't think it matters
        if (r>1):
            e1_rank=len(self.page(1,s,n).orders)
        else:
            e1_rank=0
        trivial_page_data = PageData([], matrix(self._base_ring, nrows=e1_rank, ncols=0), matrix(self._base_ring, nrows=rank, ncols=0), matrix(self._base_ring, nrows=0, ncols=rank))
        return trivial_page_data
    def differential(self,r,s,n=None):
        if r < 1:
            raise ValueError("Page index must be >= 1")
        if (r,s) in self._differentials:
            diff = self._differentials[r,s]
            if n is None:
                return diff
            if n in diff:
                return diff[n]
            source = self.page(r,s)
            target = self.page(r,s+r)
            if n in source:
                source_rank=len(source[n].orders)
            else:
                source_rank=0
            if n+1 in target:
                target_rank=len(target[n+1].orders)
            else:
                target_rank=0
            return matrix(self._base_ring, nrows=target_rank, ncols=target_rank)
        if n is None:
            result = {}
            source = self.page(r,s)
            target = self.page(r,s+r)
            keys = set(list(source.keys())+[n-1 for n in target.keys()])
            for n in keys:
                if n in source:
                    source_rank=len(source[n].orders)
                else:
                    source_rank=0
                if n+1 in target:
                    target_rank=len(target[n+1].orders)
                else:
                    target_rank=0
                result[n] = matrix(self._base_ring, nrows=target_rank, ncols=target_rank)
            return result    
        else:
            source = self.page(r,s)
            target = self.page(r,s+r)
            if n in source:
                source_rank=len(source[n].orders)
            else:
                source_rank=0
            if n+1 in target:
                target_rank=len(target[n+1].orders)
            else:
                target_rank=0
            return matrix(self._base_ring, nrows=target_rank, ncols=target_rank)

class PageData():
    def __init__(self,orders,representatives_e1, representatives,projection):
        self.orders=orders
        self.representatives_e1=representatives_e1
        self.representatives=representatives
        self.projection=projection
    def __repr__(self):
        return("{{orders: {}, representatives_e1: {}, representatives: {}, projection: {}}}".format(self.orders,self.representatives_e1,self.representatives,self.projection))





    #def _compute_page(self,r,s): 
    #    #if not computed yet
    #    self._compute_stratum(r,s) 
    #    self._compute_stratum(r,s-r+1)
    #    strata_map = #chain complex morphism built by composing r-1 copies of self.maps
    #    strata_map_homology = #induced map on homology
    #    #now compute image and store whatever data we need from that...
    #def _compute_differential(self,r,s):
    #    _compute_page(self,r,s)
    #    _compute_page(self,r,s+r)
    #    connecting_homomorphism = #obtained from projection and inclusion on level of mapping cones
    #    #compute induced map on images. For that we probably need the description of the image as cokernel from above, and lift the map over the kernel?
    #def _compute_page(self,r,s):
    #    self.connecting_homomorphisms[r] = #at [s], connecting hom. F^{s} / F^{s+r} -> F^{s+r}/F^{s+2r}[1]
    #    self._strata_maps[r] = #at [s], map of F^s/F^{s+r} -> F^{s-r+1}/F^{s+1}
    #    self.page[r] = #at s, quotient of H_*F^s/F^{s+r} by kernel of _strata_maps[r][s]
    #    self.d[r] = #at s, map of quotients represented by induced map of connecting_homomorphisms[r][s]
    #    return

def cofiber(domain, target, morphism):
    diffs = {}
    nz_degs =  set([n-1 for n in domain.nonzero_degrees()] + list(target.nonzero_degrees()))
    for deg in nz_degs:
        #if not deg+1 in domain.nonzero_degrees():
        #    diffs[deg] = target.differential(deg)
        #    continue
        #if not deg in target.nonzero_degrees():
        #    diffs[deg] = -domain.differential(deg+1)
        #    continue
        if deg+1 in morphism:
            f=morphism[deg+1]
        else:
            f=matrix(domain.base_ring(), target.free_module_rank(deg+1), domain.free_module_rank(deg+1))
        diff = block_matrix([[-domain.differential(deg+1), 0],[-f, target.differential(deg)]])
        diffs[deg] = diff
    result = ChainComplex(base_ring=domain.base_ring(), data=diffs)
    return result

    #degs_domain = domain.nonzero_degrees() 
    #degs_target = target.nonzero_degrees() 
    #if len(degs_domain) == 0:
    #    return target
    #if len(degs_target) == 0:
    #    diffs = {}
    #    for deg, diff in domain.differential().items():
    #        diffs[deg-1] = -diff
    #    return ChainComplex(base_ring=self.base_ring, data=diffs)
    #min_deg = min(min(degs_domain)-1, min(degs_target))
    #max_deg = max(max(degs_domain)-1, max(degs_target))
    #diffs = {}
    #for i in range(min_deg, max_deg):
    #    if (not i+1 in domain.differential() and not i in target.differential()):
    #        continue
    #    if not i+1 in domain.differential():
    #        diffs[i] = target.differential()[i]
    #        continue
    #    if not i in target.differential():
    #        diffs[i] = -domain.differential()[i+1]
    #        continue
    #    diff = block_matrix([[-domain.differential()[i+1], 0],[-morphism[i+1], target.differential()[i]]])
    #    diffs[i] = diff
    #return ChainComplex(base_ring=domain.base_ring(),data=diffs)

def morphism_cofiber(domain_morphism, target_morphism):
    #min_deg = min(min(domain_morphism.keys())-1, min(target_morphism.keys()))
    #max_deg = max(max(domain_morphism.keys())-1, max(target_morphism.keys()))
    result = {}

    nz_degs = set([n-1 for n in domain_morphism.keys()] + list(target_morphism.keys()))
    for i in nz_degs:
        if not i+1 in domain_morphism:
            result[i] = target_morphism[i]
            continue
        if not i in target_morphism:
            result[i] = domain_morphism[i+1]
            continue
        result[i] = block_matrix([[domain_morphism[i+1],0],[0,target_morphism[i]]]) 
    return result

def compose_morphism(morphism1,morphism2):
    result = {}
    for i,f in morphism1.items():
        if i not in morphism2:
            result[i] = matrix(ring=f.base_ring(), nrows=f.nrows(),ncols=0)
        else:
            result[i] = f * morphism2[i]
    for i,f in morphism2.items():
        if i not in morphism1:
            result[i] = matrix(ring=f.base_ring(), nrows=0,ncols=f.ncols())
    return result

class Image():
    def __init__(self,base_ring, source_orders, target_orders, morphism):
        self.base_ring = base_ring
        self.source_orders = source_orders
        self.target_orders = target_orders
        self.morphism = morphism
        self._compute_image()
    def __str__(self):
        return "Image with orders " + str(self.orders)
    def __repr__(self):
        return "Image with orders " + str(self.orders)
    def _compute_image(self):
        target_orders_matrix = matrix(self.base_ring, len(self.target_orders), len(self.target_orders))
        for j in range(len(self.target_orders)):
            target_orders_matrix[j,j] = self.target_orders[j]
        joint_matrix = block_matrix([[self.morphism, target_orders_matrix]])
        S,U,V = joint_matrix.smith_form()
        nullspace_beginning_index = 0
        while nullspace_beginning_index < min(S.nrows(), S.ncols()) and S[nullspace_beginning_index,nullspace_beginning_index] != 0:
            nullspace_beginning_index += 1
        #now S = U*joint_matrix*V
        #kernel of S consists of basis vectors from nullspace_beginning_index to S.ncols()-1
        #so kernel of joint_matrix consists of columns of V from nullspace_beginning_index to S.ncols()-1
        #so kernel of source -> target is spanned by columns of the following submatrix of V:
        kernel_generator_matrix = V.submatrix(0, nullspace_beginning_index, self.morphism.ncols(), V.ncols()-nullspace_beginning_index)

        #finally, the image is computed as the cokernel of kernel_generator_matrix

        S2, U2, V2 = kernel_generator_matrix.smith_form()
        #have S2 = U2 * kernel_generator_matrix * V2, so cokernel is given by those rows of S2 which don't correspond to unit diagonal entries.
        #the projection from source to the image is given by the corresponding submatrix of U2.
        #similarly, the map taking elements of the image to representatives in source is given by the corresponding columns of inverse of U2

        U2_inverse = integral_inverse(U2)
        

        self.orders = []
        cokernel_beginning_index = 0
        for j in range(min(S2.nrows(), S2.ncols())):
            if(S2[j,j].is_unit()):
                cokernel_beginning_index += 1
            else:
                self.orders.append(S2[j,j])
        self.projection = U2.submatrix(cokernel_beginning_index,0)
        self.representatives = U2_inverse.submatrix(0,cokernel_beginning_index)
        
def morphism_image(image0,image1,morphism):
    return image1.projection * morphism * image0.representatives
        

