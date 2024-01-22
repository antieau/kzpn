def _divide_rows(mat, diag):
    for i in range(mat.nrows()):
        for j in range(mat.ncols()):
            assert mat[i,j] % diag[i,i] == 0, "divisibility problem in _divide_rows, {} % {} (i={},j={})".format(mat[i,j],diag[i,i],i,j)
            mat[i,j] = mat[i,j] // diag[i,i]
    return mat

def _divide_cols(mat, diag):
    for i in range(mat.nrows()):
        for j in range(mat.ncols()):
            assert mat[i,j] % diag[j,j] == 0, "divisibility problem in _divide_cols, {} % {} (i={},j={})".format(mat[i,j],diag[j,j],i,j)
            mat[i,j] = mat[i,j] // diag[j,j]
    return mat

#computes mat0^{-1}*mat1, i.e. finds m with mat0*m=mat1
#as mat0.smith_form produces diag, l, r with l*mat0*r = diag,
#the above equation is equivalent to l*mat0*r*r^{-1}*m=l*mat1,
#so m = r * diag^{-1} * l * mat1, which we can compute by dividing rows
def divide_left(mat0,mat1):
    diag,l,r = mat0.smith_form()
    result = r*_divide_rows(l*mat1,diag)
    return result

#computes mat0*mat1^{-1}, i.e. finds m with m*mat1=mat0
#as mat1.smith_form produces diag, l, r with l*mat1*r = diag,
#the above equation is equivalent to m*l^{-1}*l*mat1*r=mat0*r,
#so m = mat0 * r * diag^{-1} * l, which we can compute by dividing rows
def divide_right(mat0,mat1):
    diag,l,r = mat1.smith_form()
    return _divide_cols(mat0*r,diag)*l

def integral_inverse(f):
    _,u,v=f.smith_form()
    #ufv=id, so f=u^{-1}v^{-1}, so f^{-1}=vu
    return v*u
