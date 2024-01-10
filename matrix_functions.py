"""Module for various matrix helper functions."""

from sage.all import Matrix, lcm, identity_matrix, rank, diagonal_matrix

def matrix_exponent(M):
    """A helper function.

    Input: a square matrix of full rank.
    Output: the smallest exponent N such that if y is in
    the target space, Ny is in the image.
    Warning: this does not check that M has full rank.
    """
    el_divs = M.elementary_divisors()
    N = 1
    for d in el_divs:
        N = lcm(N, d)
    return N

def square_complete_top(can0, can1, nablaP):
    """Completes the top of a partially commutative diagram of matrices of full rank."""
    R = nablaP.base_ring()
    return can0, can1, Matrix(R, (can1 ** (-1)) * nablaP * can0), nablaP

def square_complete_right(syn0, nablaN, nablaP):
    """Completes the right of a partially commutative diagram of matrices of full rank."""
    R = nablaP.base_ring()
    return syn0, Matrix(R, nablaP * syn0 * (nablaN ** (-1))), nablaN, nablaP

def square_complete_bottom(OKtoOKmodpi0, OKtoOKmodpi1, nablaP_OK):
    """Completes the bottom of a partially commutative diagram of matrices of full rank."""
    R = nablaP_OK.base_ring()
    return (
        OKtoOKmodpi0,
        OKtoOKmodpi1,
        nablaP_OK,
        Matrix(R, OKtoOKmodpi1 * nablaP_OK * (OKtoOKmodpi0 ** (-1))),
    )

def standard_projection_matrix(a, b):
    """Returns the standard projection matrix of shape a<=b."""
    if a > b:
        raise TypeError("The input should be a pair (a,b) where a<=b.")
    M = identity_matrix(b)
    return M[0:a, :]

def saturation(M):
    """The saturation of a matrix.

    Input: an r-by-s matrix M, which can be defined
    over any ring which supports smith normal form.
    Output: matrices Mtilde, P, Q, where
    Mtilde is the saturation of M,
    P is the rank-by-r projection operator P such that P*M=Mtilde,
    Q is the inclusion operator such that Q*Mtilde=M.
    """
    R = M.base_ring()
    r = rank(M)
    target_dim = M.dimensions()[0]
    if r == target_dim:
        return (
            M,
            Matrix(R, diagonal_matrix([1 for i in range(r)])),
            Matrix(R, diagonal_matrix([1 for i in range(r)])),
        )
    D, S, T = M.smith_form()
    P = standard_projection_matrix(r, target_dim)
    Dtilde = P * D
    return (
        Matrix(R, Dtilde * (T ** (-1))),
        Matrix(R, P * S),
        Matrix(R, (S ** (-1)) * P.transpose()),
    )

def saturate_square(syn0, syn1, etaN, etaP):
    """The saturation of a square.

    Takes four matrices such that syn1*etaN=etaP*syn0,
    and returns four matrices syn0,syn1tilde,etaNtilde,etaPtilde,
    where etaNtilde and etaPtilde are the saturations of etaN and etaP,
    and where syn1tilde is the induced map between the saturations
    so that syn1tilde*etaNtilde=etaPtilde&syn0.
    """
    etaNtilde, etaN_P, etaN_Q = saturation(etaN)
    etaPtilde, etaP_P, etaP_Q = saturation(etaP)
    syn1tilde = etaP_P * syn1 * etaN_Q
    return syn0, syn1tilde, etaNtilde, etaPtilde
