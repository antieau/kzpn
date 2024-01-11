"""
The precision module.

It contains various helper functions needed by the main iPython notebook examples
and the main prismatic-envelope.sage module.
"""

import math


def braces(degree, exponent):
    """Returns exponent if exponent divides degree and degree otherwise."""
    if degree % exponent == 0:
        return exponent
    return degree


def epsilon(weight, degree, exponent):
    """Returns 1 if 1<=degree<=weight*exponent and exponent does not divide degree."""
    if math.floor(degree / exponent) < math.ceil(degree / exponent) <= weight:
        return 1
    return 0


def my_valuation(prime, degree):
    """Returns valuation if degree is divisible by prime**valuation."""
    if degree == 0:
        raise ValueError("Input is zero so my_valuation is +Infinity")
    valuation = 0
    while degree % prime == 0:
        valuation += 1
        degree = degree // prime
    return valuation


# Yuri Sulyma's approach


def prism_exponent(prime, weight, exponent):
    """Returns a bound on the exponent of the prismatic cohomology.

    If F_q is the finite field with prime**residual_degree elements,
    returns a bound on the exponent of F^{[1,weight*exponent-1]}Prism_R
    where R is a finite chain ring of length exponent and residue field F_q.
    """
    valuation = 0
    for degree in range(1, weight * exponent):
        valuation += my_valuation(prime, braces(degree, exponent))
    return valuation


def nygaard_exponent(prime, weight, exponent):
    """Returns a bound on the exponent of the prismatic cohomology.

    See the documentation for prism_exponent() for more details.
    """
    valuation = 0
    for degree in range(1, weight * exponent):
        valuation += epsilon(weight, degree, exponent) + my_valuation(
            prime, braces(degree, exponent)
        )
    return valuation


def rel_exponent(weight, exponent):
    """Returns an upper bound on the p-adic valuation of the exponent."""
    return math.floor(exponent * (weight + 1) * weight / 2 - weight)


def dz_exponent(weight, exponent):
    """Returns an upper bound on the p-adic valuation of the exponent."""
    return math.floor(exponent * weight * (weight - 1) / 2)


def matrix_precision(matrix):
    """Returns the minimum absolute precision of the entries in matrix."""
    return min(entry.precision_absolute() for entry in matrix.coefficients())


def matrix_valuation(matrix):
    """Returns the minimum valuation of the entries in matrix."""
    return min(entry.my_valuation() for entry in matrix.coefficients())

def target_total_precision(p,i,n):
    # The calculated F-precision needed to compute at this weight
    Fprec=n*i

    # The target precision
    target_precision=nygaard_exponent(p,i,n)

    ####################
    # Precision losses #
    ####################

    ### From \delta
    precision_loss_delta=math.floor(math.log(Fprec-1,p))

    ### From passing from OK to OK/pi^n
    precision_loss_quotient=0
    for q in range(p,i):
        precision_loss_quotient+=n*valuation(p,math.factorial(q))

    ### From lifting nabla to Nygaard
    precision_loss_nygaard=n*math.floor(i*(i-1)/2)

    ### From computing can-phi on primitives
    precision_loss_primitives=target_precision

    ### From renormalizing the Eisenstein polynomial
    precision_loss_from_Eisenstein=1

    ### Probably this precision can be taken to be lower since we will only need the
    ### Fp-coefficient calculation
    total_precision=(target_precision+precision_loss_delta
                     +precision_loss_quotient
                     +precision_loss_nygaard
                     +precision_loss_primitives
                     +precision_loss_from_Eisenstein) 
    
    return Fprec,total_precision