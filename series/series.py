"""
This is a module for absolutely capped precision power series in multiple variables with arbitrary
positive integer weights.

To model k[[x_1,...,x_n]] we include it into k[y_1,...,y_n][[T]] x_i maps to y_i*T^w_i
where w_i is the weight of x_i. Under this inclusion, the power series ring is closed under all
arithmetic operations in the ambient power series ring.

Exposed classes:
 - WeightedPowerSeriesRingCapped,
 - WeightedPowerSeriesRingCappedHomomorphism, and
 - WeightedPowerSeriesRingCappedElement.

"""

from sage.all import *

class WeightedPowerSeriesRingCapped:
    def __init__(
        self,
        coefficient_ring,
        precision_cap,
        num_vars,
        weights,
        prefix="z",
        names=None,
    ):
        """
        Returns a new WeightedPowerSeriesRingCapped.

        coefficient_ring - a commutative ring
        precision_cap - an integer giving the precision for the prefix-adic coefficients
        num_vars - the number of variables
        weights - a tuple of *positive* integer weights, one for each variable
        prefix - a string giving the prefix for the variables
        names - a list of variable names; this overrides prefix if given

        """
        if coefficient_ring.is_ring() == False:
            raise ValueError("Argument 'coefficient_ring' must be a ring.")
        if precision_cap.is_integer() == False or precision_cap<=0:
            raise ValueError("Argument 'precision_cap' must be a positive integer.")
        self._coefficient_ring = coefficient_ring
        self._precision_cap = precision_cap
        self._num_vars = num_vars
        self._prefix = prefix
        if names:
            self._variable_names = names
        else:
            self._variable_names = [
                prefix + "{}".format(t) for t in range(self._num_vars)
            ]
        self._weights = tuple(weights)
        self._polynomial_ring=PolynomialRing(self._coefficient_ring,self._prefix,self._num_vars,order=TermOrder('wdeglex',self._weights))
        self._element_class = WeightedPowerSeriesRingCappedElement

    def coefficient_ring(self):
        return self._coefficient_ring

    def precision_cap(self):
        return self._precision_cap

    def is_ring(self):
        return True

    def is_commutative(self):
        return self._coefficient_ring.is_commutative()

    def __str__(self):
        return "Absolute precision power series ring in {} over {} with weights {} and with absolutely capped precision {}".format(
            ", ".join(self._variable_names),
            self._coefficient_ring,
            self._weights,
            self._precision_cap,
        )

    def __repr__(self):
        return self.__str__()

    def __call__(self, x):
        """
        Transforms an input polynomial into the internal representation.

        """
        pass

    def gens(self):
        self._gens=[]
        for generator in self._polynomial_ring.gens():
            # We have to coerce the SAGE integer into a Python integer.
            self._gens.append(self._element_class(self,[(int(generator.degree()),generator)]))
        return self._gens

    def zero(self):
        return self._element_class(self,[])

    def one(self):
        return self._element_class(self,[(0,self._polynomial_ring.one())])


class WeightedPowerSeriesRingCappedHomomorphism:
    def __init__(self, domain, codomain, action_on_generators):
        # A list of elements of codomian giving the image of each
        # generator of domain.
        if domain.coefficient_ring() != codomain.coefficient_ring():
            raise ValueError(
                "Homomorphisms require the coefficient rings to be identical."
            )
        if len(action_on_generators) != domain._num_vars:
            raise TypeError(
                "Cannot interpret {} as a homomorphism from {} to {}".format(
                    action_on_generators, domain, codomain
                )
            )
        if codomain.precision_cap()>domain.precision_cap():
            raise ValueError("Codomain precision cap is larger than domain precision cap.")
        for i in range(len(action_on_generators)):
            if action_on_generators[i].parent() != codomain:
                raise ValueError("The elements of 'action_on_generators' are not members of the codomain.")
            if not action_on_generators[i].is_homogeneous() or (action_on_generators[i].is_homogeneous() and action_on_generators[i].degree() != domain._weights[i]):
                # TODO: allow filtered and not just graded morphisms.
                # This done, __call__ below must be rewritten.
                raise ValueError("Generators should be sent to homogeneous polynomials of the same degree.")

        self._domain = domain
        self._codomain = codomain
        self._action_on_generators = action_on_generators
        self._polynomial_action_on_generators = [g._underlying_polynomial() for g in self._action_on_generators]
        self._polynomial_homomorphism=self._domain._polynomial_ring.hom(self._polynomial_action_on_generators,self._codomain._polynomial_ring)

    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def __call__(self, f):
        new_term_list=[]
        for deg,coeff in f._term_list:
            new_term_list.append((deg,self._polynomial_homomorphism(coeff)))
        return self._codomain._element_class(self._codomain,new_term_list)

    def __str__(self):
        return "Homomorphism from {} to {} defined by {} on generators".format(
            self._domain, self._codomain, self._action_on_generators
        )

    def __repr__(self):
        return self.__str__()



class WeightedPowerSeriesRingCappedElement:
    def __init__(self, parent, term_list):
        """
        The argument term_list is a list [(N,coefficient)] representing g*T^N
        where g is a homogoneous polynomial of degree N.
        """
        for term in term_list:
            if term[0] != term[1].degree():
                raise ValueError("Input 'term_list' must be a list of pairs (N,g) where n is a non-negative integer and g is a homogeneous polynomial of degree N.")
            # We do no automatic coercions.
            if term[1].parent() != parent._polynomial_ring:
                raise ValueError("Input coefficients must be members of parent.")
        self._term_list = term_list
        self._parent = parent

    def parent(self):
        return self._parent

    def __str__(self):
        """
        This will return a string representing the underlying element
        in whatever order the coefficient_dictionary has stored its keys in. In particular
        it is not guaranteed to be in lexicographic ordering or anything else useful.
        """
        if len(self._term_list) == 0:
            return str("0 + O(F^{})".format(self.parent().precision_cap()))
        else:
            return_str = ""
            for deg,coef in self._term_list:
                if coef != self.parent()._polynomial_ring.zero():
                    return_str += str(coef) + " + "
            if len(return_str) == 0:
                return_str = '0 + '
            return return_str + "O(F^{})".format(self.parent().precision_cap())

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        term_dict=dict(self._term_list)
        for deg,coeff in other._term_list:
            if deg in term_dict:
                term_dict[deg]+=coeff
            else:
                term_dict[deg]=coeff
        sorted_degrees=list(term_dict)
        sorted_degrees.sort()
        new_term_list=[]
        for deg in sorted_degrees:
            if not term_dict[deg].is_zero():
                new_term_list.append((deg,term_dict[deg]))
        return self.__class__(self.parent(),new_term_list)

    def __sub__(self, other):
        term_dict=dict(self._term_list)
        for deg,coeff in other._term_list:
            if deg in term_dict:
                term_dict[deg]-=coeff
            else:
                term_dict[deg]=-coeff
        sorted_degrees=list(term_dict)
        sorted_degrees.sort()
        new_term_list=[]
        for deg in sorted_degrees:
            if not term_dict[deg].is_zero():
                new_term_list.append((deg,term_dict[deg]))
        return self.__class__(self.parent(),new_term_list)

    def __neg__(self):
        return self.__class__(self.parent(),[(deg,-coeff) for deg,coeff in self._term_list])

    def __eq__(self, other):
        pass

    def __ne__(self, other):
        pass

    def is_unit(self):
        if self._term_list[0][0] != 0 or not self._term_list[0][1].constant_coefficient().is_unit():
            return False
        else:
            return True

    def inverse(self):
        newton_approximation_inverse = self.__class__(self.parent(),[(0,self.parent()._polynomial_ring(self._term_list[0][1].constant_coefficient().inverse()))])
        N=1
        while N < self.parent().precision_cap():
            N=2*N
            newton_approximation_inverse=newton_approximation_inverse + newton_approximation_inverse*(self.parent().one()-self*newton_approximation_inverse)
        return newton_approximation_inverse

    def __truediv__(self, other):
        """
        Divide self/other by an invertible power series using Newton's method.
        """
        if not other.is_unit():
            raise ValueError("Constant term of denominator must be a unit.")
        if self == self.parent().zero():
            return self.parent().zero()
        return self*other.inverse()

    def __rmul__(self,other):
        if other.parent() != self.parent().coefficient_ring():
            return NotImplemented
        else:
            # TODO: this could result in zero terms.
            return self.__class__(self.parent(),[(deg,other*coefficient) for deg,coefficient in self._term_list])

    def __mul__(self, other):
        """
        Multiplication assumes that the input term_lists are sorted by degree.
        """
        term_dict={}
        for deg,coeff in self._term_list:
            for deg1,coeff1 in other._term_list:
                if deg+deg1<self.parent().precision_cap():
                    if deg+deg1 in term_dict:
                        term_dict[deg+deg1]+=coeff*coeff1
                    else:
                        term_dict[deg+deg1]=coeff*coeff1
                else:
                    break
        sorted_degrees=list(term_dict)
        sorted_degrees.sort()
        new_term_list=[]
        for deg in sorted_degrees:
            if not term_dict[deg].is_zero():
                new_term_list.append((deg,term_dict[deg]))
        return self.__class__(self.parent(),new_term_list)


    def __pow__(self, n):
        if n == 0:
            return self.parent.one()
        else:
            # This is probably not very pythonic. But, it is the only
            # way I can figure out how to get self**n to work. The only other option
            # would be to implement __rpow__(self,n) and apply it n**self,
            # which would make for garbage code.
            # From sage/arith/power.pyx
            apow = self
            while not (n & 1): # While even...
                apow *= apow
                # Bitshift
                n >>= 1

            # Now multiply together the correct factors a^(2^i)
            res = apow
            n >>= 1
            while n:
                apow *= apow
                if n & 1:
                    res = apow * res
                n >>= 1
            return res

    def _naive_power(self,n):
        pass

    def is_homogeneous(self):
        """
        Constants and zero count as homogeneous.
        """
        non_zero_terms=0
        for deg,coeff in self._term_list:
            if not coeff.is_zero():
                non_zero_terms+=1
            if non_zero_terms>1:
                return False
        return True

    def degree(self):
        degree = -1
        for deg,coeff in self._term_list:
            if not coeff.is_zero():
                degree = deg
        return degree

    def _underlying_polynomial(self):
        out = self.parent()._polynomial_ring.zero()
        for deg,coeff in self._term_list:
            out+=coeff
        return out



if __name__ == "__main__":
    from sage.all import Zp, Zq
    from sage.rings.power_series_ring import PowerSeriesRing
    import doctest

    doctest.testmod()
