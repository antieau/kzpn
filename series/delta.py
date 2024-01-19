from series.series import WeightedPowerSeriesRingCapped
from series.series import WeightedPowerSeriesRingCappedHomomorphism
from series.series import WeightedPowerSeriesRingCappedElement

# TODO: rebase this to be a wrapper of a WeightedPowerSeriesRingCapped plus a Frobenius
# endomorphism plus delta.

class DeltaPowerSeriesCapped(WeightedPowerSeriesRingCapped):
    def __init__(self,characteristic_prime,underlying_ring,deltas,frobenii=None):
        """
        Optionally give the frobenius values. This will be useful when some pre-reduction has been done.
        """
        try:
            underlying_ring.coefficient_ring().frobenius_endomorphism()
        except TypeError:
            raise TypeError("Base ring must implement Frobenius.")
        self._underlying_ring=underlying_ring
        self._element_class=DeltaPowerSeriesCappedElement
        self._prime=characteristic_prime
        self._deltas=deltas
        if frobenii:
            self._frobenii=frobenii
        else:
            self._frobenii=[]
            for i in range(len(self._underlying_ring.gens())):
                self._frobenii.append(self._underling_ring.gens()[i]**p+p*self._deltas[i])

    def __str__(self):
        return "Delta ring structure on {} with delta given by {} on the generators".format(
            self._underlying_ring,
            self._deltas
        )

    def self.delta():
        pass

    def frobenius_endomorphism(self,element):
        if not element.parent()==self:
            raise ValueError("Element must be an element of self.")
        else:
            new_terms=[]
            for deg,coeff in element._term_list:
                if self._prime*deg < self.precision_cap():
                    new_homogeneous=self._polynomial_ring.zero()
                    for term in coeff.monomials():
                        for m in term.monomials():
                            new_homogeneous+=self.coefficient_ring().frobenius_endomorphism()(coeff.monomial_coefficient(m))*(m**self._prime)
                    new_terms.append((self._prime*deg,new_homogeneous))
                else:
                    break
            return self._element_class(self,new_terms)

class DeltaPowerSeriesCappedHomomorphism(WeightedPowerSeriesRingCappedHomomorphism):
    def __init__(self,domain,codomain,action_on_generators):
        """
        No check is made that this is compatible with the delta-ring structures.
        """
        WeightedPowerSeriesRingCappedHomomorphism.__init__(self,domain,codomain,action_on_generators)

class DeltaPowerSeriesCappedElement(WeightedPowerSeriesRingCappedElement):
    def __init__(self,parent,term_list):
        WeightedPowerSeriesRingCappedElement.__init__(self,parent,term_list)
