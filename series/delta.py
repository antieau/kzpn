"""
A module for delta-rings of the Breuil--Kisin variety.
"""

from series.series import WeightedPowerSeriesRingCappedHomomorphism


class DeltaPowerSeriesCapped:
    """
    A basic class for delta-rings based filtered power series rings.
    """
    def __init__(self, characteristic_prime, underlying_ring, deltas, frobenii=None):
        """
        Optionally give the frobenius values. This will be useful when some pre-reduction has been done.
        """
        try:
            underlying_ring.coefficient_ring().frobenius_endomorphism()
        except TypeError as e:
            raise TypeError("Base ring must implement Frobenius.") from e
        self._underlying_ring = underlying_ring
        self._prime = characteristic_prime
        self._prime_element = self._underlying_ring._unflatten(
            self._underlying_ring._polynomial_ring(self._prime)
        )
        self._deltas = deltas
        if frobenii:
            # No check to guarantee that these are compatible with the input deltas.
            self._frobenii = frobenii
        else:
            self._frobenii = []
            for i in range(len(self._underlying_ring.gens())):
                self._frobenii.append(
                    self._underlying_ring.gens()[i] ** self._prime
                    + self._prime_element * self._deltas[i]
                )

        self._frobenius = WeightedPowerSeriesRingCappedHomomorphism(
            self._underlying_ring,
            self._underlying_ring,
            self.coefficient_ring().frobenius_endomorphism(),
            self._frobenii,
        )

    def __str__(self):
        return f"Delta ring structure on {self._underlying_ring} with delta given by {self._deltas} on the generators"

    def __repr__(self):
        return self.__str__()

    def delta(self, element):
        """
        Each application of delta loses one bit of p-adic precision.
        """
        return (
            self.frobenius_endomorphism()(element) - element**self._prime
        ) // self._prime

    def frobenius_endomorphism(self):
        """
        Returns the *function*. Thus, it must be called as self.frobenius_endomorphism()(element).
        This is to maintain compatibility with SAGE.
        """
        return self._frobenius

    def w1(self, a, b):
        """
        Each application of w1 loses one bit of p-adic precision.
        """
        return (a**self._prime + b**self._prime - (a + b) ** self._prime) // self._prime

    def coefficient_ring(self):
        return self._underlying_ring.coefficient_ring()

    def precision_cap(self):
        return self._underlying_ring.precision_cap()

    def underlying_ring(self):
        return self._underlying_ring

    def prime(self):
        return self._prime

    def prime_element(self):
        return self._prime_element

    def ngens(self):
        return self.underlying_ring().ngens()


class DeltaPowerSeriesCappedHomomorphism:
    """
    A class for delta-ring homomorphisms.
    """
    def __init__(self, domain, codomain, underlying_ring_homomorphism):
        """
        No check is made that this is compatible with the delta-ring structures.
        """
        self._domain = domain
        self._codomain = codomain
        self._underlying_ring_homomorphism = underlying_ring_homomorphism

    def underlying_ring_homomorphism(self):
        return self._underlying_ring_homomorphism

    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def __call__(self, f):
        assert f.parent() == self.domain().underlying_ring()
        return self.underlying_ring_homomorphism()(f)
