from series.series import WeightedPowerSeriesRingCapped
from series.series import WeightedPowerSeriesRingCappedHomomorphism
from series.series import WeightedPowerSeriesRingCappedElement
from series.delta import DeltaPowerSeriesCapped
from series.delta import DeltaPowerSeriesCappedHomomorphism


class Prism:
    def __init__(self, underlying_delta_ring, distinguished_element):
        assert distinguished_element.parent() == underlying_delta_ring.underlying_ring()
        assert underlying_delta_ring.delta(distinguished_element).is_unit()
        self._underlying_delta_ring = underlying_delta_ring
        self._distinguished_element = distinguished_element

    def __str__(self):
        return "Oriented prism given by the {} with distinguished element {}".format(
            self._underlying_delta_ring, self._distinguished_element
        )

    def __repr__(self):
        return self.__str__()

    def prime(self):
        return self.underlying_delta_ring().prime()

    def ngens(self):
        return self.underlying_ring().ngens()

    def prime_element(self):
        return self.underlying_delta_ring().prime_element()

    def underlying_delta_ring(self):
        return self._underlying_delta_ring

    def distinguished_element(self):
        return self._distinguished_element

    def underlying_ring(self):
        return self._underlying_delta_ring.underlying_ring()

    def coefficient_ring(self):
        return self._underlying_delta_ring.coefficient_ring()

    def precision_cap(self):
        return self.underlying_ring().precision_cap()


class PrismHomomorphism:
    def __init__(self, domain, codomain, underlying_delta_ring_homomorphism):
        self._domain = domain
        self._codomain = codomain
        self._underlying_delta_ring_homomorphism = underlying_delta_ring_homomorphism

    def underlying_delta_ring_homomorphism(self):
        return self._underlying_delta_ring_homomorphism

    def underlying_ring_homomorphism(self):
        return self._underlying_delta_ring_homomorphism._underlying_ring_homomorphism

    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def __call__(self, f):
        assert f.parent() == self.domain().underlying_ring()
        return self.underlying_ring_homomorphism()(f)
