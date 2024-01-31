from series.series import WeightedPowerSeriesRingCapped
from series.series import WeightedPowerSeriesRingCappedHomomorphism
from series.series import WeightedPowerSeriesRingCappedElement
from series.delta import DeltaPowerSeriesCapped
from series.delta import DeltaPowerSeriesCappedHomomrphism

# TODO: rebase this to be a wrapper of a WeightedPowerSeriesRingCapped plus a Frobenius
# endomorphism plus delta.

class Prism():
    def __init__(self,underlying_delta_ring,distinguished_element):
        self._underlying_delta_ring=underlying_delta_ring
        self._distinguished_element=distinguished_element

    def __str__(self):
        return "Oriented prism given by the delta ring {} with distinguished element {}".format(
            self._underlying_delta_ring,
            self._distinguished_element
        )

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

class PrismHomomorphism():
    def __init__(self,domain,codomain,action_on_generators):
        pass
