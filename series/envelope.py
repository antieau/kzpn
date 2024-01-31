from series.series import WeightedPowerSeriesRingCapped
from series.series import WeightedPowerSeriesRingCappedHomomorphism
from series.series import WeightedPowerSeriesRingCappedElement
from series.delta import DeltaPowerSeriesCapped
from series.delta import DeltaPowerSeriesCappedHomomorphism
from series.prism import Prism
from series.prism import PrismHomomorphism

class FrobeniusTwistedPrismaticEnvelope:
    """
    For computing prismatic envelopes for regular sequence quotients.
    """
    def __init__(self,prism,relations):
        self._domain=prism
        self._relations=relations

        # Make new power series ring.

        new_names = self._domain.underlying_ring()._names
        new_weights = self._domain.underlying_ring()._weights
        new_ngens = self._domain.underlying_ring()._ngens
        max_relation_number = 0
        num_f_array = []
        for a in range(len(self._relations)):
            d = relations[a].filtration_weight()
            num_f=floor(log(self._domain.precision_cap()/d,p))+1
            num_f_array.append(num_f)
            max_relation_number=max(max_relation_number,num_f)
            new_weights += [d*self._domain.prime()**b for b in range(num_f+1)]
            new_ngens += num_f
            self.names+=["f"+"{}{}".format(b,a) for b in range(num_f+1)]

        codomain_underlying_ring = WeightedPowerSeriesRingCapped(self.domain.coefficient_ring(),self._domain.precision_cap(),new_ngens,new_weights,names=new_names)

        insertion_list=codomain.underlying_ring.gens()[0:len(self._domain.underlying_ring()._ngens)]

        # Make new power series ring homomorphism.

        underlying_ring_homomorphism = WeightedPowerSeriesRingCappedHomomorphism(self._domain.underlying_ring(),codomain_underlying_ring,insertion_list)
        codomain_distinguished_element = underlying_ring_homomorphism(self._domain.distinguished_element())

        # Initialize deltas, Frobenii, units, and envelope_relations.

        ### Initialize lambda_u
        domain_lambda_units=[]
        domain_lambda_units.append(-(self._domain.underlying_delta_ring().delta(self._domain.distinguished_element()).inverse()))
        # For f_0,...,f_r we need lambda_0,...,lambda_{r-1}.
        for u in range(max_relation_number-1):
            domain_lambda_units.append(
                domain_lambda_units[u]**p/(self._domain.underlying_ring().one()-self.underlying_delta_ring().delta(domain_lambda_units[u]*(self.distinguished_element()**(p**(u+1)))))
            )
        codomain_lambda_units=[underlying_ring_homomorphism(u) for u in domain_lambda_units]

        ### Initialize relations.
        codomain_relations = [underlying_ring_homomorphism(r) for r in self._relations]
        codomain_Rs = []
        For a in range(len(self._relations)):
            # The first R.
            codomain_Rs.append(codomain)

            for b in range(num_f_array[a]-1):

            # Now append zero for the top R.
            codomain_Rs.append(codomain_underlying_ring.zero())


        # Make new delta ring with new variables.

        # Make new delta ring homomorphism.

        # Make new prism

        # Make new prism homomorphism
