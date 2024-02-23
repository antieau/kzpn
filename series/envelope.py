"""
This module builds prismatic envelopes in some generality.
"""

from sage.all import floor, log, binomial

from series.series import WeightedPowerSeriesRingCapped
from series.series import WeightedPowerSeriesRingCappedHomomorphism
from series.delta import DeltaPowerSeriesCapped
from series.delta import DeltaPowerSeriesCappedHomomorphism
from series.prism import Prism
from series.prism import PrismHomomorphism


class FrobeniusTwistedPrismaticEnvelope:
    """
    For computing prismatic envelopes for regular sequence quotients.
    """

    def __init__(self, prism, relations):
        self._domain = prism
        self._relations = relations
        # It is assumed that no relation is zero and that no relation has zero filtration weight.
        self._prime = self._domain.prime()

        # Make new power series ring.
        new_names = self._domain.underlying_ring()._names.copy()
        new_weights = tuple()
        new_ngens = self._domain.underlying_ring()._ngens
        num_f_array = []
        for a in range(len(self._relations)):
            d = relations[a].filtration_weight()
            num_f = floor(log(self._domain.precision_cap() // d, self._prime))
            assert num_f > 0
            num_f_array.append(num_f)
            new_weights += tuple(d * self._domain.prime() ** b for b in range(num_f))
            new_ngens += num_f
            new_names += [f"f{b}{a}" for b in range(num_f)]

        print(f"new_names is {new_names}")
        print(f"new_weights is {new_weights}")
        codomain_underlying_ring = WeightedPowerSeriesRingCapped(
            self._domain.coefficient_ring(),
            self._domain.precision_cap(),
            new_ngens,
            self._domain.underlying_ring()._weights + new_weights,
            names=new_names,
            order=TermOrder(
                "wdeglex",
                self._domain.underlying_ring()._polynomial_ring.term_order().weights()
                + new_weights,
            ),
        )

        insertion_list = codomain_underlying_ring.gens()[
            0 : self._domain.underlying_ring()._ngens
        ]

        # Make new power series ring homomorphism.

        underlying_ring_homomorphism = WeightedPowerSeriesRingCappedHomomorphism(
            self._domain.underlying_ring(),
            codomain_underlying_ring,
            self._domain.coefficient_ring().hom(self._domain.coefficient_ring().gens()),
            insertion_list,
        )

        codomain_distinguished_element = underlying_ring_homomorphism(
            self._domain.distinguished_element()
        )

        # Initialize deltas, Frobenii, units, and envelope_relations.

        ### Create partial Frobenius and delta operations.

        def w1(a, b):
            out = codomain_underlying_ring.zero()
            for power in range(1, self._prime):
                out += (
                    -codomain_underlying_ring(
                        (binomial(self._prime, power) // self._prime)
                    )
                    * (a**power)
                    * (b ** (self._prime - power))
                )
            return out

        # We initialize the known frobenii from the base ring and set everything else to be zero.
        # This allows partial_frobenius and partial_delta to be callable but to give possibly
        # incorrect answers before complete initialization.

        codomain_frobenii = [
            underlying_ring_homomorphism(f)
            for f in self._domain.underlying_delta_ring()._frobenii
        ] + [codomain_underlying_ring.zero()] * (new_ngens - self._domain.ngens())

        codomain_deltas = [
            underlying_ring_homomorphism(f)
            for f in self._domain.underlying_delta_ring()._deltas
        ] + [codomain_underlying_ring.zero()] * (new_ngens - self._domain.ngens())

        def partial_frobenius(f):
            return f(codomain_frobenii)

        def partial_delta(f):
            return (partial_frobenius(f) - f**self._prime) // self._prime

        ### Initialize d powers.
        d_powers = []
        delta_d_powers = []
        frobenius_d_powers = []
        if max(num_f_array) > 1:
            d_powers.append(codomain_distinguished_element)
            print("a")
            delta_d_powers.append(partial_delta(d_powers[0]))
            frobenius_d_powers.append(partial_frobenius(codomain_distinguished_element))
        for u in range(max(num_f_array) - 1):
            d_powers.append(d_powers[u] ** self._prime)
            print("b")
            delta_d_powers.append(partial_delta(d_powers[u + 1]))
            frobenius_d_powers.append(frobenius_d_powers[u] ** self._prime)

        # d_times_lambda
        d_times_lambda = []
        delta_d_times_lambda = []
        inverses = []

        ### Initialize lambda_u.
        codomain_lambda_units = []
        if max(num_f_array) > 1:
            codomain_lambda_units.append(-(delta_d_powers[0].inverse()))

            d_times_lambda.append(codomain_lambda_units[0] * d_powers[1])
            print("c")
            delta_d_times_lambda.append(partial_delta(d_times_lambda[0]))
            inverses.append(
                (codomain_underlying_ring.one() - delta_d_times_lambda[0]).inverse()
            )

        # For f_0,...,f_r we need lambda_0,...,lambda_{r-1}.
        for u in range(max(num_f_array) - 2):
            codomain_lambda_units.append(
                (codomain_lambda_units[u] ** self._prime) * inverses[u]
            )
            d_times_lambda.append(codomain_lambda_units[u + 1] * d_powers[u + 2])
            print("d")
            delta_d_times_lambda.append(partial_delta(d_times_lambda[u + 1]))
            inverses.append(
                (codomain_underlying_ring.one() - delta_d_times_lambda[u + 1]).inverse()
            )

        ### Initialize relations.
        old_ngens = self._domain.ngens()
        offset = old_ngens

        # This does take a little work to coerce.
        # TODO: directly insert instead of calling underlying_ring_homomorphism, since
        # it is just insertion.
        codomain_relations = [underlying_ring_homomorphism(r) for r in self._relations]

        codomain_rs = []
        pth_powers = []

        for a in range(len(self._relations)):
            if num_f_array[a] == 1:
                codomain_rs.append(codomain_underlying_ring.zero())
                pth_powers.append(codomain_underlying_ring.zero())
                codomain_frobenii[offset] = codomain_underlying_ring.zero()
                codomain_deltas[offset] = codomain_underlying_ring.zero()
                offset += 1
            else:
                # The first R and powers.

                ### Temporary variables.
                print("e")
                codomain_rs.append(
                    partial_delta(codomain_relations[a]) / delta_d_powers[0]
                )

                pth_powers.append(
                    (-codomain_underlying_ring(self._prime) + d_times_lambda[0])
                    * codomain_underlying_ring.gens()[offset + 1]
                    + d_powers[1] * codomain_rs[offset - old_ngens]
                )
                codomain_deltas[offset] = (
                    codomain_underlying_ring.one()
                    + delta_d_powers[0] * codomain_lambda_units[0]
                ) * codomain_underlying_ring.gens()[offset] + delta_d_powers[
                    0
                ] * codomain_rs[
                    offset - old_ngens
                ]
                codomain_frobenii[offset] = frobenius_d_powers[0] * (
                    codomain_lambda_units[0]
                    * codomain_underlying_ring.gens()[offset + 1]
                    + codomain_rs[offset - old_ngens]
                )

                offset += 1

                for b in range(num_f_array[a] - 2):
                    print("f")
                    codomain_rs.append(
                        (
                            inverses[b]
                            * (partial_delta(codomain_rs[offset - old_ngens - 1])
                            + w1(
                                codomain_lambda_units[b]
                                * codomain_underlying_ring.gens()[offset],
                                codomain_rs[offset - old_ngens - 1],
                            ))
                        )
                    )
                    pth_powers.append(
                        (
                            (
                                -codomain_underlying_ring(self._prime)
                                + d_times_lambda[b + 1]
                            )
                            * codomain_underlying_ring.gens()[offset + 1]
                            + d_powers[b + 2] * codomain_rs[offset - old_ngens]
                        )
                    )
                    codomain_deltas[offset] = (
                        codomain_underlying_ring.one()
                        + delta_d_powers[b + 1] * codomain_lambda_units[b + 1]
                    ) * codomain_underlying_ring.gens()[offset + 1] + delta_d_powers[
                        b + 1
                    ] * codomain_rs[
                        offset - old_ngens - 1
                    ]
                    codomain_frobenii[offset] = frobenius_d_powers[b + 1] * (
                        codomain_lambda_units[b + 1]
                        * codomain_underlying_ring.gens()[offset + 1]
                        + codomain_rs[offset - old_ngens - 1]
                    )

                    offset += 1

                # Now append zero for the top R.
                codomain_rs.append(codomain_underlying_ring.zero())
                pth_powers.append(codomain_underlying_ring.zero())
                codomain_frobenii[offset] = codomain_underlying_ring.zero()
                codomain_deltas[offset] = codomain_underlying_ring.zero()
                offset += 1

        # Make new delta ring with new variables.
        codomain_delta_ring = DeltaPowerSeriesCapped(
            self._prime,
            codomain_underlying_ring,
            codomain_deltas,
            frobenii=codomain_frobenii,
        )

        # Make new delta ring homomorphism.
        delta_ring_morphism = DeltaPowerSeriesCappedHomomorphism(
            self._domain.underlying_delta_ring(),
            codomain_delta_ring,
            underlying_ring_homomorphism,
        )

        # Make new prism.
        self._codomain = Prism(codomain_delta_ring, codomain_distinguished_element)

        # Make new prism homomorphism.
        self._homomorphism = PrismHomomorphism(
            self._domain, self._codomain, delta_ring_morphism
        )
