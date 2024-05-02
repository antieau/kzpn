from sage.all import Zp, PowerSeriesRing
from syntomic.sage_import import sage_import

sage_import(
    "syntomic/prismatic_envelope", fromlist=["PrismaticEnvelopeF", "PrismaticEnvelopeG"]
)
sage_import("syntomic/syntomic_complex", fromlist=["SyntomicComplex"])

# User-defined inputs.
p = 2
n = 2
i = 2

# The F-precision.
prec_F = max(i * n, 1)
total_precision = 50

W = Zp(p, total_precision, print_mode="terse", type="capped-abs", show_prec=False)
A = PowerSeriesRing(W, "z", prec_F)
z = A.gen()

# An Eisenstein polynomial normalized so that E(0)=p.
E = z**2 + p*z + p

envf = PrismaticEnvelopeF(p, n, E, prec_F, total_precision)
envg = PrismaticEnvelopeG(p, E, prec_F, total_precision)
syn = SyntomicComplex(prec_F, i, envf, envg)
syn._compute_matrices()
syn._compute_chain_complex()
syn._compute_homology()
print(syn.homology.orders)