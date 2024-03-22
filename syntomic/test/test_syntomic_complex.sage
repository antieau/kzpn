from syntomic.test.testing import test

from syntomic.sage_import import sage_import

sage_import('syntomic/prismatic_envelope', fromlist=['PrismaticEnvelopeF', 'PrismaticEnvelopeG'])
sage_import('syntomic/syntomic_complex', fromlist=['SyntomicComplex'])

syntomic_testcases=[
        (2,0,2,{0: [], 1: [], 2: []}),
        (2,1,2,{0: [], 1: [2], 2: []}),
        (2,2,2,{0: [], 1: [8], 2: [2]}),
        (2,3,2,{0: [], 1: [8], 2: []}),
        (2,4,2,{0: [], 1: [2,8], 2: []}),
        (2,5,2,{0: [], 1: [2,2,8], 2: []}),
        (2,6,2,{0: [], 1: [2,32], 2: []}),
        (2,7,2,{0: [], 1: [2,4,16], 2: []}),
        (2,2,3,{0: [], 1: [4,8], 2: [2]}),
        (2,3,3,{0: [], 1: [2,64], 2: [2]}),
        (2,4,3,{0: [], 1: [16,16], 2: []}),
        (3,5,2,{0: [], 1: [3,81], 2: []}),
        (5,4,2,{0: [], 1: [5,125], 2: []}),
        (5,5,2,{0: [], 1: [5^6], 2: [5]}),
        ]

@test('p,i,n,expected_homology',syntomic_testcases)
def test_syntomic(p,i,n,expected_homology):
    prec_F = max(i*n,1)
    total_precision = 50
    W=Zp(p,total_precision,print_mode='digits',type='capped-abs',show_prec=True)
    A=PowerSeriesRing(W,'z',prec_F)
    z=A.gen()
    envf = PrismaticEnvelopeF(p,n,z+p,prec_F,total_precision)
    envg = PrismaticEnvelopeG(p,z+p,prec_F,total_precision)
    syn = SyntomicComplex(prec_F, i,envf,envg)
    syn._compute_matrices()
    syn._compute_chain_complex()
    syn._compute_homology()
    assert syn.homology.orders == expected_homology, "expected: {}, got: {}".format(expected_homology, syn.homology.orders)

