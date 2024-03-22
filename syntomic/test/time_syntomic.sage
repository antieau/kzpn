import time
from syntomic.sage_import import sage_import

sage_import('syntomic/prismatic_envelope', fromlist=['PrismaticEnvelopeF', 'PrismaticEnvelopeG'])
sage_import('syntomic/syntomic_complex', fromlist=['SyntomicComplex'])


def timed(func):
    def timed_func(*args, **kwargs):
        t=0
        try:
            start = time.time()
            func(*args, **kwargs)
            end = time.time()
            t=end-start
        except:
            t='ERROR'
        return t
    return timed_func

@timed
def compute_syntomic(p,i,n,prec_F,total_precision):
    W=Zp(p,total_precision,print_mode='digits',type='capped-abs',show_prec=True)
    A=PowerSeriesRing(W,'z',prec_F)
    z=A.gen()
    envf = PrismaticEnvelopeF(p,n,z+p,prec_F,total_precision)
    envg = PrismaticEnvelopeG(p,z+p,prec_F,total_precision)
    syn = SyntomicComplex(envf,envg)
    syn._compute_matrices(i,prec_F)
    syn._compute_chain_complex(prec_F)
    syn._compute_homology()

table = []
for total_precision in range(5,1000,5):
    row = [] 
    for prec_F in range(5,25,5):
        row.append(compute_syntomic(2,2,2,prec_F,total_precision))
    print(','.join([str(it) for it in row]))
    table.append(row)


