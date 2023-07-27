"""Commandline versus of script."""

import sys
set_verbose(-1)

input_file = open(str(sys.argv[1]),'r')
lines = input_file.readlines()
ls = [l.strip() for l in lines]

label=ls[0]
print('label is '+label)

p=ZZ(ls[1])
print('prime p is '+str(p))

i=ZZ(ls[2])
print('weight i is '+str(i))

n=ZZ(ls[3])
print('n is '+str(n))

f=ZZ(ls[4])
print('residual degree f is '+str(f))

Fprec=ZZ(ls[5])
print('Fprec is '+str(Fprec))

total_precision=ZZ(ls[6])
print('total_precision is '+str(total_precision))

# The coefficient ring W
if f==1:
    W=Zp(p,total_precision,type='capped-abs',print_mode='terse',show_prec=False)
else:
    raise NotImplementedError("Extensions with residual degree more than 1 are not implemented.")

# The Breuil-Kisin ring A
A.<z>=PowerSeriesRing(W,Fprec)

# The Eisenstein polynomial E
print('Eisenstein coefficients are '+ls[7])
E=A(eval(ls[7]))
print('Eisenstein polynomial is '+str(E))

input_file.close()

load('prismatic_envelope.sage')
load('precision.py')



    
#%%capture
# Suppresses some Python warnings and SAGE variable injections

# The syntomic matrices from new
syn0,syn1,nablaN,nablaP=syntomic_matrices(p,i,n,E,total_precision,Fprec,debug=False)

# The K-groups (cohomology of the p-adic syntomic complex)
# New
coh_dict,final_precision=syntomic_cohomology(syn0,syn1,nablaN,nablaP)

print('Elementary divisors of K_{}(R;Z_p)'.format(2*i-2)+' are {}'.format(coh_dict['h2'][1]))
print('Elementary divisors of K_{}(R;Z_p)'.format(2*i-1)+' are {}'.format(coh_dict['h1'][1]))
