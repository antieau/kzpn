# Imports
load('prismatic_envelope.sage')
load('precision.py')

# User-defined input

# The prime p
p=2

# The motivic weight i
i=2

# The power of the uniformizer n
n=2

# The residual degree f
# The present code only supports totally ramified extensions of Qp,
# i.e., where f=1
f=1

# Some precision calculations

# The calculated F-precision needed to compute at this weight
Fprec=n*i

# The target precision
target_precision=nygaard_exponent(p,i,n)

####################
# Precision losses #
####################

### From \delta
precision_loss_delta=math.floor(math.log(Fprec-1,p))

### From passing from OK to OK/pi^n
precision_loss_quotient=0
for q in range(p,i):
    precision_loss_quotient+=n*valuation(p,math.factorial(q))
    
### From lifting nabla to Nygaard
precision_loss_nygaard=n*math.floor(i*(i-1)/2)

### From computing can-phi on primitives
precision_loss_primitives=target_precision

### From renormalizing the Eisenstein polynomial
precision_loss_from_Eisenstein=1

total_precision=(target_precision+precision_loss_delta
                 +precision_loss_quotient
                 +precision_loss_nygaard
                 +precision_loss_primitives
                 +precision_loss_from_Eisenstein)

# The coefficient ring W
if f==1:
    W=Zp(p,total_precision,type='capped-abs',print_mode='terse',show_prec=False)
else:
    W=Zq(p**f,total_precision,names='a',type='capped-abs',print_mode='series',show_prec=False)
    
# The Breuil-Kisin ring A
A.<z>=PowerSeriesRing(W,Fprec)

# The Eisenstein polynomial E
E=A(z+p)

# The normalized Eisenstein polynomial
# The normalization is to bring the Eisenstein polynomial into the form E(0)=p
E=eisenstein_normalization(p,E)

# The syntomic matrices from new
syn0,syn1,nablaN,nablaP=syntomic_matrices(p,i,n,E,total_precision,Fprec,debug=False)

# The K-groups (cohomology of the p-adic syntomic complex)
# New
coh_dict,final_precision=syntomic_cohomology(syn0,syn1,nablaN,nablaP)

print(f"For f={1} and the Eisenstein polynomial {E}, with i={i} and n={n} we have:")
print('Elementary divisors of K_{}(R;Z_p)'.format(2*i-2)+' are {}'.format(coh_dict['h2'][1]))
print('Elementary divisors of K_{}(R;Z_p)'.format(2*i-1)+' are {}'.format(coh_dict['h1'][1]))
