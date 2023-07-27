# Imports
load('precision.py')

import sys
target_lmfdb=sys.argv[1]
load(target_lmfdb)

for field in range(10):
    # Create a directory a.b.c.d
    # sage.misc.misc.sage_makedirs('data/'+labels[field])
    
    # Open a file a.b.c.d/input.txt
    input_file = open('data/'+labels[field]+'_input.txt','w')
    
    # Populate a.b.c.d/input.txt with prime, weight bounds, n bounds, precision bound, and Eisenstein polynomial    
    input_file.write(labels[field]+'\n')
    
    # The prime p
    p=data[field][0]
    input_file.write(str(p)+'\n')
    
    # The possible weights
    # SET BY USER
    i_min=p-1
    i_max=p-1+10
    input_file.write(str(i_min)+'\n')
    input_file.write(str(i_max)+'\n')
    
    # The possible n
    # SET BY USER
    n_min=2
    n_max=9
    input_file.write(str(n_min)+'\n')
    input_file.write(str(n_max)+'\n')
    
    # The residual degree, currently 1
    f=1
    input_file.write(str(f)+'\n')
    
    # The maximum Fprecision
    Fprec_max = n_max*i_max
    input_file.write(str(Fprec_max)+'\n')

    # The target precision
    target_precision=nygaard_exponent(p,i_max,n_max)

    ####################
    # Precision losses #
    ####################

    ### From \delta
    precision_loss_delta=math.floor(math.log(Fprec_max-1,p))

    ### From passing from OK to OK/pi^n
    precision_loss_quotient=0
    for q in range(p,i_max):
        precision_loss_quotient+=n_max*valuation(p,math.factorial(q))

    ### From lifting nabla to Nygaard
    precision_loss_nygaard=n_max*math.floor(i_max*(i_max-1)/2)

    ### From computing can-phi on primitives
    precision_loss_primitives=target_precision
    
    ### From renormalizing the Eisenstein polynomial
    precision_loss_from_Eisenstein=1

    total_precision=(target_precision+precision_loss_delta
                     +precision_loss_quotient
                     +precision_loss_nygaard
                     +precision_loss_primitives
                     +precision_loss_from_Eisenstein)
    input_file.write(str(total_precision)+'\n')
    
    # Write an Eisenstein polynomial
    if len(data[field][1])==2:
        input_file.write('['+str(p)+',1]'+'\n')
    else:
        if f==1:
            W=Zp(p,total_precision,type='capped-abs',print_mode='terse',show_prec=False)
            E_coefficients=[W(c) for c in data[field][1]]
            constant_term=E_coefficients[0]
            E_normalized_coefficients=[c*(p/constant_term) for c in E_coefficients]
            input_file.write(str(E_normalized_coefficients)+'\n')
        else:
            raise NotImplementedError("Extensions with residual degree more than 1 are not implemented.")
    
    # Close the file
    input_file.close()

for field in range(10):
    # Load data from a.b.c.d_input.txt
    input_file = open('data/'+labels[field]+'_input.txt','r')
    lines = input_file.readlines()
    input_file.close()
    stripped_lines = [l.strip() for l in lines]
    p=ZZ(stripped_lines[1])
    i_min=ZZ(stripped_lines[2])
    i_max=ZZ(stripped_lines[3])
    n_min=ZZ(stripped_lines[4])
    n_max=ZZ(stripped_lines[5])
    f=ZZ(stripped_lines[6])
    Fprec_max=ZZ(stripped_lines[7])
    target_precision=ZZ(stripped_lines[8])
    E=eval(stripped_lines[9])
    for n in range(n_min,n_max+1):
        for i in range(i_min,i_max+1):
            # Create a file a.b.c.d/n/i_in.txt for each n,i
            input_file = open('data/'+labels[field]+'.'+str(n)+'.'+str(i)+'_in.txt','w')
            
            # Populate    
            input_file.write(labels[field]+'\n')

            # The prime p
            input_file.write(str(p)+'\n')

            # The weight
            input_file.write(str(i)+'\n')

            # n
            input_file.write(str(n)+'\n')

            # The residual degree, currently 1
            input_file.write(str(f)+'\n')

            # The maximum Fprecision
            Fprec = n*i
            input_file.write(str(Fprec)+'\n')

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

            total_precision=(target_precision+precision_loss_delta
                             +precision_loss_quotient
                             +precision_loss_nygaard
                             +precision_loss_primitives)
            input_file.write(str(total_precision)+'\n')

            # Write an Eisenstein polynomial
            W=Zp(p,total_precision,type='capped-abs',print_mode='terse',show_prec=False)
            new_E=[W(coeff) for coeff in E]
            input_file.write(str(new_E)+'\n')

            # Close the file
            input_file.close()

