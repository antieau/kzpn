def name_vector(gen_names, vector):
    if(len(vector) != len(gen_names)):
        raise ValueError("Dimension mismatch, names: {}, vector: {}".format(len(gen_names), len(vector)))
    result=""
    first_term=True
    for i in range(len(vector)):
        if vector[i]==0:
            continue
        if not first_term:
            result += '+'
        first_term=False
        if vector[i]!=1: 
            result += str(vector[i])+'*'
        result += gen_names[i]
    if first_term:
        result = '0'
    return result

