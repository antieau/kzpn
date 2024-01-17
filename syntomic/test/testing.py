def test(varnames,testcases):
    return lambda func: test_args(func,varnames,testcases)

def str_args(args,maxlength):
    argstring=""
    argstring += ','.join(["{}={}".format(k,v) for k,v in args.items()])
    if len(argstring)>maxlength:
        argstring = argstring[0:maxlength-1] + '..'
    return argstring


def test_args(func,varnames,testcases):
    for arg_list in testcases:
        varnames_list = varnames.split(',')
        args = {varnames_list[i]: arg_list[i] for i in range(len(varnames_list))}
        try:
            func(**args)
        except Exception as e:
            print('{}({}): \033[91m{}\033[00m'.format(func.__name__, str_args(args,40),'FAIL'))
            raise e
            return func
        print('{}({}): \033[92m{}\033[00m'.format(func.__name__, str_args(args,40),'PASS'))
    return func
