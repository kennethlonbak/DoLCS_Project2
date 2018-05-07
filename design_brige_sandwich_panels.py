import pylab as py
import sandwich_beam as SWB
import ABD_matrix
from collections import OrderedDict as Dict
from scipy.optimize import minimize

''' -----------------------  Reading material files ---------------------------------------------------------------- '''
def matfile2dict(filename):
    '''
    Reading a material file and returing at dict with the properties
    :param filename:
    :return:
    '''
    dict_out = {}
    with open(filename,"r") as file:
        for line in file:
            if line[0] == "#":
                pass
            elif not(dict_out):
                name_vec = line.split()
                for name in name_vec:
                    dict_out[name] = {}
            else:
                line = line.split()
                mat_name = line[0]
                values = line[1:]
                for value, name in zip(values,name_vec):
                    dict_out[name][mat_name] = float(value)
    return dict_out

def matfile2SIdict(filename):
    dict_out = matfile2dict(filename)

    for name, field in dict_out.items():
        for mat_name in field:
            if "E" in mat_name or "G" in mat_name or "sig" in mat_name or "tau" in mat_name:
                dict_out[name][mat_name] *= 1e6
    return dict_out

def test_matfile2dict():
    print("------ Divinycell -------")
    divin = matfile2dict("Divinycell.dat")
    for name, value in divin.items():
        print("%s: %s"%(name,value))
    print("\n------ Material -------")
    mat = matfile2dict("material.dat")
    for name, value in mat.items():
        print("%s: %s" % (name, value))

    print("\n\n------ Divinycell SI -------")
    divin = matfile2SIdict("Divinycell.dat")
    for name, value in divin.items():
        print("%s: %s" % (name, value))
    print("\n------ Material SI -------")
    mat = matfile2SIdict("material.dat")
    for name, value in mat.items():
        print("%s: %s" % (name, value))

''' ---------------------------------- Solving ABD for matrix layup ------------------------------------------------ '''
def fiber2prop(fib_prop, angle, thickness):
    fib2ABD_inp = {"fiber_nr": len(angle)}
    z_start = 0.0
    fib2ABD_inp["thickness"] = sum(thickness)
    # Generating input for ABD
    for i in range(1,len(angle)+1):
        fib2ABD_inp[i] = fib_prop.copy()
        fib2ABD_inp[i]["angle"] = angle[i - 1] * py.pi / 180
        fib2ABD_inp[i]["thickness"] = thickness[i - 1]
        fib2ABD_inp[i]["z_start"] = z_start - fib2ABD_inp["thickness"]/2
        fib2ABD_inp[i]["z_end"] = z_start + thickness[i - 1]- fib2ABD_inp["thickness"]/2

        z_start += thickness[i - 1]

    lami_prop = ABD_matrix.fib2ABD(fib2ABD_inp)
    return lami_prop

def test_fiber2prop():
    fib_prop = matfile2SIdict("material.dat")["fiber"]
    angle = [0.0,90.0]
    thickness = [1e-3]*2
    lami_prop = fiber2prop(fib_prop, angle, thickness)
    for name, value in lami_prop.items():
        print("%s: %s" % (name, value))



if __name__ == '__main__':
    pass
    #test_matfile2dict()
    #test_fiber2prop()
    #test_SW_BEAM2response()
    #test_param2weight_and_param2con()
    #param2opt()
    #inspect_opt_results()

