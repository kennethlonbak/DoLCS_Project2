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

''' --------------------------------- Beam response ---------------------------------------------------------------- '''
def SW_BEAM2response(L, W, press, fib_ang_col, t_f, t_c, core_prop, fib_prop, deff_max, fib_ang_dif=90.0, use_clamped =True, **kwargs):
    out = {}

    # Assigning L
    out["L"] = L
    out["press"] = press

    # set t_f and E_f
    if isinstance(t_f,list):
        fib_thickness = t_f+t_f[::-1]
        out["t_f"] = sum(fib_thickness)
        fib_angles = [fib_ang_col, fib_ang_col + fib_ang_dif] * int(len(t_f) / 2)
        fib_angles += fib_angles[::-1]
    else:
        fib_thickness = [t_f / 4] * 4
        out["t_f"] = t_f
        fib_angles = [fib_ang_col, fib_ang_col + fib_ang_dif, fib_ang_col + fib_ang_dif, fib_ang_col]

    lami_prop = fiber2prop(fib_prop, fib_angles, fib_thickness)
    out["E_f"] = lami_prop["E_x"]

    # set t_c and E_c (E_c is calculated as the mean of the compressive and tensile modulie)
    out["t_c"] = t_c
    out["E_c"] = (core_prop["E1_c"] + core_prop["E1_t"]) / 2.0

    # Set G_c
    out["G_c"] = core_prop["G23"]

    # Set rho_c and rho_f
    out["rho_f"] = lami_prop[1]["rho"]
    out["rho_c"] = core_prop["rho"]

    # Adding inputs
    out["fiber"] = lami_prop
    out = {}
    # Calculate max deflection in L and W
    out["deff"] = SWB.BEAM_BEC_UNI_PRESS_MAX_DEFLE(**out)

    # Calculate max stress in L and W for face and core
    out["sig_max_f"], out["sig_max_c"] = SWB.BEAM_BEC_UNI_PRESS_MAX_STESS_F_C(**out)

    # Calculte max shear in L and L for face and core
    out["tau_max_f"], out["tau_max_c"] = SWB.BEAM_UNI_PRESS_MAX_SHEAR_STESS_F_C(**out)

    # Calculate weight of panel
    out["m"] = PANEL_WEIGHT(W=W,**out)

    out["con"] = con2dict(out["deff"], out["W"]["deff"],
                          out["sig_max_f"], out["sig_max_c"],
                          out["W"]["sig_max_f"], out["W"]["sig_max_c"],
                          out["tau_max_f"], out["tau_max_c"],
                          out["W"]["tau_max_f"], out["W"]["tau_max_c"],
                          deff_max, lami_prop[1], core_prop)

    out["con"]["is_con_okay"] = check_con(out["con"])
    if out["con"]["is_con_okay"]:
        out["con"]["closest"] = out["con"]["is_con_okay"].closest
    else:
        out["con"]["con_broken"] = out["con"]["is_con_okay"].con_broken

    return out

def con2dict(L_deff, W_deff,
             L_sig_max_f, L_sig_max_c,
             W_sig_max_f, W_sig_max_c,
             L_tau_max_f, L_tau_max_c,
             W_tau_max_f, W_tau_max_c,
             deff_max, lami_prop, core_prop):
    con = Dict()
    # Compare max deflection with the allow
    con["deff"] = L_deff / deff_max

    # Compare max stress with the allow
    con["stress"] = Dict()
    # Face
    con["stress"]["face"] = Dict()
    # Compressive
    con["stress"]["face"]["com"] = L_sig_max_f / lami_prop["sig1_max_com"]
    # Tensile
    con["stress"]["face"]["ten"] = L_sig_max_f / lami_prop["sig1_max_ten"]

    # Core
    con["stress"]["core"] = Dict()
    con["stress"]["core"]["com"] = L_sig_max_c / core_prop["sig1_max_com"]

    con["stress"]["core"]["ten"] = L_sig_max_c / core_prop["sig1_max_ten"]

    # Compare max shear stress with the allow
    con["shear"] = Dict()

    con["shear"]["face"] = L_tau_max_f / lami_prop["tau_max"]

    con["shear"]["core"] = L_tau_max_c / core_prop["tau_max"]
    return con

class PANEL_WEIGHT:
    def __init__(self,L,W,t_c,t_f,rho_c,rho_f,**kwargs):
        self.inps = [L,W,t_c,t_f,rho_c,rho_f]
        # Calculate core volume
        self.core_vol = L*W*t_c

        # Calculate face volume
        self.face_vol = 2*L*W*t_f

        # Calculate core weight
        self.core_weight = rho_c*self.core_vol

        # Calculate face weigth
        self.face_weight = rho_f*self.face_vol

        # Add up total weight
        self.weight = self.core_weight + self.face_weight

    def __repr__(self):
        return repr(self.weight)

class check_con:
    def __init__(self, con):
        is_okay, con_broken, closest = check_con_rec(con,[],[],[],{"value":0.0})
        self.is_okay = all(is_okay)
        if self.is_okay:
            self.closest = closest
        else:
            self.con_broken = con_broken

    def __repr__(self):
        return repr(self.is_okay)

    def __bool__(self):
        return self.is_okay

def check_con_rec(con, parrent, is_okay, con_broken, closest):
    for name, value in con.items():
        if isinstance(value,dict):
            is_okay, con_broken, closest = check_con_rec(value, parrent + [name], is_okay, con_broken, closest)
        elif value > 1:
            is_okay.append(False)
            con_broken.append({})
            con_broken[-1]["name"] = parrent + [name]
            con_broken[-1]["value"] = value
        elif value > closest["value"]:
            closest["name"] = parrent + [name]
            closest["value"] = value
            is_okay.append(True)
        else:
            is_okay.append(True)
    return is_okay, con_broken, closest

''' --------------------------------- Calculate current brige design ----------------------------------------------- '''
def woodprop2brigeprop():
    pass
    wood_prop = matfile2SIdict("material.dat")["wood"]
    rho = wood_prop["rho"]

    # Setting length properties
    Planks = {}
    Planks["L"] = 4.0
    Planks["W"] = 4.2
    Planks["T"] = 0.025

    Long = {}
    Long["L"] = 4.0
    Long["T"] = 0.225
    Long["W"] = 0.175
    Long["#"] = 7

    Trans = {}
    Trans["L"] = 4.9
    Trans["T"] = 0.45
    Trans["W"] = 0.2

    # Calculating Volume
    V = {}
    V["Planks"] = Planks["L"]*Planks["W"]*Planks["T"]
    V["Long"] = Long["L"]*Long["W"]*Long["T"]*Long["#"]
    V["Panel"] = V["Planks"]+V["Long"]
    V["Trans"] = Trans["L"]*Trans["W"]*Trans["T"]
    V["tot"] = V["Planks"]+V["Long"]+V["Trans"]

    print("----- Volume -----")
    for name, value in V.items():
        print("%s: %s"%(name, value))
    print()

    # Calculating Mass
    W = {}
    W["Planks"] = V["Planks"]*rho
    W["Long"] = V["Long"]*rho
    W["Panel"] = W["Planks"] + W["Long"]
    W["Trans"] = V["Trans"]*rho
    W["tot"] = W["Planks"]+W["Long"]+W["Trans"]

    print("----- Weight -----")
    for name, value in W.items():
        print("%s: %s" % (name, value))
    print()

    # Calculating deck height
    T_panel = Planks["T"]+Long["T"]
    print("Panel thickness: %s+%s = %s"%(Planks["T"],Long["T"],T_panel))
    T_tot = Planks["T"]+Long["T"]+Trans["T"]
    print("Bridge total deck thickness: %s+%s+%s = %s"%(Planks["T"],Long["T"],Trans["T"],T_tot))





''' --------------------------------- Calculate design 1 ----------------------------------------------------------- '''


''' --------------------------------- Calculate design 2 ----------------------------------------------------------- '''


''' --------------------------------- Calculate design 3 ----------------------------------------------------------- '''


''' ---------------------------------- Transvers beam -------------------------------------------------------------- '''


if __name__ == '__main__':
    pass
    #test_matfile2dict()
    #test_fiber2prop()
    woodprop2brigeprop()




