import pylab as py
import sandwich_beam as SWB
import ABD_matrix

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


''' ---------------------------------- Calculating panel respone --------------------------------------------------- '''
def SW_BEAM2response(L, W, press, fib_ang_com, t_f, t_c, core_grade, fib_prop, deff_max, fib_ang_bet=90.0,**kwargs):
    SW_inp = {}
    SW_inp["L"] = {}

    # Assigning L
    SW_inp["L"]["L"] = L
    SW_inp["L"]["press"] = press

    # set t_f and E_f
    SW_inp["L"]["t_f"] = t_f
    fib_angles = [fib_ang_com,fib_ang_com+fib_ang_bet,fib_ang_com+fib_ang_bet,fib_ang_com]
    #fib_angles = [fib_ang_com,fib_ang_com+fib_ang_bet,fib_ang_com,fib_ang_com+fib_ang_bet]
    fib_thickness = [t_f/4]*4
    lami_prop = fiber2prop(fib_prop, fib_angles, fib_thickness)
    SW_inp["L"]["E_f"] = lami_prop["E_x"]

    # set t_c and E_c (E_c is calculated as the mean of the compressive and tensile modulie)
    SW_inp["L"]["t_c"] = t_c
    SW_inp["L"]["E_c"] = (core_grade["E1_c"]+core_grade["E1_t"])/2.0

    # Set G_c
    SW_inp["L"]["G_c"] = core_grade["G23"]

    # Set rho_c and rho_f
    SW_inp["L"]["rho_f"] = lami_prop[1]["rho"]
    SW_inp["L"]["rho_c"] = core_grade["rho"]

    # Assigning width input
    SW_inp["W"] = SW_inp["L"].copy()
    SW_inp["W"]["E_f"] = lami_prop["E_y"]
    SW_inp["W"]["L"] = W

    out = {}
    out["L"] = {}
    out["W"] = {}
    # Calculate max deflection in L and W
    out["L"]["deff"] = SWB.BEAM_BEC_UNI_PRESS_MAX_DEFLE(**SW_inp["L"])
    out["W"]["deff"] = SWB.BEAM_SSE_UNI_PRESS_MAX_DEFLE(**SW_inp["W"])

    # Calculate max stress in L and W for face and core
    out["L"]["sig_max_f"], out["L"]["sig_max_c"] = SWB.BEAM_BEC_UNI_PRESS_MAX_STESS_F_C(**SW_inp["L"])
    out["W"]["sig_max_f"], out["W"]["sig_max_c"] = SWB.BEAM_SSE_UNI_PRESS_MAX_STESS_F_C(**SW_inp["W"])

    # Calculte max shear in L and L for face and core
    out["L"]["tau_max_f"], out["L"]["tau_max_c"] = SWB.BEAM_UNI_PRESS_MAX_SHEAR_STESS_F_C(**SW_inp["L"])
    out["W"]["tau_max_f"], out["W"]["tau_max_c"] = SWB.BEAM_UNI_PRESS_MAX_SHEAR_STESS_F_C(**SW_inp["W"])

    # Calculate weight of panel
    out["m"] = PANEL_WEIGHT(W=W,**SW_inp["L"])

    out["con"] = {}
    # Compare max deflection with the allow
    out["con"]["deff"] = {}
    out["con"]["deff"]["L"] = out["L"]["deff"]/deff_max
    out["con"]["deff"]["W"] = out["W"]["deff"]/deff_max

    # Compare max stress with the allow
    out["con"]["stress"] = {}
    # Face
    out["con"]["stress"]["face"] = {}
    out["con"]["stress"]["face"]["L"] = {}
    out["con"]["stress"]["face"]["W"] = {}
    out["con"]["stress"]["face"]["L"]["com"] = out["L"]["sig_max_c"]/lami_prop[1]["sig1_max_com"]
    out["con"]["stress"]["face"]["W"]["com"] = out["W"]["sig_max_c"]/lami_prop[1]["sig2_max_com"]

    out["con"]["stress"]["face"]["L"]["ten"] = out["L"]["sig_max_c"]/lami_prop[1]["sig1_max_ten"]
    out["con"]["stress"]["face"]["W"]["ten"] = out["W"]["sig_max_c"]/lami_prop[1]["sig2_max_ten"]

    # Core
    out["con"]["stress"]["core"] = {}
    out["con"]["stress"]["core"]["L"] = {}
    out["con"]["stress"]["core"]["W"] = {}
    out["con"]["stress"]["core"]["L"]["com"] = out["L"]["sig_max_c"]/core_grade["sig1_max_com"]
    out["con"]["stress"]["core"]["W"]["com"] = out["W"]["sig_max_c"]/core_grade["sig2_max_com"]

    out["con"]["stress"]["core"]["L"]["ten"] = out["L"]["sig_max_c"]/core_grade["sig1_max_ten"]
    out["con"]["stress"]["core"]["W"]["ten"] = out["W"]["sig_max_c"]/core_grade["sig2_max_ten"]

    # Compare max shear stress with the allow
    out["con"]["shear"] = {}
    out["con"]["shear"]["face"] = {}
    out["con"]["shear"]["core"] = {}

    out["con"]["shear"]["face"] = {}
    out["con"]["shear"]["face"]["L"] = out["L"]["tau_max_f"]/lami_prop[1]["tau_max"]
    out["con"]["shear"]["face"]["W"] = out["W"]["tau_max_f"]/lami_prop[1]["tau_max"]

    out["con"]["shear"]["core"] = {}
    out["con"]["shear"]["core"]["L"] = out["L"]["tau_max_c"]/core_grade["tau_max"]
    out["con"]["shear"]["core"]["W"] = out["W"]["tau_max_c"]/core_grade["tau_max"]

    out["con"]["is_con_okay"] = check_con(out["con"])
    if out["con"]["is_con_okay"]:
        out["con"]["closest"] = out["con"]["is_con_okay"].closest
    else:
        out["con"]["con_broken"] = out["con"]["is_con_okay"].con_broken

    return out





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


def test_SW_BEAM2response():
    L = 4.2 # [m]
    W = 4.0 # [m]
    weight_pressure = 5e3 # [kg/m^2]
    g = 9.82 # [m/s^2]
    press = weight_pressure*g # [N/m^2]
    fib_prop = matfile2SIdict("material.dat")["fiber"]
    core_grade = matfile2SIdict("Divinycell.dat")["H35"]
    t_tot = 0.25
    t_f_norm = 0.15
    t_f = t_tot*t_f_norm
    t_c = t_tot - 2*t_f
    fib_ang_com = 0.0
    fib_ang_bet = 90.0
    deff_max = 40e-3

    response = SW_BEAM2response(L, W, press, fib_ang_com, t_f, t_c, core_grade, fib_prop, deff_max, fib_ang_bet)
    for name, value in response.items():
        print("%s: %s" % (name, value))

''' ---------------------------- Optimizing Panel ------------------------------------------------------------------ '''
def minWeight(t_tot,t_f_norm):
    pass


if __name__ == '__main__':
    #test_matfile2dict()
    #test_fiber2prop()
    test_SW_BEAM2response()


