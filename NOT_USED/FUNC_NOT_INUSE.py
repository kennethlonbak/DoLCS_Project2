def SW_BEAM2response(L, W, press, fib_ang_com, t_f, t_c, core_grade, fib_prop, deff_max, fib_ang_bet=90.0, use_clamped =True,**kwargs):
    SW_inp = {}
    SW_inp["L"] = {}

    # Assigning L
    SW_inp["L"]["L"] = L
    SW_inp["L"]["press"] = press

    # set t_f and E_f
    if isinstance(t_f,list):
        fib_thickness = t_f+t_f[::-1]
        SW_inp["L"]["t_f"] = sum(fib_thickness)
        fib_angles = [fib_ang_com, fib_ang_com + fib_ang_bet]*int(len(t_f)/2)
        fib_angles += fib_angles[::-1]
    else:
        fib_thickness = [t_f / 4] * 4
        SW_inp["L"]["t_f"] = t_f
        fib_angles = [fib_ang_com, fib_ang_com + fib_ang_bet, fib_ang_com + fib_ang_bet, fib_ang_com]

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
    # Adding inputs
    out["inp"] = SW_inp.copy()
    out["fiber"] = lami_prop
    out["L"] = {}
    out["W"] = {}
    # Calculate max deflection in L and W
    out["L"]["deff"] = SWB.BEAM_BEC_UNI_PRESS_MAX_DEFLE(**SW_inp["L"])
    if use_clamped:
        out["W"]["deff"] = SWB.BEAM_BEC_UNI_PRESS_MAX_DEFLE(**SW_inp["W"])
    else:
        out["W"]["deff"] = SWB.BEAM_SSE_UNI_PRESS_MAX_DEFLE(**SW_inp["W"])

    # Calculate max stress in L and W for face and core
    out["L"]["sig_max_f"], out["L"]["sig_max_c"] = SWB.BEAM_BEC_UNI_PRESS_MAX_STESS_F_C(**SW_inp["L"])
    if use_clamped:
        out["W"]["sig_max_f"], out["W"]["sig_max_c"] = SWB.BEAM_BEC_UNI_PRESS_MAX_STESS_F_C(**SW_inp["W"])
    else:
        out["W"]["sig_max_f"], out["W"]["sig_max_c"] = SWB.BEAM_SSE_UNI_PRESS_MAX_STESS_F_C(**SW_inp["W"])

    # Calculte max shear in L and L for face and core
    out["L"]["tau_max_f"], out["L"]["tau_max_c"] = SWB.BEAM_UNI_PRESS_MAX_SHEAR_STESS_F_C(**SW_inp["L"])
    out["W"]["tau_max_f"], out["W"]["tau_max_c"] = SWB.BEAM_UNI_PRESS_MAX_SHEAR_STESS_F_C(**SW_inp["W"])

    # Calculate weight of panel
    out["m"] = PANEL_WEIGHT(W=W,**SW_inp["L"])

    out["con"] = con2dict(out["L"]["deff"], out["W"]["deff"],
                          out["L"]["sig_max_f"], out["L"]["sig_max_c"],
                          out["W"]["sig_max_f"], out["W"]["sig_max_c"],
                          out["L"]["tau_max_f"], out["L"]["tau_max_c"],
                          out["W"]["tau_max_f"], out["W"]["tau_max_c"],
                          deff_max, lami_prop[1], core_grade)

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
    con["deff"] = Dict()
    con["deff"]["L"] = L_deff / deff_max
    con["deff"]["W"] = W_deff / deff_max

    # Compare max stress with the allow
    con["stress"] = Dict()
    # Face
    con["stress"]["face"] = Dict()
    con["stress"]["face"]["L"] = Dict()
    con["stress"]["face"]["W"] = Dict()
    # Compressive
    con["stress"]["face"]["L"]["com"] = L_sig_max_f / lami_prop["sig1_max_com"]
    con["stress"]["face"]["W"]["com"] = W_sig_max_f / lami_prop["sig2_max_com"]
    # Tensile
    con["stress"]["face"]["L"]["ten"] = L_sig_max_f / lami_prop["sig1_max_ten"]
    con["stress"]["face"]["W"]["ten"] = W_sig_max_f / lami_prop["sig2_max_ten"]

    # Core
    con["stress"]["core"] = Dict()
    con["stress"]["core"]["L"] = Dict()
    con["stress"]["core"]["W"] = Dict()
    con["stress"]["core"]["L"]["com"] = L_sig_max_c / core_prop["sig1_max_com"]
    con["stress"]["core"]["W"]["com"] = W_sig_max_c / core_prop["sig2_max_com"]

    con["stress"]["core"]["L"]["ten"] = L_sig_max_c / core_prop["sig1_max_ten"]
    con["stress"]["core"]["W"]["ten"] = W_sig_max_c / core_prop["sig2_max_ten"]

    # Compare max shear stress with the allow
    con["shear"] = Dict()
    con["shear"]["face"] = Dict()
    con["shear"]["core"] = Dict()

    con["shear"]["face"] = Dict()
    con["shear"]["face"]["L"] = L_tau_max_f / lami_prop["tau_max"]
    con["shear"]["face"]["W"] = W_tau_max_f / lami_prop["tau_max"]

    con["shear"]["core"] = Dict()
    con["shear"]["core"]["L"] = L_tau_max_c / core_prop["tau_max"]
    con["shear"]["core"]["W"] = W_tau_max_c / core_prop["tau_max"]
    return con

def conDict2conArray(con):
    array = con_rec(con,[])
    return py.array(array)

def con_rec(con_in, array):
    for name, value in con_in.items():
        if isinstance(value,Dict):
            array = con_rec(con_in[name], array)
        elif isinstance(value,float):
            array.append(value)
    return array


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
    core_grade = matfile2SIdict("Divinycell.dat")["H60"]
    t_tot = 0.25
    t_f_norm = 0.4
    t_f = t_tot*t_f_norm
    t_c = t_tot - 2*t_f
    fib_ang_com = 0.0
    fib_ang_bet = 90.0
    deff_max = 40e-3

    response = SW_BEAM2response(L, W, press, fib_ang_com, t_f, t_c, core_grade, fib_prop, deff_max, fib_ang_bet)
    for name, value in response.items():
        print("%s: %s" % (name, value))

''' ---------------------------- Optimizing Panel ------------------------------------------------------------------ '''
def param2weight(t_tot, t_f_norm, t_f_ratio, fib_ang_com, fib_ang_bet, L, W, press, core_grade, fib_prop, deff_max, **kwargs):
    t_f = t_tot * t_f_norm
    t_c = t_tot - 2 * t_f
    rho_f = fib_prop["rho"]
    rho_c = core_grade["rho"]
    weight = PANEL_WEIGHT(L, W, t_c, t_f, rho_c, rho_f).weight
    return weight

def x2weight(x, L, W, press, core_grade, fib_prop, deff_max):
    weigh = param2weight(*x, L, W, press, core_grade, fib_prop, deff_max)
    #print(x,weigh)
    return weigh

def param2con(t_tot, t_f_norm, t_f_ratio, fib_ang_com, fib_ang_bet, L, W, press, core_grade, fib_prop, deff_max, return_out= False,**kwargs):
    t_f = t_tot * t_f_norm
    t_c = t_tot - 2 * t_f
    t_f = [t_f/2*t_f_ratio,t_f/2*(1-t_f_ratio)]

    out = SW_BEAM2response(L, W, press, fib_ang_com, t_f, t_c, core_grade, fib_prop, deff_max, fib_ang_bet)#,use_clamped=True)
    con_array = conDict2conArray(out["con"])
    if return_out:
        return out
    else:
        return con_array

def x2con(x, L, W, press, core_grade, fib_prop, deff_max,**kwargs):
    con1 = param2con(*x, L, W, press, core_grade, fib_prop, deff_max)
    con = 1-con1
    return con

def test_param2weight_and_param2con():
    L = 4.2 # [m]
    W = 4.0 # [m]
    weight_pressure = 5e3 # [kg/m^2]
    g = 9.82 # [m/s^2]
    press = weight_pressure*g # [N/m^2]
    fib_prop = matfile2SIdict("material.dat")["fiber"]
    core_grade = matfile2SIdict("Divinycell.dat")["H60"]
    t_tot = 0.25
    t_f_norm = 0.4
    fib_ang_com = 0.0
    fib_ang_bet = 90.0
    t_f_ratio = 0.5
    deff_max = 40e-3

    weight = param2weight(t_tot, t_f_norm, t_f_ratio, fib_ang_com, fib_ang_bet, L, W, press, core_grade, fib_prop, deff_max)
    con_array = param2con(t_tot, t_f_norm, t_f_ratio, fib_ang_com, fib_ang_bet, L, W, press, core_grade, fib_prop, deff_max)
    print("weight=%s"%weight)
    print("con_array=%s" % con_array)

def x2dict(x,**kwargs):
    out = {}
    name_list = "t_tot, t_f_norm, t_f_ratio, fib_ang_com, fib_ang_bet, name_grade".split(", ")
    for val, name in zip(x,name_list):
        out[name] = val
    return out

def dict2x(**kwargs):
    out = []
    name_list = "t_tot, t_f_norm, t_f_ratio, fib_ang_com, fib_ang_bet, name_grade".split(", ")
    for name in name_list:
        out.append(kwargs[name])
    return py.array(out)

def param2opt():
    L = 4.2  # [m]
    W = 4.0  # [m]
    weight_pressure = 5e3  # [kg/m^2]
    g = 9.82  # [m/s^2]
    press = weight_pressure * g  # [N/m^2]
    fib_prop = matfile2SIdict("material.dat")["fiber"]
    core_grade_dict = matfile2SIdict("Divinycell.dat")
    name_grade = "H250"
    core_grade = core_grade_dict[name_grade]
    return_out = True
    deff_max = 40e-3
    t_tot = 0.25
    t_f_norm = 0.4
    fib_ang_com = 0.0
    fib_ang_bet = 90.0
    t_f_ratio = 0.5

    x0 = [t_tot, t_f_norm, t_f_ratio, fib_ang_com, fib_ang_bet]
    #args = (L, W, press, core_grade, fib_prop, deff_max)
    #constraints = [{"type": "ineq", "fun": x2con, "args": args}]
    bounds = (
        (0.01,0.25),
        (0.001,0.1),
        (0.9,1.0),
        (0,45),
        (45,90)
    )
    name_grade = "H250"
    core_grade = core_grade_dict[name_grade]
    weight_opt = None
    close_bond = 1
    locals_opt = locals().copy()

    # Brut force search
    for name_grade in ["H45"]:#core_grade_dict.keys():
        core_grade = core_grade_dict[name_grade]
        for t_f_norm in [0.0273668341709]:#py.linspace(*bounds[1], 200):
            for t_f_ratio in [0.927638190955]:#py.linspace(*bounds[2], 200):
                for fib_ang_com in [0.0]:#py.linspace(*bounds[3], 10):
                    for fib_ang_bet in [90.0]:#py.linspace(*bounds[4], 10):
                        x = dict2x(**locals())
                        out = param2con(**locals())
                        if out["con"]["is_con_okay"]:
                            if not(weight_opt) or out["m"].weight < weight_opt or (out["m"].weight == weight_opt and out["con"]["closest"]["value"] < close_bond):
                                weight_opt = out["m"].weight
                                close_bond = out["con"]["closest"]["value"]
                                locals_opt = locals().copy()
                                print(x2dict(**locals_opt))
                                #print(out["m"].weight, out["con"]["closest"]["value"] , locals_opt["name_grade"],
                                #      locals_opt["t_tot"],locals_opt["t_f_norm"],locals_opt["t_f_ratio"],
                                #      locals_opt["fib_ang_com"],locals_opt["fib_ang_bet"])

    if True:
        x0 = [locals_opt["t_tot"],locals_opt["t_f_norm"],locals_opt["t_f_ratio"],locals_opt["fib_ang_com"],locals_opt["fib_ang_bet"]]
        args = (L, W, press, locals_opt["core_grade"], fib_prop, deff_max)
        constraints = [{"type": "ineq", "fun": x2con, "args": args}]

        opt_res = minimize(fun = x2weight, x0 = x0, args=args, method="SLSQP", bounds=bounds, constraints=constraints, tol=1e-20, options={"maxiter":1e5,"disp": True} )
        con = x2con(opt_res.x, *args)
        x = opt_res.x
        print(x2dict(**locals()))

def inspect_opt_results():
    L = 4.2  # [m]
    W = 4.0  # [m]
    weight_pressure = 5e3  # [kg/m^2]
    g = 9.82  # [m/s^2]
    press = weight_pressure * g  # [N/m^2]
    fib_prop = matfile2SIdict("material.dat")["fiber"]
    core_grade_dict = matfile2SIdict("Divinycell.dat")
    deff_max = 40e-3
    name_grade = 'H45'
    core_grade = core_grade_dict[name_grade]
    return_out = True

    opt1 = {'t_tot': '0.25', 't_f_norm': '0.0273668341709', 't_f_ratio': '0.927638190955', 'fib_ang_com': '0.0', 'fib_ang_bet': '90.0'}
    opt2 = {'t_tot': 0.25, 't_f_norm': 0.02697526258222031, 't_f_ratio': 0.94052940393950557, 'fib_ang_com': 6.4283154067665111e-12, 'fib_ang_bet': 90.0}
    opt2 = {'t_tot': 0.25, 't_f_norm': 0.02998526258222031, 't_f_ratio': 1, 'fib_ang_com': 0.0, 'fib_ang_bet': 90.0}
    out = param2con(**locals(),**opt2)
    return


# 680.982857143 0.957495161506 H60 0.25 0.0277551020408 0.0888888888889 0.0 45.0
# 675.763603604 0.997710784896 H60 0.25 0.0274174174174 0.19797979798 0.0 90.0
# 674.981788945 0.993234653043 H60 0.25 0.0273668341709 0.595979899497 0.0 90.0
# 667.292623116 0.997470675351 H60 0.25 0.0268693467337 0.9 0.0 90.0
# 627.340365829 0.987625221499 H45 0.25 0.0273668341709 0.927638190955 0.0 90.0