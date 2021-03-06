import pylab as py
import scipy.integrate as integrate

def UNIT_WIDTH_SW_FLEX_RIGI(E_f,t_f,E_c,t_c,**kwargs):
    '''
    SandWich FLEXural RIGIdity
    Ref: An Introduction to Sandwich Structures (Dan Zenkert), p. 3.2, eq. 3.4

    :return: D
    '''
    # Centrioed face distance/thickness
    t_cf = t_c+t_f

    # Face sheets
    D_f = E_f*t_f**3/12

    # Coupling ?
    D_0 = E_f*t_f*t_cf**2/2

    # Core
    D_c = E_c*t_c**3/12

    # Total Flexural Rigidity
    D = 2*D_f+D_0+D_c
    return D

def SW_FLEX_RIGI(w, E_f, t_f, E_c, t_c, **kwargs):
    '''
    SandWich FLEXural RIGIdity
    Ref: An Introduction to Sandwich Structures (Dan Zenkert), p. 3.2, eq. 3.4

    :return: D
    '''
    L_x = t_c + 2*t_f
    L_y = w
    return FLEX_RIGI(L_x, L_y, E_x_SW, E_f=E_f, E_c = E_c, t_c=t_c)


def SW_SHEAR_STIFF(G_c,t_f,t_c,**kwargs):
    '''
    SandWich SHEAR STIFFness
    Approximate equation for the shear stiffness, using t_f << t_c and E_c << E_f
    Ref: An Introduction to Sandwich Structures (Dan Zenkert), p. 4.3, eq. 4.6

    :param G_c: Shear Modulus of the core
    :param t_f: Thickness of face sheets (Laminar)
    :param t_c: Thickness of core
    :param kwargs:
    :return:
    '''
    # Centrioed face distance/thickness
    t_cf = t_c + t_f

    # Shear Stiffness
    S = G_c*t_cf**2/t_c

    return S

def FLEX_RIGI(L_x, L_y, E_xy,  **kwargs):
    int_fun = lambda x,y: E_xy(x,y,**kwargs)*(x**2)
    bond = [[-L_x / 2, L_x / 2],[-L_y / 2, L_y / 2]]
    int_out = integrate.nquad(int_fun, bond)
    D = int_out[0]
    #print("D_error: %s"%int_out[1])
    return D

def E_x_SW(x, y, t_c, E_c, E_f, **kwargs):
    if abs(x) > t_c/2:
        return E_f
    else:
        return E_c

''' ----------- Simply Supported ----------------------------------------------------------------------------------- '''
def BEAM_SSE_UNI_PRESS_MAX_DEFLE(L, press, E_f, t_f, E_c, t_c, G_c, **kwargs):
    '''
    calculates MAXimum DEFLEction for BEAM that has Simply Supported Edges subjected to UNIform PRESsure
    Ref: An Introduction to Sandwich Structures (Dan Zenkert), p. 4.26

    :param x: (float or numpy array) x-coordinate for the beam
    :param L: (float) Length of the beam
    :param press: (float) Uniform pressure
    :param kwargs: (dict) Needs to containe D and S or information to compute it!
    :return: (float or numpy array, depending on x) Total deflection at x
    '''

    # Getting Flex-Rigi
    D = SW_FLEX_RIGI(E_f, t_f, E_c, t_c)

    # Getting Shear-Stiffness
    S = SW_SHEAR_STIFF(G_c,t_f,t_c)

    # Relative x distance
    xl = 1/2

    # Bending deflection
    d_b = press * L ** 4 / (24 * D) * (xl ** 4 - 2 * xl ** 3 + xl)

    # Shear deflection
    d_s = press * L ** 2 / (2 * S) * (xl - xl ** 2)

    # Total deflection
    d = d_b + d_s
    return d

def BEAM_SSE_UNI_PRESS_BM_SF(x, L, press, **kwargs):
    '''
    calculates Bending Moment and Shear Force for a BEAM that has Simply Supported Edges subjected to UNIform PRESsure
    Ref: An Introduction to Sandwich Structures (Dan Zenkert), p. 4.26

    :param x: (float or numpy array) x-coordinate for the beam
    :param L: (float) Length of the beam
    :param press: (float) Uniform pressure
    :param kwargs: (dict) Needs to containe D and S or information to compute it!
    :return: (float or numpy array, depending on x) Total deflection at x
    '''

    # Relative x distance
    xl = x/L

    # Bending Moment
    M = press * L ** 2 * xl * (1 - xl)

    # Shear Force
    T = press * L * (1 / 2 - xl)

    return M, T

def BEAM_SSE_UNI_PRESS_MAX_STESS_F_C(press, L, t_c, t_f, E_f, E_c, **kwargs):
    '''
    calculates MAXimum STRESS in the Face and Core for a BEAM that has Simply Supported Edges subjected to UNIform PRESsure
    Ref: An Introduction to Sandwich Structures (Dan Zenkert), combining eq. 3.8 (p. 3.3) with Mx from p. 4.26

    :return: (float or numpy array, depending on x) Total deflection at x
    '''
    # Getting Flex-Rigi
    D = SW_FLEX_RIGI(E_f,t_f,E_c,t_c)

    sig_max_f = press * L ** 2 * (t_c / 2 + t_f) * E_f / (4 * D)

    sig_max_c = press * L ** 2 * (t_c / 2) * E_c / (4 * D)

    return sig_max_f, sig_max_c


''' ----------- Clamped Edges -------------------------------------------------------------------------------------- '''
def BEAM_BEC_UNI_PRESS_MAX_DEFLE(L, press, E_f, t_f, E_c, t_c, G_c, **kwargs):
    '''
    calculates MAXimum DEFLEction for BEAM that has Both Edges Clamped subjected to UNIform PRESsure
    Ref: An Introduction to Sandwich Structures (Dan Zenkert), p. 4.27

    :param x: (float or numpy array) x-coordinate for the beam
    :param L: (float) Length of the beam
    :param press: (float) Uniform pressure
    :param kwargs: (dict) Needs to containe D and S or information to compute it!
    :return: (float or numpy array, depending on x) Total deflection at x
    '''

    # Getting Flex-Rigi
    D = SW_FLEX_RIGI(E_f, t_f, E_c, t_c)

    # Getting Shear-Stiffness
    S = SW_SHEAR_STIFF(G_c, t_f, t_c)

    # Relative x distance
    xl = 1/2

    # Bending deflection
    d_b = press * L ** 4 / (24 * D) * (xl ** 4 - 2 * xl ** 3 + xl ** 2)

    # Shear deflection
    d_s = press * L ** 2 / (2 * S) * (xl - xl ** 2)

    # Total deflection
    d = d_b + d_s
    return d

def BEAM_BEC_UNI_PRESS_BM_SF(x, L, press, **kwargs):
    '''
    calculates Bending Moment and Shear Force for a BEAM that has Both Edges Clamped subjected to UNIform PRESsure
    Ref: An Introduction to Sandwich Structures (Dan Zenkert), p. 4.27

    :param x: (float or numpy array) x-coordinate for the beam
    :param L: (float) Length of the beam
    :param press: (float) Uniform pressure
    :param kwargs: (dict) Needs to containe D and S or information to compute it!
    :return: (float or numpy array, depending on x) Total deflection at x
    '''

    # Relative x distance
    xl = x/L

    # Bending Moment
    M = press * L ** 2 / 2 * (xl - xl ** 2 - 1 / 6)

    # Shear Force
    T = press * L * (1 / 2 - xl)

    return M, T

def BEAM_BEC_UNI_PRESS_MAX_STESS_F_C(press, L, t_c, t_f, E_f, E_c, **kwargs):
    '''
    calculates MAXimum STRESS in the Face and Core for a BEAM that has Both Edges Clamped subjected to UNIform PRESsure
    Ref: An Introduction to Sandwich Structures (Dan Zenkert), combining eq. 3.8 (p. 3.3) with Mx from p. 4.27 eval at
    x/L = 1/2

    :return:
    '''
    # Getting Flex-Rigi
    D = SW_FLEX_RIGI(E_f,t_f,E_c,t_c)

    sig_max_f = press * L ** 2 * (t_c / 2 + t_f) * E_f / (24 * D)

    sig_max_c = press * L ** 2 * (t_c / 2) * E_c / (24 * D)

    return sig_max_f, sig_max_c

def BEAM_UNI_PRESS_MAX_SHEAR_STESS_F_C(press, L, t_c, t_f, E_f, E_c, **kwargs):
    '''
    calculates MAXimum STRESS in the Face and Core for a BEAM subjected to UNIform PRESsure
    Ref: An Introduction to Sandwich Structures (Dan Zenkert), combining eq. 3.14 and 3.15 (p. 3.3) with Tx from p. 4.26
    and 4.27 eval at x/L = 0

    :return:
    '''
    # Getting Flex-Rigi
    D = SW_FLEX_RIGI(E_f,t_f,E_c,t_c)

    Tx_max = press*L/2

    d = t_c+t_f

    tau_max_f = Tx_max/D*(E_f*t_f*d/2)

    tau_max_c = Tx_max/D*(E_f*t_f*d/2+E_c*t_c**2/8)

    return tau_max_f, tau_max_c

if __name__ == '__main__':
    t_tot = 2.0
    t_f_norm = 0.2
    t_f = t_tot*t_f_norm
    t_c = t_tot-2*t_f
    E_f = 2.0
    E_c = 1.0
    w = 1.0

    D_general = SW_FLEX_RIGI(w, E_f, t_f, E_c, t_c)
    D_unit = UNIT_WIDTH_SW_FLEX_RIGI(E_f,t_f,E_c,t_c)

    vari = locals().copy()
    for name, value in vari.items():
        print("%s: %s" % (name, value))
