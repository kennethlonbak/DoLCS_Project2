import pylab as py

def SW_FLEX_RIGI(E_f,t_f,E_c,t_c,**kwargs):
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

def BEAM_SSE_UNI_PRESS_DEFORM(x, L, pres, **kwargs):
    '''
    calculates deflection for BEAM that has Simply Supported Edges subjected to UNIform PRESsure

    :param x: (float or numpy array) x-coordinate for the beam
    :param L: (float) Length of the beam
    :param pres: (float) Uniform pressure
    :param kwargs: (dict) Needs to containe D and S or information to compute it!
    :return: (float or numpy array, depending on x) Total deflection at x
    '''

    # Getting Flex-Rigi
    if "D" in kwargs:
        D = kwargs["D"]
    else:
        D = SW_FLEX_RIGI(**kwargs)

    # Getting Shear-Stiffness
    if "S" in kwargs:
        S = kwargs["S"]
    else:
        S = SW_SHEAR_STIFF(**kwargs)

    # Relative x distance
    xl = x/L

    # Bending deflection
    d_b = pres*L**4/(24*D)*(xl**4-2*xl**3+xl)

    # Shear deflection
    d_s = pres*L**2/(2*S)*(xl-xl**2)

    # Total deflection
    d = d_b + d_s
    return d


def BEAM_BEC_UNI_PRESS_DEFORM(x, L, pres, **kwargs):
    '''
    calculates deflection for BEAM that has Both Edges Clamped subjected to UNIform PRESsure

    :param x: (float or numpy array) x-coordinate for the beam
    :param L: (float) Length of the beam
    :param pres: (float) Uniform pressure
    :param kwargs: (dict) Needs to containe D and S or information to compute it!
    :return: (float or numpy array, depending on x) Total deflection at x
    '''

    # Getting Flex-Rigi
    if "D" in kwargs:
        D = kwargs["D"]
    else:
        D = SW_FLEX_RIGI(**kwargs)

    # Getting Shear-Stiffness
    if "S" in kwargs:
        S = kwargs["S"]
    else:
        S = SW_SHEAR_STIFF(**kwargs)

    # Relative x distance
    xl = x / L

    # Bending deflection
    d_b = pres * L ** 4 / (24 * D) * (xl ** 4 - 2 * xl ** 3 + xl**2)

    # Shear deflection
    d_s = pres * L ** 2 / (2 * S) * (xl - xl ** 2)

    # Total deflection
    d = d_b + d_s
    return d