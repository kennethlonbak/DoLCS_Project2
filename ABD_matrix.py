import pylab as py

def fib2ABD(fiber_layup):
    '''
    :param fiber_layup: type(fiber_layup)=dict, Mandatory: numbers from 1 to n_layer where n_layer is the number of layers. n_layer
    :subparam fiber_layup[i_layer]: type(fiber_layup[i_layer])=dict, Mandatory fields: E1, E2, nu12, G12, angle, thickness, z_start, z_end. Optional fields: fiber_type,
    :return fiber_layup["ABD"] = ABD: ["thickness"] = sum(fiber_layup[i_layer]["thickness"])
    :subfield fiber_layup[i_layer]: ["S_l"] = S_l, ["Q_l"] = Q_l, ["Q_G"] = Q_G, ["A"] = A, ["B"] = B, ["D"] = D
    '''

    A = py.zeros((3,3))
    B = py.zeros((3,3))
    D = py.zeros((3,3))

    for i_fib in range(1, fiber_layup["fiber_nr"] + 1):
        # Setting up local compliance matrix (S) for each fiber
        fiber_layup[i_fib]["S_l"] = get_compliance_matrix(fiber_layup[i_fib])

        # Getting local stiffness matrix (Q) - (Q = inv(S))
        fiber_layup[i_fib]["Q_l"] = py.inv(fiber_layup[i_fib]["S_l"])

        # Make transformation matrix for each fiber (T_lg, local->global - inv(T_lg)=T_gl, Global->local)
        fiber_layup[i_fib]["T_L2G"] = get_transform_Local2Global(fiber_layup[i_fib]["angle"])
        fiber_layup[i_fib]["T_G2L"] = py.inv(fiber_layup[i_fib]["T_L2G"])

        # Make global stiffness matrix
        fiber_layup[i_fib]["Q_G"] = py.dot(fiber_layup[i_fib]["T_L2G"], py.dot(fiber_layup[i_fib]["Q_l"], fiber_layup[i_fib]["T_L2G"].T))

        # Make global A, B, D
        A += fiber_layup[i_fib]["Q_G"]*(fiber_layup[i_fib]["z_end"]-fiber_layup[i_fib]["z_start"])
        B += fiber_layup[i_fib]["Q_G"] *1.0/2.0* (fiber_layup[i_fib]["z_end"]**2 - fiber_layup[i_fib]["z_start"]**2)
        D += fiber_layup[i_fib]["Q_G"] *1.0/3.0* (fiber_layup[i_fib]["z_end"]**3 - fiber_layup[i_fib]["z_start"]**3)

    # Collect ABD matrix
    fiber_layup["A"] = A
    fiber_layup["B"] = B
    fiber_layup["D"] = D
    ABD = py.zeros((6,6))
    ABD[:3,:3] = A
    ABD[:3, 3:] = B
    ABD[3:, :3] = B
    ABD[3:, 3:] = D
    fiber_layup["ABD"] = ABD
    fiber_layup["abd"] = py.inv(fiber_layup["ABD"])
    fiber_layup["E_x"] = 1.0/(fiber_layup["abd"][0,0]*fiber_layup["thickness"])
    fiber_layup["E_y"] = 1.0/(fiber_layup["abd"][1,1]*fiber_layup["thickness"])
    fiber_layup["G_xy"] = 1.0/(fiber_layup["abd"][2,2]*fiber_layup["thickness"])
    return(fiber_layup)


# Compliance and stiffness matrices ---------------------------------------------------------------------------------- #
def get_compliance_matrix(fiber_layup):
    S_11 = 1.0/fiber_layup["E1"]
    S_22 = 1.0/fiber_layup["E2"]
    S_12 = S_21 = -fiber_layup["nu12"]/fiber_layup["E1"]
    S_33 = 1.0/fiber_layup["G12"]
    S_l = py.array([[S_11,S_12,0],[S_21,S_22,0],[0,0,S_33]])
    return(S_l)

def get_transform_Local2Global(angle):
    c = py.cos(angle)
    s = py.sin(angle)
    T_L2G = py.array([[c**2 ,s**2   ,-2*c*s],
                      [s**2 ,c**2   ,2*c*s],
                      [c*s  ,-c*s   ,c**2-s**2]])
    return(T_L2G)

def get_test_layup():
    fiber_layup = {}
    th = py.array([1.0,2.0,3.0,4.0,5.0,6.0])*1e-3
    angle = py.array([0.0,10.0,20.0,30.0,40.0,50.0])
    E1 = 40000e6
    E2 = 11500e6
    v12 = 0.3
    G12 = 4500e6
    thickness = py.sum(th)
    th_cur = -thickness/2.0
    fiber_layup["thickness"] = thickness
    fiber_layup["fiber_nr"] = len(th)
    for i_lay in range(1,len(th)+1):
        fiber_layup[i_lay] = {}
        fiber_layup[i_lay]["E1"] = E1
        fiber_layup[i_lay]["E2"] = E2
        fiber_layup[i_lay]["nu12"] = v12
        fiber_layup[i_lay]["G12"] = G12
        fiber_layup[i_lay]["angle"] = angle[i_lay-1]*py.pi/180
        fiber_layup[i_lay]["thickness"] = th[i_lay-1]
        fiber_layup[i_lay]["z_start"] = th_cur
        th_cur += th[i_lay - 1]
        fiber_layup[i_lay]["z_end"] = th_cur

    return(fiber_layup)

def get_test_layup2():
    fiber_layup = {}
    th = py.array([1.0]*6)*1e-3
    angle = py.array([0.0]*len(th))
    E1 = 40000e6
    E2 = 11500e6
    v12 = 0.3
    G12 = 4500e6
    thickness = py.sum(th)
    th_cur = -thickness/2.0
    fiber_layup["thickness"] = thickness
    fiber_layup["fiber_nr"] = len(th)
    for i_lay in range(1,len(th)+1):
        fiber_layup[i_lay] = {}
        fiber_layup[i_lay]["E1"] = E1
        fiber_layup[i_lay]["E2"] = E2
        fiber_layup[i_lay]["nu12"] = v12
        fiber_layup[i_lay]["G12"] = G12
        fiber_layup[i_lay]["angle"] = angle[i_lay-1]*py.pi/180
        fiber_layup[i_lay]["thickness"] = th[i_lay-1]
        fiber_layup[i_lay]["z_start"] = th_cur
        th_cur += th[i_lay - 1]
        fiber_layup[i_lay]["z_end"] = th_cur

    return(fiber_layup)

def get_baseline_ABD():
    ABD_vec = [552392310.15491,
179551899.37054,
347479895.21026,
149248741.97482,
98110485.74160,
199677460.77180,
-996682.06764,
252950.03565,
490781.99633,
34503.75497,
376844.94749,
252950.03565,
21299.34707,
5925.46855,
13117.08398,
4259.69535,
3447.19648,
6665.08293]
    i_0 = 0
    A = py.array([[ABD_vec[i_0+0],ABD_vec[i_0+1],ABD_vec[i_0+3]],[ABD_vec[i_0+1],ABD_vec[i_0+2],ABD_vec[i_0+4]],[ABD_vec[i_0+3],ABD_vec[i_0+4],ABD_vec[i_0+5]]])
    i_0 = 6
    B = py.array([[ABD_vec[i_0+0],ABD_vec[i_0+1],ABD_vec[i_0+3]],[ABD_vec[i_0+1],ABD_vec[i_0+2],ABD_vec[i_0+4]],[ABD_vec[i_0+3],ABD_vec[i_0+4],ABD_vec[i_0+5]]])
    i_0 = 12
    D = py.array([[ABD_vec[i_0+0],ABD_vec[i_0+1],ABD_vec[i_0+3]],[ABD_vec[i_0+1],ABD_vec[i_0+2],ABD_vec[i_0+4]],[ABD_vec[i_0+3],ABD_vec[i_0+4],ABD_vec[i_0+5]]])

    ABD = py.zeros((6, 6))
    ABD[:3, :3] = A
    ABD[:3, 3:] = B
    ABD[3:, :3] = B
    ABD[3:, 3:] = D
    return (ABD)

def get_baseline_ABD2():
    ABD_vec = [246374951.87989,
21249839.59964,
70832798.66547,
0.00000,
0.00000,
27000000.00000,
0.00000,
0.00000,
0.00000,
0.00000,
0.00000,
0.00000,
739.12486,
63.74952,
212.49840,
0.00000,
0.00000,
81.00000]
    i_0 = 0
    A = py.array(
        [[ABD_vec[i_0 + 0], ABD_vec[i_0 + 1], ABD_vec[i_0 + 3]], [ABD_vec[i_0 + 1], ABD_vec[i_0 + 2], ABD_vec[i_0 + 4]],
         [ABD_vec[i_0 + 3], ABD_vec[i_0 + 4], ABD_vec[i_0 + 5]]])
    i_0 = 6
    B = py.array(
        [[ABD_vec[i_0 + 0], ABD_vec[i_0 + 1], ABD_vec[i_0 + 3]], [ABD_vec[i_0 + 1], ABD_vec[i_0 + 2], ABD_vec[i_0 + 4]],
         [ABD_vec[i_0 + 3], ABD_vec[i_0 + 4], ABD_vec[i_0 + 5]]])
    i_0 = 12
    D = py.array(
        [[ABD_vec[i_0 + 0], ABD_vec[i_0 + 1], ABD_vec[i_0 + 3]], [ABD_vec[i_0 + 1], ABD_vec[i_0 + 2], ABD_vec[i_0 + 4]],
         [ABD_vec[i_0 + 3], ABD_vec[i_0 + 4], ABD_vec[i_0 + 5]]])

    ABD = py.zeros((6, 6))
    ABD[:3, :3] = A
    ABD[:3, 3:] = B
    ABD[3:, :3] = B
    ABD[3:, 3:] = D
    return (ABD)

def extract_ABC_vec(ABD):
    ABD_vec = []
    name_vec = []
    # A values
    for i in range(3):
        for j in range(i, 3):
            ABD_vec.append(ABD[i, j])
            name_vec.append(r"$A_{%d%d}$"%(i+1,j+1))
    # B values
    for i in range(3):
        for j in range(i, 3):
            ABD_vec.append(ABD[3+i, j])
            name_vec.append(r"$B_{%d%d}$"%(i+1,j+1))
    # C values
    for i in range(3):
        for j in range(i, 3):
            ABD_vec.append(ABD[3 + i,3+ j])
            name_vec.append(r"$D_{%d%d}$"%(i+1,j+1))
    return(ABD_vec,name_vec)


if __name__ == "__main__":
    test_function()


