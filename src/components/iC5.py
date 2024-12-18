class iC5(object):
    """Multiparameter equation of state for isopentane"""
    name = "isopentane"
    CASNumber = "78-78-4"
    formula = "(CH3)2-CH-CH2-CH3"
    synonym = "R-601a"
    _refPropName = "IPENTANE"
    _coolPropName = "Isopentane"
    A = -1.42415E-01
    B = -8.88522E-01
    C = -4.56753E-02
    D = -3.21314E-03
    rhocL = 3.271
    rhoc = 235.9986594
    Tc = 460.35
    Pc = 3378.0*1e3 # "kPa"
    M = 72.14878  # g/mol
    Tt = 112.65
    Tb = 300.98
    f_acent = 0.2274
    #momentoDipolar = unidades.DipoleMoment(0.11, "Debye")
    id = 7

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [2.5822330405, 1.1609103419],
           "ao_exp": [7.4056, 9.5772, 15.765, 12.119],
           "titao": [442/Tc, 1109/Tc, 2069/Tc, 4193/Tc]}

    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [15.449907693, -101.298172792],
           "ao_sinh": [11.7618, 33.1688],
           "sinh": [0.635392636, 4.169371131],
           "ao_cosh": [20.1101],
           "cosh": [1.977271641]}

    f = 72.151/8.3143
    CP3 = {"ao": 0.396504*f,
           "an": [0.260678e-2*f, 0.93677e-5*f, -0.158286e-7*f,  0.76525e-11*f],
           "pow": [1, 2, 3, 4]}

    f = 4.184/8.3159524
    CP4 = {"ao": 21.3861*f,
           "ao_sinh": [2.1524504e8*f], "sinh": [1.70158e3],
           "ao_cosh": [2.8330244e7*f], "cosh": [7.75899e2]}


    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isopentane of Kunz and "
                    "Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 1000000.0, "rhomax": 13.3,

        "nr1": [0.11017531966644e1, -0.30082368531980e1, 0.99411904271336,
                -0.14008636562629, 0.11193995351286, 0.29548042541230e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.36370108598133, -0.48236083488293e-1, -0.35100280270615,
                -0.10185043812047, -0.35242601785454e-1, 0.19756797599888e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}


    #eq = lemmon, GERG, polt, starling