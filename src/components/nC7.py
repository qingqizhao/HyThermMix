class nC7(object):
    """Multiparameter equation of state for n-heptane"""
    name = "heptane"
    CASNumber = "142-82-5"
    formula = "CH3-(CH2)5-CH3"
    synonym = ""
    _refPropName = "HEPTANE"
    _coolPropName = "n-Heptane"
    A = -2.47853E-01
    B = -8.41157E-01
    C = -4.54250E-02
    D = -2.21548E-03
    rhocL = 2.315324434
    rhoc = 232.0001389
    Tc = 540.13
    Pc = 2736.0*1e3 # "kPa"
    M = 100.202  # g/mol
    Tt = 182.55
    Tb = 371.53
    f_acent = 0.349
    #momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 11
    _Tr = 525.389862
    _rhor = 235.977855
    _w = 0.350780196

    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [15.063786601, -97.345252349],
           "ao_sinh": [13.7266, 43.55610], "sinh": [169.789/Tc, 1760.46/Tc],
           "ao_cosh": [30.4707], "cosh": [836.195/Tc]}

    CP1 = {"ao": 4,
           "ao_sinh": [13.7266, 43.55610], "sinh": [169.789, 1760.46],
           "ao_cosh": [30.4707], "cosh": [836.195]}

    CP3 = {"ao": 1.157528,
           "an": [0.070489617, -2.3419686e-5, -1.4768221e-9, -2.0117611e-12],
           "pow": [1, 2, 3, 4]}

    CP4 = {"ao": 30.4029/8.3159524*4.184,
           "ao_sinh": [2.5273083e8/8.3159524*4.184], "sinh": [1.6693200e3],
           "ao_cosh": [3.9046536e7/8.3159524*4.184], "cosh": [7.86001e2]}

    CP5 = {"ao": 4,
           "ao_exp": [15.29994054, 31.86604737, 14.10640675],
           "exp": [401.5547607, 1813.365387, 5041.869289]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for n-heptane of Kunz and "
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

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 7.75,

        "nr1": [0.10543747645262e1, -0.26500681506144e1, 0.81730047827543,
                -0.30451391253428, 0.122538687108, 0.27266472743928e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.49865825681670, -0.71432815084176e-3, -0.54236895525450,
                -0.13801821610756, -0.61595287380011e-2, 0.48602510393022e-3],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}


    #eq = shortSpan, GERG, polt, starling,  # ratanapisit