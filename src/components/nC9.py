class nC9(object):
    """Multiparameter equation of state for n-nonane"""
    name = "nonane"
    CASNumber = "111-84-2"
    formula = "CH3-(CH2)7-CH3"
    synonym = ""
    _refPropName = "NONANE"
    _coolPropName = "n-Nonane"
    A = -4.38367E-01
    B = -9.11909E-01
    C = -6.50389E-02
    D = -3.26398E-03
    rhocL = 1.81
    rhoc = 232.141731
    Tc = 594.55
    Pc = 2281.0*1e3 # "kPa"
    M = 128.2551  # g/mol
    Tt = 219.7
    Tb = 423.91
    f_acent = 0.4433
    #momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 13

    Fi1 = {"ao_log": [1, 16.349],
           "pow": [0, 1],
           "ao_pow": [10.7927224829, -8.2418318753],
           "ao_exp": [24.926, 24.842, 11.188, 17.483],
           "titao": [1221/Tc, 2244/Tc, 5008/Tc, 11724/Tc]}

    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [16.313913248, -102.160247463],
           "ao_sinh": [18.0241, 53.3415], "sinh": [0.263819696, 2.848860483],
           "ao_cosh": [38.1235], "cosh": [1.370586158]}


    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for nonane of Kunz and "
                    "Wagner (2008).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi":  "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 800000.0, "rhomax": 6.06,

        "nr1": [0.11151e1, -0.27020e1, 0.83416, -0.38828, 0.13760, 0.28185e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.62037, 0.015847, -0.61726, -0.15043, -0.012982, 0.0044325],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    #eq = lemmon, GERG