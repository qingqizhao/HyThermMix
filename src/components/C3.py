class C3(object):
    """Multiparameter equation of state for propane"""
    name = "propane"
    CASNumber = "74-98-6"
    formula = "CH3CH2CH3"
    synonym = "R-290"
    _refPropName = "PROPANE"
    _coolPropName = "n-Propane"
    # rhoc = unidades.Density(220.4781)
    A = -2.38359E-01
    B = -9.72242E-01
    C = -4.71452E-02
    D = -3.18179E-03
    rhocL = 5.000043088
    rhoc = 220.4781
    # Tc = unidades.Temperature(369.89)
    Tc = 369.825
    # Pc = unidades.Pressure(4251.2, "kPa")
    Pc = 4251.2 * 1e3  # "kPa"
    M = 44.09562  # g/mol
    # Tt = unidades.Temperature(85.525)
    Tt = 85.525
    # Tb = unidades.Temperature(231.036)
    Tb = 231.036
    f_acent = 0.1521
    # momentoDipolar = unidades.DipoleMoment(0.084, "Debye")
    id = 4
    # _Tr = unidades.Temperature(354.964211)
    _Tr = 354.964211
    # _rhor = unidades.Density(221.906745)
    _rhor = 221.906745
    _w = 0.149041513

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-4.970583, 4.29352],
           "ao_exp": [3.043, 5.874, 9.337, 7.922],
           "titao": [393 / Tc, 1237 / Tc, 1984 / Tc, 4351 / Tc]}

    Fi2 = {"ao_log": [1, 3.02939],
           "pow": [0, 1],
           "ao_pow": [31.602908195, -84.463284382],
           "ao_exp": [], "titao": [],
           "ao_sinh": [6.60569, 19.1921], "sinh": [479.856 / Tc, 955.312 / Tc],
           "ao_cosh": [3.197, -8.37267], "cosh": [200.893 / Tc, 1027.29 / Tc]}

    CP1 = {"ao": 4.02939,
           "ao_sinh": [6.60569, 19.1921], "sinh": [479.856, 955.312],
           "ao_cosh": [3.197, -8.37267], "cosh": [200.893, 1027.29]}

    Fi3 = {"ao_log": [1, 3.02256195],
           "pow": [0, 1],
           "ao_pow": [10.14394256, -4.79513693],
           "ao_exp": [2.90591124, 4.68495401, 10.2971154, 8.08977905],
           "titao": [1.0515052038, 3.0961635368, 5.0845797877, 11.4329447982]}

    Fi4 = {"ao_log": [1, 3.021394],
           "pow": [0, 1],
           "ao_pow": [-4.992402, 4.291476],
           "ao_exp": [2.889980, 4.474243, 8.139803, 10.48251],
           "titao": [1.048309, 3.053170, 11.42280, 5.042815]}

    CP5 = {"ao": -5.4041204338,
           "an": [3.1252450099e6, -1.1415253638e5, 1.4971650720e3,
                  3.9215452897e-2, -2.1738913926e-5, 4.8274541303e-9],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [3.1907016349], "exp": [1500]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propane of Kunz and "
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

        "Tmin": 85.48, "Tmax": 500.0, "Pmax": 100000.0, "rhomax": 17.41,

        "nr1": [1.0403973107358, -2.8318404081403, 0.84393809606294,
                -0.076559591850023, 0.094697373057280, 0.24796475497006e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.2774376042287, -0.043846000648377, -0.2699106478435,
                -0.069313413089860, -0.029632145981653, 0.014040126751380],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1] * 6}

    eq = GERG