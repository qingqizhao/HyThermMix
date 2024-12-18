class nC8(object):
    """Multiparameter equation of state for n-octane"""
    name = "octane"
    CASNumber = "111-65-9"
    formula = "CH3-(CH2)6-CH3"
    synonym = ""
    _refPropName = "OCTANE"
    _coolPropName = "n-Octane"
    A = -4.98731E-01
    B = -9.09925E-01
    C = -5.66863E-02
    D = -2.66169E-03
    rhocL = 2.056404127
    rhoc = 234.8999588
    Tc = 569.32
    Pc = 2497.0 * 1e3  # "kPa"
    M = 114.2285  # g/mol
    Tt = 216.37
    Tb = 398.77
    f_acent = 0.395
    # momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 12
    _Tr = 565.427917
    _rhor = 234.605116
    _w = 0.402698435

    CP1 = {"ao": 4,
           "ao_sinh": [15.6865, 48.1731], "sinh": [158.9220, 1693.07],
           "ao_cosh": [33.8029], "cosh": [815.064]}

    Fi1 = {"ao_log": [1, 3.0],
           "ao_pow": [15.864687161, -97.370667555], "pow": [0, 1],
           "ao_sinh": [15.6865, 48.1731], "sinh": [158.9220 / Tc, 1693.07 / Tc],
           "ao_cosh": [33.8029], "cosh": [815.064 / Tc]}

    CP3 = {"ao": 3.018753,
           "an": [0.07297005, -0.14171168e-4, -0.1225317e-7, 0.12912645e-11],
           "pow": [1, 2, 3, 4]}

    f = 8.3159524 / 4.184
    CP4 = {"ao": 34.0847 * f,
           "ao_sinh": [2.603664e8 * f], "sinh": [1.6115500e3],
           "ao_cosh": [4.1241363e7 * f], "cosh": [7.6884700e2]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for n-octane of Kunz and "
                    "Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 6.69,

        "nr1": [0.10722544875633e1, -0.24632951172003e1, 0.65386674054928,
                -0.36324974085628, 0.12713269626764, 0.30713572777930e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.5265685698754, 0.19362862857653e-1, -0.58939426849155,
                -0.14069963991934, -0.78966330500036e-2, 0.33036597968109e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1] * 6}

    # eq = shortSpan, GERG, polt, starling, sun