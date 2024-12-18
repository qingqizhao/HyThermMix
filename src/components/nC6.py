class nC6(object):
    """Multiparameter equation of state for n-hexane"""
    name = "hexane"
    CASNumber = "110-54-3"
    formula = "CH3-(CH2)4-CH3"
    synonym = ""
    _refPropName = "HEXANE"
    _coolPropName = "n-Hexane"
    A = -2.79791E-01
    B = -8.84493E-01
    C = -5.20562E-02
    D = -3.13857E-03
    rhocL = 2.705877875
    rhoc = 233.18
    Tc = 507.82
    Pc = 3034.0 * 1e3  # "kPa"
    M = 86.17536  # g/mol
    Tt = 177.83
    Tb = 341.86
    f_acent = 0.299
    # momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 10
    _Tr = 487.762087
    _rhor = 235.700888
    _w = 0.298052404

    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [14.345969349, -96.165722367],
           "ao_exp": [], "titao": [],
           "ao_sinh": [11.6977, 38.6164], "sinh": [182.326 / Tc, 1826.59 / Tc],
           "ao_cosh": [26.8142], "cosh": [859.207 / Tc]}

    CP1 = {"ao": 4,
           "ao_sinh": [11.6977, 38.6164], "sinh": [182.326, 1826.59],
           "ao_cosh": [26.8142], "cosh": [859.207]}

    CP3 = {"ao": 2.5200507,
           "an": [0.05280653, -5.7861557e-6, -1.0899040e-8, -1.8988742e-13],
           "pow": [1, 2, 3, 4]}

    CP4 = {"ao": 26.6225 / 8.3159524 * 4.184,
           "ao_sinh": [2.3738446e8 / 8.3159524 * 4.184], "sinh": [1.71849e3],
           "ao_cosh": [3.5806766e7 / 8.3159524 * 4.184], "cosh": [8.02069e2]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for n-hexane of Kunz and "
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

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 8.85,

        "nr1": [0.10553238013661e1, -0.26120615890629e1, 0.76613882967260,
                -0.29770320622459, 0.11879907733358, 0.27922861062617e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.46347589844105, 0.011433196980297, -0.48256968738131,
                -0.093750558924659, -0.0067273247155994, -0.0051141583585428],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1] * 6}

    # eq = shortSpan, GERG, polt, starling, sun