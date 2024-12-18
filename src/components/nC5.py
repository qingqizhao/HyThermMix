class nC5(object):
    """Multiparameter equation of state for n-pentane"""
    name = "pentane"
    CASNumber = "109-66-0"
    formula = "CH3-(CH2)3-CH3"
    synonym = "R-601"
    _refPropName = "PENTANE"
    _coolPropName = "n-Pentane"
    A = -2.31712E-01
    B = -9.08403E-01
    C = -4.81325E-02
    D = -2.26804E-03
    rhocL = 3.215577588
    rhoc = 232
    Tc = 469.7
    Pc = 3370.0 * 1e3  # "kPa"
    M = 72.14878  # g/mol
    Tt = 143.47
    Tb = 309.21
    f_acent = 0.251
    # momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 8
    _Tr = 449.271155
    _rhor = 233.873368
    _w = 0.247058753

    Fi1 = {"ao_log": [1, 3.0],
           "ao_pow": [14.536611217, --89.919548319],
           "ao_sinh": [8.95043, 33.4032], "sinh": [0.380391739, 3.777411113],
           "ao_cosh": [21.836], "cosh": [1.789520971]}

    CP0 = {"ao": 4,
           "ao_sinh": [8.95043, 33.4032], "sinh": [178.67, 1774.25],
           "ao_cosh": [21.836], "cosh": [840.538]}

    CP1 = {"ao": 10.288132,
           "an": [-0.2695377e-1, 0.20951065e-3, -0.27910773e-6, 0.12266269e-9],
           "pow": [1, 2, 3, 4]}

    f = 4.184 / 8.3159524
    CP2 = {"ao": 22.5012 * f,
           "ao_sinh": [2.057417e8 * f], "sinh": [1.71958e3],
           "ao_cosh": [2.972927e7 * f], "cosh": [8.02069e2]}

    CP3 = {"ao": 4,
           "ao_exp": [9.751560716, 22.71445741, 11.65392685],
           "exp": [404.8796661, 1785.491483, 4504.430788]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for pentane of Kunz and "
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

        "Tmin": 143.47, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 10.57,

        "nr1": [0.10968643098001e1, -0.29988888298061e1, 0.99516886799212,
                -0.16170708558539, 0.11334460072775, 0.26760595150748e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.40979881986931, -0.40876423083075e-1, -0.38169482469447,
                -0.10931956843993, -0.32073223327990e-1, 0.16877016216975e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1] * 6}

    # eq = shortSpan, GERG, polt, starling, sun, ratanapisit