class nC10(object):
    """Multiparameter equation of state for n-decane"""
    name = "decane"
    CASNumber = "124-18-5"
    formula = "CH3-(CH2)8-CH3"
    synonym = ""
    _refPropName = "DECANE"
    _coolPropName = "n-Decane"
    A = -1.30465E-01
    B = -7.29921E-01
    C = -3.77324E-02
    D = -1.78622E-03
    rhocL = 1.64
    rhoc = 233.3419552
    Tc = 617.7
    Pc = 2103.0*1e3 # "kPa"
    M = 142.28168  # g/mol
    Tt = 243.5
    Tb = 447.27
    f_acent = 0.4884
    #momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 14

    Fi1 = {"ao_log": [1, 18.109],
           "pow": [0, 1],
           "ao_pow": [13.9361966549, -10.5265128286],
           "ao_exp": [25.685, 28.233, 12.417, 10.035],
           "titao": [1193/Tc, 2140/Tc, 4763/Tc, 10862/Tc]}

    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [15.870791919, -108.858547525],
           "ao_sinh": [21.0069, 58.3657], "sinh": [0.267034159, 2.833479035],
           "ao_cosh": [43.4931], "cosh": [1.353835195]}


    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propane of Kunz and "
                    "Wagner (2008).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 675.0, "Pmax": 800000.0, "rhomax": 5.41,

        "nr1": [0.10461e1, -0.24807e1, 0.74372, -0.52579, 0.15315, 0.32865e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [.84178, .55424e-1, -.73555, -.18507, -.20775e-1, .12335e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    #eq = lemmon, GERG