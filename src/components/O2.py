class O2(object):
    """Multiparameter equation of state for oxygen"""
    name = "oxygen"
    CASNumber = "7782-44-7"
    formula = "O2"
    synonym = "R-732"
    _refPropName = "OXYGEN"
    _coolPropName = "Oxygen"
    # rhoc = unidades.Density(436.143644)
    A = -1.21504E-01
    B = -9.52554E-01
    C = 3.72548E-02
    D = 8.72080E-03
    rhocL = 13.63
    rhoc = 436.143644
    # Tc = unidades.Temperature(154.581)
    Tc = 154.595
    # Pc = unidades.Pressure(5043.0, "kPa")
    Pc = 5043.0 * 1e3  # "kPa"
    M = 31.9988  # g/mol
    # Tt = unidades.Temperature(54.361)
    Tt = 54.361
    # Tb = unidades.Temperature(90.1878)
    Tb = 90.1878
    f_acent = 0.0222
    # momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 47
    # _Tr = unidades.Temperature(150.875090)
    _Tr = 150.875090
    # _rhor = unidades.Density(439.519141)
    _rhor = 439.519141
    _w = 0.023479051

    # Use the Cp0 expresion from reference in paper
    # Wagner, W., Ewers, J., Schmidt, R.
    # An Equation for the Ideal-Gas Heat Capacity of Molecular Oxygen for
    # Temperatures from 30 K to 3000 K
    # Ber. Bunsenges. Phys. Chem. 86 (1982) 538-540
    # doi: 10.1002@bbpc.19820860613
    # Fifth term not implemented
    CP1 = {"ao": 3.50042,
           "an": [1.06778, 1.66961e-8],
           "pow": [-1.5, 2],
           "ao_exp": [1.01258], "exp": [2242.45]}

    CP2 = {"ao": 3.521876773671,
           "an": [-0.4981998537119e4, 0.2302477799952e3, -0.3455653235107e1,
                  -0.4354202160244e-4, 0.1346353450132e-7, .1620598259591e-10],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [1.031468515726], "exp": [2239.18105]}

    Fi1 = {"ao_log": [1, 2.50146],
           "pow": [0, 1],
           "ao_pow": [10.001843586, -14.996095135],
           "ao_sinh": [1.07558], "sinh": [14.461722565],
           "ao_cosh": [1.01334], "cosh": [7.223325463]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for oxygen of Kunz and "
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

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 82000.0, "rhomax": 43.348,

        "nr1": [0.88878286369701, -0.24879433312148e1, 0.59750190775886,
                0.96501817061881e-2, 0.71970428712770e-1, 0.22337443000195e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.18558686391474, -0.038129368035760, -0.15352245383006,
                -0.026726814910919, -0.025675298677127, 0.95714302123668e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.623, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1] * 20,

        "nr3": [],
        "nr4": []}

    eq = GERG