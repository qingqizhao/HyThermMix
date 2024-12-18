class CO(object):
    """Multiparameter equation of state for carbon monoxide"""
    name = "carbon monoxide"
    CASNumber = "630-08-0"
    formula = "CO"
    synonym = ""
    _refPropName = "CO"
    _coolPropName = "CarbonMonoxide"
    A = -8.82515E-03
    B = -6.53553E-01
    C = 7.10728E-02
    D = 1.06831E-03
    rhocL = 10.85
    rhoc = 303.909585
    Tc = 132.86
    Pc = 3494.0 * 1e3  # "kPa"
    M = 28.0101  # g/mol
    Tt = 68.16
    Tb = 81.64
    f_acent = 0.0497
    # momentoDipolar = unidades.DipoleMoment(0.1, "Debye")
    id = 48

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1, -1.5],
           "ao_pow": [-3.3728318564, 3.3683460039, -9.111274701235156e-5],
           "ao_exp": [1.0128],
           "titao": [3089 / Tc]}

    Fi2 = {"ao_log": [1, 2.50055],
           "pow": [0, 1],
           "ao_pow": [10.813340744, -19.834733959],
           "ao_exp": [], "titao": [],
           "ao_sinh": [1.02865], "sinh": [11.6698028],
           "ao_cosh": [0.00493], "cosh": [5.302762306]}

    CP3 = {"ao": 0.36028218e1,
           "an": [-0.20871594e5, 0.89208708e3, -0.14157993e2, -0.34021345e-3,
                  0.44616091e-6, -0.15154703e-9],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [0.90426143], "exp": [30000]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon monoxide of Kunz "
                    "and Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 500., "Pmax": 100000.0, "rhomax": 33.84,

        "nr1": [0.92310041400851, -0.248858452058e1, 0.58095213783396,
                0.028859164394654, 0.070256257276544, 0.21687043269488e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 0.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.13758331015182, -0.51501116343466e-1, -0.14865357483379,
                -0.03885710088681, -0.029100433948943, 0.14155684466279e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1] * 6}

    # eq = lemmon, mccarty, GERG