class GERG(object):
    """Multiparameter equation of state GERG 2008
    ref http://dx.doi.org/10.1021/je300655b"""
    import os
    import pickle
    from .components.CH4 import CH4
    from .components.N2 import N2
    from .components.CO2 import CO2
    from .components.C2 import C2
    from .components.C3 import C3
    from .components.nC4 import nC4
    from .components.iC4 import iC4
    from .components.nC5 import nC5
    from .components.iC5 import iC5
    from .components.nC6 import nC6
    from .components.nC7 import nC7
    from .components.nC8 import nC8
    from .components.H2 import H2
    from .components.O2 import O2
    from .components.CO import CO
    from .components.H2O import H2O
    from .components.He import He
    from .components.Ar import Ar
    from .components.H2S import H2S
    from .components.nC9 import nC9
    from .components.nC10 import nC10

    Tref = 298.15
    Pref = 101325.

    kwargs = {"componente": [],
              "fraccion": [],
              "T": 0.0,
              "rho": 0.0,
              "P": 0.0,
              "v": 0.0,
              "h": None,
              "s": None,
              "u": None,
              "x": None,
              "mezcla": None}

    componentes = [CH4, N2, CO2, C2, C3, nC4,
                   iC4, nC5, iC5, nC6, nC7, nC8,
                   H2, O2, CO, H2O, He, Ar,
                   H2S, nC9, nC10]

    import os
    import pickle

    this_dir = os.path.dirname(__file__)

    # Paths to data files
    f_Fij_path = os.path.join(this_dir, '..', 'data', 'mEoS_Fij.pkl')
    f_Prop_c_path = os.path.join(this_dir, '..', 'data', 'Tc.pickle')

    # Load Fij data
    with open(f_Fij_path, 'rb') as f_Fij:
        Fij = pickle.load(f_Fij)

    # Load Prop_c data
    with open(f_Prop_c_path, 'rb') as f_Prop_c:
        Prop_c = pickle.load(f_Prop_c)

    fir_ij = {
        "0-1": {
            "nr1": [-0.98038985517335e-2, 0.42487270143005e-3],
            "d1": [1, 4],
            "t1": [0.000, 1.850],

            "nr2": [-.34800214576142e-1, -.13333813013896, -.11993694974627e-1,
                    0.69243379775168e-1, -0.31022508148249, 0.24495491753226,
                    0.22369816716981],
            "d2": [1, 2, 2, 2, 2, 2, 3],
            "t2": [7.850, 5.400, 0.000, 0.750, 2.800, 4.450, 4.250],
            "n2": [1, 1, 0.25, 0, 0, 0, 0],
            "e2": [0.5] * 7,
            "b2": [1, 1, 2.5, 3, 3, 3, 3],
            "g2": [0.5] * 7},

        "0-2": {
            "nr1": [-.10859387354942, .80228576727389e-1, -.93303985115717e-2],
            "d1": [1, 2, 3],
            "t1": [2.6, 1.95, 0],

            "nr2": [0.40989274005848e-1, -0.24338019772494, 0.23855347281124],
            "d2": [1, 2, 3],
            "t2": [3.95, 7.95, 8],
            "n2": [1, 0.5, 0],
            "e2": [0.5] * 3,
            "b2": [1, 2, 3],
            "g2": [0.5] * 3},

        "0-3": {
            "nr1": [-0.80926050298746e-3, -0.75381925080059e-3],
            "d1": [3, 4],
            "t1": [0.65, 1.55],

            "nr2": [-0.41618768891219e-1, -0.23452173681569, 0.14003840584586,
                    .63281744807738e-1, -.34660425848809e-1, -.23918747334251,
                    0.19855255066891e-2, 0.61777746171555e1,
                    -0.69575358271105e1, 0.10630185306388e1],
            "d2": [1, 2, 2, 2, 2, 2, 2, 3, 3, 3],
            "t2": [3.1, 5.9, 7.05, 3.35, 1.2, 5.8, 2.7, 0.45, 0.55, 1.95],
            "n2": [1, 1, 1, 0.875, 0.75, 0.5, 0, 0, 0, 0],
            "e2": [0.5] * 10,
            "b2": [1, 1, 1, 1.25, 1.5, 2, 3, 3, 3, 3],
            "g2": [0.5] * 10},

        "0-4": {
            "nr1": [.13746429958576e-1, -.74425012129552e-2,
                    -.45516600213685e-2, -0.54546603350237e-2,
                    0.23682016824471e-2],
            "d1": [3, 3, 4, 4, 4],
            "t1": [1.85, 3.95, 0, 1.85, 3.85],

            "nr2": [0.18007763721438, -0.44773942932486, 0.19327374888200e-1,
                    -0.30632197804624],
            "d2": [1, 1, 1, 2],
            "t2": [5.25, 3.85, 0.2, 6.5],
            "n2": [0.25, 0.25, 0, 0],
            "e2": [0.5, 0.5, 0.5, 0.5],
            "b2": [0.75, 1, 2, 3],
            "g2": [0.5, 0.5, 0.5, 0.5]},

        "1-2": {
            "nr1": [0.28661625028399, -0.10919833861247],
            "d1": [2, 3],
            "t1": [1.85, 1.4],

            "nr2": [-0.11374032082270e1, 0.76580544237358, 0.42638000926819e-2,
                    0.17673538204534],
            "d2": [1, 1, 1, 2],
            "t2": [3.2, 2.5, 8, 3.75],
            "n2": [0.25, 0.25, 0, 0],
            "e2": [0.5, 0.5, 0.5, 0.5],
            "b2": [0.75, 1., 2., 3.],
            "g2": [0.5, 0.5, 0.5, 0.5]},

        "1-3": {
            "nr1": [-0.47376518126608, 0.48961193461001, -0.57011062090535e-2],
            "d1": [2, 2, 3],
            "t1": [0, 0.05, 0],

            "nr2": [-0.19966820041320, -0.69411103101723, 0.69226192739021],
            "d2": [1, 2, 2],
            "t2": [3.65, 4.9, 4.45],
            "n2": [1, 1, 0.875],
            "e2": [0.5, 0.5, 0.5],
            "b2": [1, 1, 1.25],
            "g2": [0.5, 0.5, 0.5]},

        "0-12": {
            "nr1": [-.25157134971934, -.62203841111983e-2, .88850315184396e-1,
                    -0.35592212573239e-1],
            "d1": [1, 3, 3, 4],
            "t1": [2, -1, 1.75, 1.4],
            "nr2": []},

        "0-5": {
            "nr1": [.25574776844118e1, -.79846357136353e1, 0.47859131465806e1,
                    -0.73265392369587, 0.13805471345312e1, 0.28349603476365,
                    -0.49087385940425, -0.10291888921447, 0.11836314681968,
                    0.55527385721943e-4],
            "d1": [1, 1, 1, 2, 2, 3, 3, 4, 4, 4],
            "t1": [1.0, 1.550, 1.700, 0.250, 1.350, 0.0, 1.250, 0.0, 0.7, 5.4],
            "nr2": []},

        "0-6": {
            "nr1": [.25574776844118e1, -.79846357136353e1, 0.47859131465806e1,
                    -0.73265392369587, 0.13805471345312e1, 0.28349603476365,
                    -0.49087385940425, -0.10291888921447, 0.11836314681968,
                    0.55527385721943e-4],
            "d1": [1, 1, 1, 2, 2, 3, 3, 4, 4, 4],
            "t1": [1.0, 1.55, 1.70, 0.25, 1.35, 0.0, 1.25, 0.0, 0.70, 5.40],
            "nr2": []},

        "3-4": {
            "nr1": [0.25574776844118e1, -.79846357136353e1, 0.47859131465806e1,
                    -0.73265392369587, 0.13805471345312e1, 0.28349603476365,
                    -0.49087385940425, -0.10291888921447, 0.11836314681968,
                    0.55527385721943e-4],
            "d1": [1, 1, 1, 2, 2, 3, 3, 4, 4, 4],
            "t1": [1.0, 1.550, 1.700, 0.250, 1.350, 0.0, 1.250, 0.0, 0.7, 5.4],
            "nr2": []},

        "3-5": {
            "nr1": [.25574776844118e1, -.79846357136353e1, .47859131465806e1,
                    -0.73265392369587, 0.13805471345312e1, 0.28349603476365,
                    -0.49087385940425, -0.10291888921447, 0.11836314681968,
                    0.55527385721943e-4],
            "d1": [1, 1, 1, 2, 2, 3, 3, 4, 4, 4],
            "t1": [1.0, 1.550, 1.700, 0.250, 1.350, 0.0, 1.250, 0.0, 0.7, 5.4],
            "nr2": []},

        "3-6": {
            "nr1": [.25574776844118e1, -0.79846357136353e1, 0.47859131465806e1,
                    -0.73265392369587, 0.13805471345312e1, 0.28349603476365,
                    -0.49087385940425, -0.10291888921447, 0.11836314681968,
                    0.55527385721943e-4],
            "d1": [1, 1, 1, 2, 2, 3, 3, 4, 4, 4],
            "t1": [1.0, 1.550, 1.7, 0.250, 1.350, 0.0, 1.250, 0.0, 0.70, 5.4],
            "nr2": []},

        "4-5": {
            "nr1": [0.25574776844118e1, -.79846357136353e1, 0.47859131465806e1,
                    -0.73265392369587, 0.13805471345312e1, 0.28349603476365,
                    -0.49087385940425, -0.10291888921447, 0.11836314681968,
                    0.55527385721943e-4],
            "d1": [1, 1, 1, 2, 2, 3, 3, 4, 4, 4],
            "t1": [1.0, 1.550, 1.7, 0.250, 1.350, 0.0, 1.250, 0.0, 0.7, 5.4],
            "nr2": []},

        "4-6": {
            "nr1": [.25574776844118e1, -0.79846357136353e1, 0.47859131465806e1,
                    -0.73265392369587, 0.13805471345312e1, 0.28349603476365,
                    -0.49087385940425, -0.10291888921447, 0.11836314681968,
                    0.55527385721943e-4],
            "d1": [1, 1, 1, 2, 2, 3, 3, 4, 4, 4],
            "t1": [1.0, 1.550, 1.7, 0.250, 1.350, 0.0, 1.250, 0.0, 0.7, 5.4],
            "nr2": []},

        "5-6": {
            "nr1": [.25574776844118e1, -0.79846357136353e1, 0.47859131465806e1,
                    -0.73265392369587, 0.13805471345312e1, 0.28349603476365,
                    -0.49087385940425, -0.10291888921447, 0.11836314681968,
                    0.55527385721943e-4],
            "d1": [1, 1, 1, 2, 2, 3, 3, 4, 4, 4],
            "t1": [1.0, 1.550, 1.7, 0.250, 1.350, 0.0, 1.250, 0.0, 0.7, 5.400],
            "nr2": []}}

    def __init__(self, **kwargs):
        """Constructor
            To define it need specified composition:
                -componente: array with index of fluids
                -fraccion: molar fraction
           It need to specified state define two properties from this:
                -T: temperature, K
                -rho: density, kg/m3 - mol/dm3
                -P: pressure, Pa
                -v: specific volume, m3/kg
                -h: enthalpy, J/kg
                -s: entropy, J/kgK
                -u: internal energy, J/kg
                -x: quality
            """
        self.kwargs = GERG.kwargs.copy()
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)

        if self.calculable:
            self.status = 1
            self.calculo()
            self.msg = ""

    @property
    def calculable(self):
        if self.kwargs["componente"] and self.kwargs["fraccion"]:
            self._definition = True
        else:
            self._definition = False

        thermo = 0
        for key in ("T", "P", "rho", "v"):
            if self.kwargs[key]:
                thermo += 1
        for key in ("h", "s", "u", "x"):
            if self.kwargs[key] is not None:
                thermo += 1
        return self._definition and thermo >= 2

    def calculo(self):
        import math
        import numpy

        from numpy import exp, log, zeros, r_
        from scipy.constants import R
        from scipy.optimize import fsolve
        from scipy.constants import R, calorie, liter, atm, Btu, lb

        T = self.kwargs["T"]
        rho = self.kwargs["rho"]
        P = self.kwargs["P"]
        v = self.kwargs["v"]
        h = self.kwargs["h"]
        s = self.kwargs["s"]
        u = self.kwargs["u"]
        x = self.kwargs["x"]

        self.comp = []
        for i in self.kwargs["componente"]:
            c = self.componentes[i]
            self.comp.append(c)
        self.id = self.kwargs["componente"]
        self.xi = self.kwargs["fraccion"]
        self.zi = self.xi
        self.phase = None

        # Critic properties for mixture,
        # eq. 7.9, 7.10 pag.125, Tabla 7.10 pag 136
        bt = self.Prop_c["beta_t"]
        bt = bt + [[0 for i in range(21)]]
        bv = self.Prop_c["beta_v"]
        bv = bv + [[0 for i in range(21)]]
        gt = self.Prop_c["gamma_t"]
        gt = gt + [[0 for i in range(21)]]
        gv = self.Prop_c["gamma_v"]
        gv = gv + [[0 for i in range(21)]]
        c_T = zeros((len(self.comp), len(self.comp)))
        c_rho = zeros((len(self.comp), len(self.comp)))
        for i, cmpi in enumerate(self.comp):
            for j, cmpj in enumerate(self.comp):
                if i > j and bt[self.id[j]][self.id[i]] != 0:
                    bt[self.id[i]][self.id[j]] = 1. / bt[self.id[j]][self.id[i]]
                    bv[self.id[i]][self.id[j]] = 1. / bv[self.id[j]][self.id[i]]
                    gt[self.id[i]][self.id[j]] = gt[self.id[j]][self.id[i]]
                    gv[self.id[i]][self.id[j]] = gv[self.id[j]][self.id[i]]
        for i, cmpi in enumerate(self.comp):
            for j, cmpj in enumerate(self.comp):
                c_T[i, j] = 2 * bt[self.id[i]][self.id[j]] * gt[self.id[i]][self.id[j]] * \
                            (cmpi.Tc * cmpj.Tc) ** 0.5
                c_rho[i, j] = 2 * bv[self.id[i]][self.id[j]] * gv[self.id[i]][self.id[j]] / 8. * \
                              (1. / cmpi.rhocL ** (1. / 3) + 1. / cmpj.rhocL ** (1. / 3)) ** 3

        f_T = zeros((len(self.comp), len(self.comp)))
        f_rho = zeros((len(self.comp), len(self.comp)))
        dFT_ik = zeros((len(self.comp), len(self.comp)))
        dFT_ki = zeros((len(self.comp), len(self.comp)))
        dFrho_ik = zeros((len(self.comp), len(self.comp)))
        dFrho_ki = zeros((len(self.comp), len(self.comp)))
        for i, x_i in enumerate(self.xi):
            for j, x_j in enumerate(self.xi):
                f_T[i, j] = x_i * x_j * (x_i + x_j) / (bt[self.id[i]][self.id[j]] ** 2 * x_i + x_j)
                f_rho[i, j] = x_i * x_j * (x_i + x_j) / (bv[self.id[i]][self.id[j]] ** 2 * x_i + x_j)
                dFT_ik[i, j] = x_j * (x_j + x_i) / (bt[self.id[i]][self.id[j]] ** 2 * x_j + x_i) + \
                               x_j * x_i / (bt[self.id[i]][self.id[j]] ** 2 * x_j + x_i) * \
                               (1 - bt[self.id[i]][self.id[j]] ** 2 * (x_j + x_i) / (bt[self.id[i]][self.id[j]] ** 2 * \
                                                                                     x_i + x_j))
                dFrho_ik[i, j] = x_j * (x_j + x_i) / (bv[self.id[i]][self.id[j]] ** 2 * x_i + x_j) + \
                                 x_j * x_i / (bv[self.id[i]][self.id[j]] ** 2 * x_i + x_j) * \
                                 (1 - bv[self.id[i]][self.id[j]] ** 2 * (x_j + x_i) / (bv[self.id[i]][self.id[j]] ** 2 * \
                                                                                       x_i + x_j))
                dFT_ki[j, i] = x_j * (x_j + x_i) / (bt[self.id[i]][self.id[j]] ** 2 * x_j + x_i) + x_j * x_i / \
                               (bt[self.id[i]][self.id[j]] ** 2 * x_j + x_i) * \
                               (1 - (x_j + x_i) / (bt[self.id[i]][self.id[j]] ** 2 * x_j + x_i))
                dFrho_ki[j, i] = x_j * (x_j + x_i) / (bv[self.id[i]][self.id[j]] ** 2 * x_j + x_i) + x_j * x_i / \
                                 (bv[self.id[i]][self.id[j]] ** 2 * x_j + x_i) * \
                                 (1 - (x_j + x_i) / (bv[self.id[i]][self.id[j]] ** 2 * x_j + x_i))

        sumai_v = 0
        sumaij_v = 0
        sumai_T = 0
        sumaij_T = 0
        m = 0
        for i, componentei in enumerate(self.comp):
            sumai_v += self.xi[i] ** 2 / componentei.rhocL
            sumai_T += self.xi[i] ** 2 * componentei.Tc
            m += self.xi[i] * componentei.M
            for j, componentej in enumerate(self.comp):
                if j > i:
                    sumaij_v += c_rho[i, j] * f_rho[i, j]
                    sumaij_T += c_T[i, j] * f_T[i, j]

        self.rhoc = 1. / (sumai_v + sumaij_v) * m # eq.(16)
        self.Tc = sumai_T + sumaij_T  # eq.(17)
        self.M = m / 1000  # kg/mol
        self.R = R / self.M

        Tcxi = []
        rhocxi = []
        for i, componentei in enumerate(self.comp):
            sumav1 = sumat1 = 0
            for k in range(i):
                sumav1 += c_rho[k, i] * dFrho_ki[k, i]
                sumat1 += c_T[k, i] * dFT_ki[k, i]
            sumav2 = sumat2 = 0
            for k in range(i + 1, len(self.xi)):
                sumav2 += c_rho[i, k] * dFrho_ik[i, k]
                sumat2 += c_T[i, k] * dFT_ik[i, k]

            Tcxi.append(2 * self.xi[i] * componentei.Tc + sumat1 + sumat2)
            rhocxi.append(2 * self.xi[i] / componentei.rhocL + sumav1 + sumav2)
        self.Tcxi = Tcxi
        self.rhocxi = rhocxi


        if v and not rho:
            rho = 1. / v

        if T and x is not None:
            pass
        else:
            if T and P:
                rhoo = 1000.
                rho = fsolve(lambda rho: self._solve(rho, T)["P"] - P, rhoo)
            elif T and rho:
                pass
            elif T and h is not None:
                rho = fsolve(lambda rho: self._solve(rho, T)["h"] - h, 200)
            elif T and s is not None:
                rho = fsolve(lambda rho: self._solve(rho, T)["s"] - s, 200)
            elif T and u is not None:
                rho = fsolve(lambda rho: self._solve(rho, T)["u"] - u, 200)
            elif P and rho:
                T = fsolve(lambda T: self._solve(rho, T)["P"] - P, 600)
            elif P and h is not None:
                rho, T = fsolve(lambda par: (
                    self._solve(par[0], par[1])["P"] - P, self._solve(
                        par[0], par[1])["h"] - h), [200, 600])
            elif P and s is not None:
                rho, T = fsolve(lambda par: (
                    self._solve(par[0], par[1])["P"] - P, self._solve(
                        par[0], par[1])["s"] - s), [200, 600])
            elif P and u is not None:
                rho, T = fsolve(lambda par: (
                    self._solve(par[0], par[1])["P"] - P, self._solve(
                        par[0], par[1])["u"] - u), [200, 600])
            elif rho and h is not None:
                T = fsolve(lambda T: self._solve(rho, T)["h"] - h, 600)
            elif rho and s is not None:
                T = fsolve(lambda T: self._solve(rho, T)["s"] - s, 600)
            elif rho and u is not None:
                T = fsolve(lambda T: self._solve(rho, T)["u"] - u, 600)
            elif h is not None and s is not None:
                rho, T = fsolve(lambda par: (
                    self._solve(par[0], par[1])["h"] - h, self._solve(
                        par[0], par[1])["s"] - s), [200, 600])
            elif h is not None and u is not None:
                rho, T = fsolve(lambda par: (
                    self._solve(par[0], par[1])["h"] - h, self._solve(
                        par[0], par[1])["u"] - u), [200, 600])
            elif s is not None and u is not None:
                rho, T = fsolve(lambda par: (
                    self._solve(par[0], par[1])["s"] - s, self._solve(
                        par[0], par[1])["u"] - u), [200, 600])
            else:
                raise IOError

        fio, fiot, fiott, fiod, fiodd, fiodt, fir, firt, firtt, fird, firdd, \
        firdt, firdtt, nfioni, nfirni = self._eq(rho, T)

        # Tabla 7.1 pag 127
        tau = self.Tc / T
        delta = rho / self.rhoc

        self.T = T
        self.rho = rho
        print('rho', self.rho)
        self.v = 1. / rho
        self.P = (1 + delta * fird) * self.R * T * rho
        print('P', self.P)

        self.Z = 1 + delta * fird
        print('compression factor', self.Z)
        self.s = self.R * (tau * (fiot + firt) - fio - fir)
        print('entropy', self.s)
        sr = (tau * firt - fir)
        print('residual entropy', sr)
        visco = self.dimensionless_viscosity(sr)
        CEmix_visco = self.mix_CE_visco()
        self.visco = CEmix_visco * visco
        print("viscosity: ", self.visco)
        self.u = self.R * T * tau * (fiot + firt)
        print('internal energy', self.u)
        self.h = self.R * T * (1 + tau * (fiot + firt) + delta * fird)
        print('enthalpy', self.h)
        self.cp = self.R * (-tau ** 2 * (fiott + firtt) + (1 + delta * fird - delta * tau * firdt) ** 2 / (
                    1 + 2 * delta * fird + delta ** 2 * firdd))
        print('isobaric heat capacity', self.cp)
        self.cv = -self.R * tau ** 2 * (fiott + firtt)
        print('isochoric heat capacity', self.cv)
        self.g = self.R * T * (1 + fio + fir + delta * fird)
        print('Gibbs free energy', self.g)
        self.w = (self.R * T * (1 + 2 * delta * fird + delta ** 2 * firdd - (
                    1 + delta * fird - delta * tau * firdt) ** 2 / tau ** 2 / (fiott + firtt))) ** 0.5
        print('speed of sound', self.w)
        self.mrhog = -1
        self.mrhol = -1
        self.mKi = -1
        self.mxi = -1
        self.myi = -1
        self.mgasvisco = -1
        self.mliqvisco = -1
        self.mQo = -1
        f, FI = self.fug(rho, T)
        self.f = f
        self.FI = FI
        print("fugacity: ", self.FI)
        Ki, xi, yi, Q = self.flash()

    def dimensionless_viscosity(self, sr):
        from numpy import exp, log, zeros, r_
        A = [None] * len(self.comp)
        B = [None] * len(self.comp)
        C = [None] * len(self.comp)
        D = [None] * len(self.comp)
        As = 0
        Bs = 0
        Cs = 0
        Ds = 0
        for i, componentei in enumerate(self.comp):
            A[i] = componentei.A
            B[i] = componentei.B
            C[i] = componentei.C
            D[i] = componentei.D
        for i, componentei in enumerate(self.comp):
            As += A[i] * self.xi[i]
            Bs += B[i] * self.xi[i]
            Cs += C[i] * self.xi[i]
            Ds += D[i] * self.xi[i]
        lnvisco = As + Bs * sr + Cs * sr**2 + Ds * sr**3
        visco = exp(lnvisco)
        return visco
    def mix_CE_visco(self):
        from numpy import exp, log, zeros, r_
        CE_visco = [None] * len(self.comp)
        for i, componentei in enumerate(self.comp):
            Tc = componentei.Tc
            rhoc = componentei.rhocL
            M = componentei.M
            CE_visco[i] = self.CE_viscosity(Tc, rhoc, M)
        fi = zeros((len(self.comp), len(self.comp)))
        for i, componentei in enumerate(self.comp):
            for j, componentej in enumerate(self.comp):
                fi[i][j] = (1 + (CE_visco[i]/CE_visco[j])**(1/2)*(componentej.M/componentei.M)**(1/4))**2/(8*(1 + componentei.M/componentej.M))**(1/2)
        mixCE_visco = 0
        sum1 = zeros(len(self.comp))
        for i, componentei in enumerate(self.comp):
            for j, componentej in enumerate(self.comp):
                sum1[i] += self.xi[i] * fi[i][j]
        for i, componentei in enumerate(self.comp):
            mixCE_visco += self.xi[i] * CE_visco[i] / sum1[i]
        return mixCE_visco

    def CE_viscosity(self, Tc, rhocL, M):
        from numpy import exp, log, zeros, r_
        import math
        kB = 1.3807e-23 #Boltzmann's constant
        NA = 6.02214076e23 #Avogadro's number
        epson = Tc/1.2593*kB
        Tr = kB * self.T/epson #reduced temperature
        A = 1.16145
        B = 0.14874
        C = 0.52487
        D = 0.77320
        E = 2.16178
        F = 2.43787
        G = -6.435e-4
        W = -0.76830
        P = 7.27371
        S = 18.0323
        R = -6.435e-4
        omiga = (A / Tr**B) + C / exp(D * Tr) + E / exp(F * Tr) + G * Tr**B * math.sin(S * Tr**W - P) # reduced collision integral
        omiga2 = (A / Tr**B) + ((C/exp(D*Tr))) + (E/exp(F*Tr)) + R*Tr**B*math.sin(S*Tr**W-P)
        delta = (0.3189 / rhocL )**(1/3)
        CE_visco = 2.6693e-5 * 0.006747 * 0.10575054361 * (M * self.T )**(1/2) / delta**2 / omiga2
        return CE_visco


    def fug(self, rho, T, nfirni=None):
        from scipy.constants import R, calorie, liter, atm, Btu, lb
        from numpy import exp, log, zeros, r_

        # if not nfirni:
        tau = self.Tc / T
        delta = rho / self.rhoc
        fio, fiot, fiott, fiod, fiodd, fiodt, fir, firt, firtt, fird, firdd, \
        firdt, firdtt, nfioni, nfirni = self._eq(rho, T)
        f = []
        FI = []
        for xi, dn in zip(self.xi, nfirni):
            R_atml = R / liter / atm
            f.append(xi * rho / self.M * R_atml * T * exp(dn))  # Fugacity of Component
            FI.append(exp(dn - log(1 + delta * fird)))
        return f, FI

    def _eq(self, rho, T):
        tau = self.Tc / T
        delta = rho / self.rhoc
        fio, fiot, fiott, fiod, fiodd, fiodt, nfioni = self._phi0(tau, delta)
        fir, firt, firtt, fird, firdd, firdt, firdtt, nfirni = self._phir(tau, delta)
        return (fio, fiot, fiott, fiod, fiodd, fiodt, fir, firt, firtt, fird,
                firdd, firdt, firdtt, nfioni, nfirni)

    def _solve(self, rho, T):
        # print(rho)
        tau = self.Tc / T
        delta = rho / self.rhoc
        fio, fiot, fiott, fiod, fiodd, fiodt, nfioni = self._phi0(tau, delta)
        fir, firt, firtt, fird, firdd, firdt, firdtt, nfirni = self._phir(tau, delta)
        propiedades = {}
        propiedades["P"] = (1 + delta * fird) * self.R * T * rho
        propiedades["s"] = self.R * (tau * (fiot + firt) - fio - fir)
        propiedades["u"] = self.R * T * tau * (fiot + firt)
        propiedades["h"] = self.R * T * (1 + tau * (fiot + firt) + delta * fird)
        return propiedades

    def _phi0(self, tau, delta):
        """Contribución ideal de la energía libre de Helmholtz eq. 7.5"""
        from numpy import exp, log, zeros, r_
        # Table 7.5, 7.6 page 110
        fio = 0
        fiot = 0
        fiott = 0
        fiod = 0
        fiodd = 0
        fiodt = 0
        nfioni = []  # ðnao/ðni
        for i, componente in enumerate(self.comp):
            deltai = delta
            taui = componente.Tc * tau / self.Tc
            deltar = self.rhoc / componente.rhocL
            taur = componente.Tc / self.Tc
            fio_, fiot_, fiott_, fiod_, fiodd_, fiodt_ = self._f0(
                componente, taui, deltai)
            fio += self.xi[i] * (fio_ + log(self.xi[i]))  # page110 eq.(7.20a)
            fiot += self.xi[i] * fiot_ * taur  # eq.(7.20e)
            fiott += self.xi[i] * fiott_ * (taur ** 2)  # eq.(7.20f)
            fiod += self.xi[i] * fiod_ * deltar  # eq.(7.20b)
            fiodd += self.xi[i] * fiodd_ * (deltar ** 2)  # eq.(7.20c)
            fiodt += self.xi[i] * fiodt_ * deltar * taur  # eq.(7.20d)
            nfioni.append(fio_ + 1 + log(self.xi[i]))
        return fio, fiot, fiott, fiod, fiodd, fiodt, nfioni

    def _f0(self, componente, taui, deltai):
        from numpy import exp, log, zeros, r_
        import math
        RR = 8.314510 / 8.314472
        cp = componente.GERG["cp"]
        if len(cp['ao_sinh']) == 2 and len(cp['ao_cosh']) == 2:
            fio_ = log(deltai) + RR * (cp['ao_pow'][0] + cp['ao_pow'][1] * taui + \
                                       cp['ao_log'][1] * log(taui) + \
                                       cp['ao_sinh'][0] * log(abs(math.sinh(cp['sinh'][0] * taui))) + \
                                       cp['ao_sinh'][1] * log(abs(math.sinh(cp['sinh'][1] * taui))) - \
                                       cp['ao_cosh'][0] * log(math.cosh(cp['cosh'][0] * taui)) - \
                                       cp['ao_cosh'][1] * log(math.cosh(cp['cosh'][1] * taui)))
            fiot_ = RR * (cp['ao_pow'][1] + \
                          cp['ao_log'][1] * (1 / taui) + \
                          cp['ao_sinh'][0] * cp['sinh'][0] / math.tanh(cp['sinh'][0] * taui) + \
                          cp['ao_sinh'][1] * cp['sinh'][1] / math.tanh(cp['sinh'][1] * taui) - \
                          cp['ao_cosh'][0] * cp['cosh'][0] * math.tanh(cp['cosh'][0] * taui) - \
                          cp['ao_cosh'][1] * cp['cosh'][1] * math.tanh(cp['cosh'][1] * taui))
            fiott_ = RR * (-cp['ao_log'][1] * (1 / taui) ** 2 - \
                           cp['ao_sinh'][0] * (cp['sinh'][0] ** 2) / (math.sinh(cp['sinh'][0] * taui)) ** 2 - \
                           cp['ao_sinh'][1] * (cp['sinh'][1] ** 2) / (math.sinh(cp['sinh'][1] * taui)) ** 2 - \
                           cp['ao_cosh'][0] * (cp['cosh'][0] ** 2) / (math.cosh(cp['cosh'][0] * taui)) ** 2 - \
                           cp['ao_cosh'][1] * (cp['cosh'][1] ** 2) / (math.cosh(cp['cosh'][1] * taui)) ** 2)
            fiod_ = 1 / deltai
            fiodd_ = -(1 / deltai) ** 2
            fiodt_ = 0

        elif len(cp['ao_sinh']) == 2 and len(cp['ao_cosh']) == 1:
            fio_ = log(deltai) + RR * (cp['ao_pow'][0] + cp['ao_pow'][1] * taui + \
                                       cp['ao_log'][1] * log(taui) + \
                                       cp['ao_sinh'][0] * log(abs(math.sinh(cp['sinh'][0] * taui))) + \
                                       cp['ao_sinh'][1] * log(abs(math.sinh(cp['sinh'][1] * taui))) - \
                                       cp['ao_cosh'][0] * log(abs(math.cosh(cp['cosh'][0] * taui))))
            fiot_ = RR * (cp['ao_pow'][1] + \
                          cp['ao_log'][1] * (1 / taui) + \
                          cp['ao_sinh'][0] * cp['sinh'][0] / math.tanh(cp['sinh'][0] * taui) + \
                          cp['ao_sinh'][1] * cp['sinh'][1] / math.tanh(cp['sinh'][1] * taui) - \
                          cp['ao_cosh'][0] * cp['cosh'][0] * math.tanh(cp['cosh'][0] * taui))
            fiott_ = RR * (-cp['ao_log'][1] * (1 / taui) ** 2 - \
                           cp['ao_sinh'][0] * (cp['sinh'][0] ** 2) / (math.sinh(cp['sinh'][0] * taui)) ** 2 - \
                           cp['ao_sinh'][1] * (cp['sinh'][1] ** 2) / (math.sinh(cp['sinh'][1] * taui)) ** 2 - \
                           cp['ao_cosh'][0] * (cp['cosh'][0] ** 2) / (math.cosh(cp['cosh'][0] * taui)) ** 2)
            fiod_ = 1 / deltai
            fiodd_ = -(1 / deltai) ** 2
            fiodt_ = 0

        elif len(cp['ao_sinh']) == 1 and len(cp['ao_cosh']) == 1:
            fio_ = log(deltai) + RR * (cp['ao_pow'][0] + cp['ao_pow'][1] * taui + \
                                       cp['ao_log'][1] * log(taui) + \
                                       cp['ao_sinh'][0] * log(abs(math.sinh(cp['sinh'][0] * taui))) - \
                                       cp['ao_cosh'][0] * log(abs(math.cosh(cp['cosh'][0] * taui))))
            fiot_ = RR * (cp['ao_pow'][1] + \
                          cp['ao_log'][1] * (1 / taui) + \
                          cp['ao_sinh'][0] * cp['sinh'][0] / math.tanh(cp['sinh'][0] * taui) - \
                          cp['ao_cosh'][0] * cp['cosh'][0] * math.tanh(cp['cosh'][0] * taui))
            fiott_ = RR * (-cp['ao_log'][1] * (1 / taui) ** 2 - \
                           cp['ao_sinh'][0] * (cp['sinh'][0] ** 2) / (math.sinh(cp['sinh'][0] * taui)) ** 2 - \
                           cp['ao_cosh'][0] * (cp['cosh'][0] ** 2) / (math.cosh(cp['cosh'][0] * taui)) ** 2)
            fiod_ = 1 / deltai
            fiodd_ = -(1 / deltai) ** 2
            fiodt_ = 0
        elif len(cp['ao_sinh']) == 1 and len(cp['ao_cosh']) == 0:
            fio_ = log(deltai) + RR * (cp['ao_pow'][0] + cp['ao_pow'][1] * taui + \
                                       cp['ao_log'][1] * log(taui) + \
                                       cp['ao_sinh'][0] * log(abs(math.sinh(cp['sinh'][0] * taui))))
            fiot_ = RR * (cp['ao_pow'][1] + \
                          cp['ao_log'][1] * (1 / taui) + \
                          cp['ao_sinh'][0] * cp['sinh'][0] / math.tanh(cp['sinh'][0] * taui))
            fiott_ = RR * (-cp['ao_log'][1] * (1 / taui) ** 2 - \
                           cp['ao_sinh'][0] * (cp['sinh'][0] ** 2) / (math.sinh(cp['sinh'][0] * taui)) ** 2)
            fiod_ = 1 / deltai
            fiodd_ = -(1 / deltai) ** 2
            fiodt_ = 0
        elif len(cp['ao_sinh']) == 0 and len(cp['ao_cosh']) == 0:
            fio_ = log(deltai) + RR * (cp['ao_pow'][0] + cp['ao_pow'][1] * taui + \
                                       cp['ao_log'][1] * log(taui))
            fiot_ = RR * (cp['ao_pow'][1] + \
                          cp['ao_log'][1] * (1 / taui))
            fiott_ = RR * (-cp['ao_log'][1] * (1 / taui) ** 2)
            fiod_ = 1 / deltai
            fiodd_ = -(1 / deltai) ** 2
            fiodt_ = 0
        return fio_, fiot_, fiott_, fiod_, fiodd_, fiodt_

    def _phir(self, tau, delta):
        from numpy import exp, log, zeros, r_
        """Contribución residual de la energía libre de Helmholtz eq. 7.7"""
        # Table 7.5
        fir = 0
        firt = 0
        firtt = 0
        fird = 0
        firdd = 0
        firdt = 0
        firdtt = 0
        firxi = []
        firxixj = zeros((len(self.comp), len(self.comp)))
        firdxi = []
        firtxi = []
        fir0 = [None] * len(self.comp)

        for i, componente in enumerate(self.comp):
            fir_, firt_, firtt_, fird_, firdd_, firdt_ = \
                self._Helmholtz(componente, tau, delta)
            fir += self.xi[i] * fir_  # eq(7.21a)
            firt += self.xi[i] * firt_  # eq(7.21e)
            firtt += self.xi[i] * firtt_  # eq(7.21f)
            fird += self.xi[i] * fird_  # eq(7.21b)
            firdd += self.xi[i] * firdd_  # eq(7.21c)
            firdt += self.xi[i] * firdt_  # eq(7.21d)

            fir0[i] = fir_

            firxi.append(fir_)  # eq(7.21g)
            firdxi.append(fird_)  # eq(7.21j)
            firtxi.append(firt_)  # eq(7.21k)

        # Contribución residual cruzada eq 7.8
        for i, x_i in enumerate(self.xi):
            for j, x_j in enumerate(self.xi):
                if j > i:
                    (fir_, firt_, firtt_, fird_, firdd_, firdt_, firdtt_, B_,
                     C_) = self._phijr(i, j, tau, delta)
                    fir += x_i * x_j * self.Fij[self.id[i]][self.id[j]] * fir_
                    firt += x_i * x_j * self.Fij[self.id[i]][self.id[j]] * firt_
                    firtt += x_i * x_j * self.Fij[self.id[i]][self.id[j]] * firtt_
                    fird += x_i * x_j * self.Fij[self.id[i]][self.id[j]] * fird_
                    firdd += x_i * x_j * self.Fij[self.id[i]][self.id[j]] * firdd_
                    firdt += x_i * x_j * self.Fij[self.id[i]][self.id[j]] * firdt_

                if i != j and j > i:
                    (fir_, firt_, firtt_, fird_, firdd_, firdt_, firdtt_, B_,
                     C_) = self._phijr(i, j, tau, delta)  # 写了
                    firxi[i] += x_j * self.Fij[self.id[i]][self.id[j]] * fir_
                    firxixj[i, j] = self.Fij[self.id[i]][self.id[j]] * fir_
                    firdxi[i] += x_j * self.Fij[self.id[i]][self.id[j]] * fird_
                    firtxi[i] += x_j * self.Fij[self.id[i]][self.id[j]] * firt_

                if i != j and j < i:
                    (fir_, firt_, firtt_, fird_, firdd_, firdt_, firdtt_, B_,
                     C_) = self._phijr(j, i, tau, delta)  # 写了
                    firxi[i] += x_j * self.Fij[self.id[j]][self.id[i]] * fir_
                    firxixj[i, j] = self.Fij[self.id[j]][self.id[i]] * fir_
                    firdxi[i] += x_j * self.Fij[self.id[j]][self.id[i]] * fird_
                    firtxi[i] += x_j * self.Fij[self.id[j]][self.id[i]] * firt_

        suma = 0
        suma_rho = 0
        suma_T = 0

        rhocxi = [None] * len(self.comp)
        for i, componente in enumerate(self.comp):
            rhocxi[i] = -self.rhoc ** 2 * self.rhocxi[i]
        for i, x_i in enumerate(self.xi):
            suma += x_i * firxi[i]
            suma_rho += x_i * rhocxi[i]
            suma_T += x_i * self.Tcxi[i]

        n_rhocni = []
        n_Tcni = []
        for i, x_i in enumerate(self.xi):
            n_rhocni.append(rhocxi[i] - suma_rho)
            n_Tcni.append(self.Tcxi[i] - suma_T)

        n_firni = []  # ðar/ðni
        nfirni = []  # ðnar/ðni
        for i, componente in enumerate(self.comp):
            n_firni.append(delta * fird * (1 - 1. / self.rhoc * n_rhocni[i]) + \
                           tau * firt / self.Tc * n_Tcni[i] + firxi[i] - suma)
            nfirni.append(fir + n_firni[i])
        return fir, firt, firtt, fird, firdd, firdt, firdtt, nfirni

    def _Helmholtz(self, componente, taui, deltai):
        from numpy import exp, log, zeros, r_
        nr1 = componente.GERG["nr1"]
        d1 = componente.GERG["d1"]
        t1 = componente.GERG["t1"]
        nr2 = componente.GERG["nr2"]
        d2 = componente.GERG["d2"]
        t2 = componente.GERG["t2"]
        c2 = componente.GERG["c2"]
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        sum4 = 0.0
        sum5 = 0.0
        sum6 = 0.0
        sum7 = 0.0
        sum8 = 0.0
        sum9 = 0.0
        sum10 = 0.0
        sum11 = 0.0
        sum12 = 0.0
        for i in range(len(nr1)):
            sum1 = nr1[i] * (deltai ** d1[i]) * (taui ** t1[i]) + sum1
        for i in range(len(nr2)):
            sum2 = nr2[i] * (deltai ** d2[i]) * (taui ** t2[i]) * exp(-(deltai ** c2[i])) + sum2
        fir_ = sum1 + sum2

        for i in range(len(nr1)):
            sum3 = nr1[i] * t1[i] * (deltai ** d1[i]) * (taui ** (t1[i] - 1)) + sum3
        # print('sum3', sum3)
        for i in range(len(nr2)):
            sum4 = nr2[i] * t2[i] * (deltai ** d2[i]) * (taui ** (t2[i] - 1)) * exp(-deltai ** c2[i]) + sum4
        firt_ = sum3 + sum4

        for i in range(len(nr1)):
            sum5 = nr1[i] * t1[i] * (t1[i] - 1) * (deltai ** d1[i]) * (taui ** (t1[i] - 2)) + sum5
        for i in range(len(nr2)):
            sum6 = nr2[i] * t2[i] * (t2[i] - 1) * (deltai ** d2[i]) * (taui ** (t2[i] - 2)) * exp(
                -(deltai ** c2[i])) + sum6
        firtt_ = sum5 + sum6

        for i in range(len(nr1)):
            sum7 = nr1[i] * d1[i] * (deltai ** (d1[i] - 1)) * (taui ** t1[i]) + sum7
        for i in range(len(nr2)):
            sum8 = nr2[i] * (deltai ** (d2[i] - 1)) * (d2[i] - c2[i] * (deltai ** c2[i])) * (taui ** t2[i]) * exp(
                -(deltai ** c2[i])) + sum8
        fird_ = sum7 + sum8

        for i in range(len(nr1)):
            sum9 = nr1[i] * d1[i] * (d1[i] - 1) * (deltai ** (d1[i] - 2)) * (taui ** t1[i]) + sum9
        for i in range(len(nr2)):
            sum10 = nr2[i] * (deltai ** (d2[i] - 2)) * ((d2[i] - c2[i] * (deltai ** c2[i])) * \
                                                        (d2[i] - 1 - c2[i] * (deltai ** c2[i])) - (c2[i] ** 2) * (
                                                                    deltai ** c2[i])) * \
                    (taui ** t2[i]) * exp(-(deltai ** c2[i])) + sum10
        firdd_ = sum9 + sum10

        for i in range(len(nr1)):
            sum11 = nr1[i] * d1[i] * t1[i] * (deltai ** (d1[i] - 1)) * (taui ** (t1[i] - 1)) + sum11
        for i in range(len(nr2)):
            sum12 = nr2[i] * t2[i] * (deltai ** (d2[i] - 1)) * (d2[i] - c2[i] * (deltai ** c2[i])) * \
                    (taui ** (t2[i] - 1)) * exp(-(deltai ** c2[i])) + sum12
        firdt_ = sum11 + sum12

        return fir_, firt_, firtt_, fird_, firdd_, firdt_

    def _phijr(self, i, j, tau, delta):
        from numpy import exp, log, zeros, r_
        i = self.id[i]
        j = self.id[j]
        txt = str(i) + "-" + str(j)
        constants = self.fir_ij.get(txt, 0)
        fir = 0
        fird = 0
        firdd = 0
        firt = 0
        firtt = 0
        firdt = 0
        firdtt = 0
        B = 0
        C = 0
        if constants:
            delta_0 = 1e-50
            nr1 = constants["nr1"]
            d1 = constants["d1"]
            t1 = constants["t1"]
            for i in range(len(constants.get("nr1", []))):
                # Polinomial terms
                fir += nr1[i] * delta ** d1[i] * tau ** t1[i]
                fird += nr1[i] * d1[i] * delta ** (d1[i] - 1) * tau ** t1[i]
                firdd += nr1[i] * d1[i] * (d1[i] - 1) * delta ** (d1[i] - 2) * tau ** t1[i]
                firt += nr1[i] * t1[i] * delta ** d1[i] * tau ** (t1[i] - 1)
                firtt += nr1[i] * t1[i] * (t1[i] - 1) * delta ** d1[i] * tau ** (t1[i] - 2)
                firdt += nr1[i] * t1[i] * d1[i] * delta ** (d1[i] - 1) * tau ** (t1[i] - 1)
                firdtt += nr1[i] * t1[i] * d1[i] * (t1[i] - 1) * delta ** (d1[i] - 1) * tau ** (t1[i] - 2)
                B += nr1[i] * d1[i] * delta_0 ** (d1[i] - 1) * tau ** t1[i]
                C += nr1[i] * d1[i] * (d1[i] - 1) * delta_0 ** (d1[i] - 2) * tau ** t1[i]

            nr2 = constants["nr2"]
            if nr2:
                d2 = constants["d2"]
                t2 = constants["t2"]
                n2 = constants["n2"]
                e2 = constants["e2"]
                b2 = constants["b2"]
                g2 = constants["g2"]
                for i in range(len(constants.get("nr2", []))):

                    fir += nr2[i] * delta ** d2[i] * tau ** t2[i] * \
                           exp(-n2[i] * (delta - e2[i]) ** 2 - b2[i] * (delta - g2[i]))
                    fird += nr2[i] * delta ** d2[i] * tau ** t2[i] * \
                            exp(-n2[i] * (delta - e2[i]) ** 2 - b2[i] * (delta - g2[i])) * \
                            (d2[i] / delta - 2 * n2[i] * (delta - e2[i]) - b2[i])
                    firdd += nr2[i] * delta ** d2[i] * tau ** t2[i] * \
                             exp(-n2[i] * (delta - e2[i]) ** 2 - b2[i] * (delta - g2[i])) * \
                             ((d2[i] / delta - 2 * n2[i] * (delta - e2[i]) - b2[i]) ** 2 - d2[i] / delta ** 2 - 2 * n2[
                                 i])
                    firt += nr2[i] * t2[i] * delta ** d2[i] * tau ** (t2[i] - 1) * \
                            exp(-n2[i] * (delta - e2[i]) ** 2 - b2[i] * (delta - g2[i]))
                    firtt += nr2[i] * delta ** d2[i] * tau ** (t2[i] - 2) * t2[i] * (t2[i] - 1) * \
                             exp(-n2[i] * (delta - e2[i]) ** 2 - b2[i] * (delta - g2[i]))
                    firdt += nr2[i] * t2[i] * delta ** d2[i] * tau ** (t2[i] - 1) * \
                             exp(-n2[i] * (delta - e2[i]) ** 2 - b2[i] * (delta - g2[i])) * \
                             (d2[i] / delta - 2 * n2[i] * (delta - e2[i]) - b2[i])
                    firdtt += nr2[i] * delta ** d2[i] * tau ** t2[i] * \
                              exp(-n2[i] * (delta - e2[i]) ** 2 - b2[i] * (tau - g2[i]) ** (2)) * \
                              ((t2[i] / tau - 2 * b2[i] * (tau - g2[i])) ** (2) - t2[i] / tau ** 2 - 2 * b2[i]) * \
                              (d2[i] / delta - 2 * n2[i] * (delta - e2[i]))
                    B += nr2[i] * delta_0 ** d2[i] * tau ** t2[i] * \
                         exp(-n2[i] * (delta_0 - e2[i]) ** 2 - b2[i] * (tau - g2[i]) ** (2)) * \
                         (d2[i] / delta_0 - 2 * n2[i] * (delta_0 - e2[i]))
                    C += nr2[i] * tau ** t2[i] * \
                         exp(-n2[i] * (delta_0 - e2[i]) ** 2 - b2[i] * (tau - g2[i]) ** (2)) * \
                         (-2 * n2[i] * delta_0 ** d2[i] + 4 * n2[i] ** 2 * delta_0 ** d2[i] * (delta_0 - e2[i]) ** 2 -
                          4 * d2[i] * n2[i] * delta_0 ** 2 * (delta_0 - e2[i]) + d2[i] * 2 * delta_0)

        return fir, firt, firtt, fird, firdd, firdt, firdtt, B, C

    def flash(self):
        """Cálculo de los coeficientes de reparto entre fases"""
        from numpy import exp, log, zeros, r_
        from scipy.optimize import fsolve
        from scipy.constants import R
        # Estimación inicial de K mediante correlación wilson Eq 5.61 Pag 82
        Ki = []
        for componente in self.comp:
            Ki.append(componente.Pc / self.P * exp(5.373 * (1. + componente.f_acent) * (1. - componente.Tc / self.T)))
        Kinew = [None] * len(Ki)

        Qmin = 1. / (1 - max(Ki))
        Qmax = 1. / (1 - min(Ki))

        Qo = 0.5

        def f(Q):
            sum = 0
            for i, x in enumerate(self.zi):
                sum += x * (Ki[i] - 1.) / (1. - Q + Q * Ki[i])
            return sum

        kz = 0
        zk = 0
        for i, x in enumerate(self.zi):
            kz += Ki[i] * self.zi[i]

        for i, x in enumerate(self.zi):
            zk += self.zi[i] / Ki[i]
        print('Ki', Ki)
        if max(Ki) < 1:
            xi = self.zi
            yi = self.zi
            Q = 1
            print('liquid only')
        if min(Ki) > 1:
            xi = self.zi
            yi = self.zi
            Q = 0
            print('vapor only')
        else:
            # two phases
            Qo = 0.6
            # self.zi = self.xi
            while True:
                Qmin = 1. / (1 - max(Ki))
                Qmax = 1. / (1 - min(Ki))
                xi = []  # liquid
                yi = []  # vapor
                # print('self.zi', self.zi)
                for i, fraccion in enumerate(self.zi):
                    xi.append(fraccion / (1 + Qo * (Ki[i] - 1)))
                    yi.append(fraccion * Ki[i] / (1 + Qo * (Ki[i] - 1)))
                rhorg, Trg, Tcxig, rhocxig, m = self.reducingF(yi)
                rhog0 = 200.
                self.Tc = Trg
                self.rhoc = rhorg
                self.M = m / 1000  # kg/mol
                self.R = R / self.M
                self.xi = yi
                self.rhocxi = rhocxig
                self.Tcxi = Tcxig
                rhog = fsolve(lambda rhog: self._solve(rhog, self.T)["P"] - self.P, rhog0)
                rhorg, Trg, Tcxig, rhocxig, m = self.reducingF(yi)
                # rhog = 2.
                self.Tc = Trg
                self.rhoc = rhorg
                self.M = m / 1000  # kg/mol
                self.R = R / self.M
                self.xi = yi
                self.rhocxi = rhocxig
                self.Tcxi = Tcxig
                self.phase = 0
                fio, fiot, fiott, fiod, fiodd, fiodt, fir, firt, firtt, fird, firdd, \
                firdt, firdtt, nfioni, nfirni = self._eq(rhog, self.T)

                # Tabla 7.1 pag 127
                tau = self.Tc / self.T
                sr_g = (tau * firt - fir)
                fiv, Fiv = self.fug(rhog, self.T)
                print('Fiv', Fiv)

                rhorl, Trl, Tcxil, rhocxil, m = self.reducingF(xi)
                rhol0 = 1000.
                self.Tc = Trl
                self.rhoc = rhorl
                self.M = m / 1000  # kg/mol
                self.R = R / self.M
                self.xi = xi
                self.rhocxi = rhocxil
                self.Tcxi = Tcxil


                rhol = fsolve(lambda rhol: self._solve(rhol, self.T)["P"] - self.P, rhol0)

                rhorl, Trl, Tcxil, rhocxil, m = self.reducingF(xi)
                self.Tc = Trl
                self.rhoc = rhorl
                self.M = m / 1000  # kg/mol
                self.R = R / self.M
                self.xi = xi
                fio, fiot, fiott, fiod, fiodd, fiodt, fir, firt, firtt, fird, firdd, \
                firdt, firdtt, nfioni, nfirni = self._eq(rhol, self.T)

                # Tabla 7.1 pag 127
                tau = self.Tc / self.T
                sr_l = (tau * firt - fir)
                self.rhocxi = rhocxil
                self.Tcxi = Tcxil
                self.phase = 1
                fil, Fil = self.fug(rhol, self.T)
                print('Fil', Fil)

                # K value iteration method
                for i in range(len(Ki)):
                    Kinew[i] = Fil[i] / Fiv[i]
                print('Kinew', Kinew)
                sumk = 0
                for i in range(len(Ki)):
                    sumk += (Kinew[i] / Ki[i] - 1) ** 2
                if sumk < 0.0000000001:
                    print('two phase')
                    print('vapor fraction', Qo)
                    self.mQo = Qo
                    print('gas density', rhog)
                    self.mrhog = rhog
                    print('liquid density', rhol)
                    self.mrhol = rhol
                    print('gas fraction', yi)
                    self.myi = yi
                    print('liquid fraction', xi)
                    self.mxi = xi
                    print('K-facor', Kinew)
                    self.mKi = Kinew
                    self.xi = yi
                    viscog = self.dimensionless_viscosity(sr_g)
                    CEmix_viscog = self.mix_CE_visco()
                    viscog = viscog * CEmix_viscog
                    print("gas viscosity: ", viscog, "mPa.s")
                    self.mgasvisco = viscog
                    self.xi = xi
                    viscol = self.dimensionless_viscosity(sr_l)
                    CEmix_viscol = self.mix_CE_visco()
                    viscol = viscol * CEmix_viscol
                    print("liquid viscosity: ", viscol, "mPa.s")
                    self.mliqvisco = viscol
                    break
                else:
                    Ki = Kinew.copy()
                    Qmin = 1. / (1 - max(Ki))
                    Qmax = 1. / (1 - min(Ki))
                    if Qo < Qmin:
                        Qo = 0
                        Ki = [None] * len(Fiv)
                        for i in range(len(Fil)):
                            Ki[i] = self.FI[i] / Fiv[i]
                        xi = [None] * len(Ki)  # liquid
                        yi = []  # vapor
                        for i, fraccion in enumerate(self.zi):
                            xi[i] = self.zi[i]
                            yi.append(fraccion * Ki[i] / (1 + Qo * (Ki[i] - 1)))

                    elif Qo > Qmax:
                        Ki = [None] * len(Fiv)
                        for i in range(len(Fil)):
                            Ki[i] = Fil[i] / self.FI[i]
                        Qo = 1
                        xi = []  # liquid
                        yi = [None] * len(Ki)  # vapor
                        for i, fraccion in enumerate(self.zi):
                            xi.append(fraccion / (1 + Qo * (Ki[i] - 1)))
                            yi[i] = self.zi[i]
                    if abs(rhog - rhol) < 0.00000001:
                        print('single phase')
                        print('rho', self.rho)
                        break

                    #determine single region
                    kz = 0
                    zk = 0
                    # print('Qo', Qo)
                    for i, x in enumerate(self.zi):
                        kz += Ki[i] * self.zi[i]
                    for i, x in enumerate(self.zi):
                        zk += self.zi[i] / Ki[i]
                    # print('zk', zk)
                    if kz < 1:
                        xi = self.zi
                        yi = self.zi
                        Q = 1
                        print('liquid only')
                        break
                    elif zk < 1:
                        xi = self.zi
                        yi = self.zi
                        Q = 0
                        print('vapor only')
                        break
                Qo = self.NewtonVF(Ki, Qo)
                if Qo > 1:
                    Qo = 1
                if Qo < 0:
                    Qo = 0
        return Kinew, xi, yi, Qo

    def NewtonVF(self, Ki, Q):
        fQ = 0
        dfQ = 0
        Qold = Q
        for i in range(1, 1000):
            fQ = self.fQ(Ki, Qold)
            dfQ = self.dfQ(Ki, Qold)
            Qnew = Qold - fQ / dfQ
            if abs(Qnew - Qold) < 0.000001:
                return Qnew
            else:
                Qold = Qnew
        if i > 1000:
            Qnew = -0.999
            print('Method did not converge')

    def fQ(self, Ki, Q):
        fQ = 0
        for i in range(len(Ki)):
            fQ = fQ + self.zi[i] * (Ki[i] - 1) / ((Ki[i] - 1) * Q + 1)
        return fQ

    def dfQ(self, Ki, Q):
        dfQ = 0
        for i in range(len(Ki)):
            dfQ = dfQ + self.zi[i] * (Ki[i] - 1) ** 2 / (((Ki[i] - 1) * Q + 1) ** 2)
        dfQ = -dfQ
        return dfQ

    def reducingF(self, xy):
        from numpy import exp, log, zeros, r_
        bt = self.Prop_c["beta_t"]
        bt = bt + [[0 for i in range(21)]]
        bv = self.Prop_c["beta_v"]
        bv = bv + [[0 for i in range(21)]]
        gt = self.Prop_c["gamma_t"]
        gt = gt + [[0 for i in range(21)]]
        gv = self.Prop_c["gamma_v"]
        gv = gv + [[0 for i in range(21)]]
        c_T = zeros((len(self.comp), len(self.comp)))
        c_rho = zeros((len(self.comp), len(self.comp)))
        for i, cmpi in enumerate(self.comp):
            for j, cmpj in enumerate(self.comp):
                if i > j and bt[self.id[j]][self.id[i]] != 0:
                    bt[self.id[i]][self.id[j]] = 1. / bt[self.id[j]][self.id[i]]
                    bv[self.id[i]][self.id[j]] = 1. / bv[self.id[j]][self.id[i]]
                    gt[self.id[i]][self.id[j]] = gt[self.id[j]][self.id[i]]
                    gv[self.id[i]][self.id[j]] = gv[self.id[j]][self.id[i]]
        for i, cmpi in enumerate(self.comp):
            for j, cmpj in enumerate(self.comp):
                c_T[i, j] = 2 * bt[self.id[i]][self.id[j]] * gt[self.id[i]][self.id[j]] * \
                            (cmpi.Tc * cmpj.Tc) ** 0.5
                c_rho[i, j] = 2 * bv[self.id[i]][self.id[j]] * gv[self.id[i]][self.id[j]] / 8. * \
                              (1. / cmpi.rhocL ** (1. / 3) + 1. / cmpj.rhocL ** (1. / 3)) ** 3
        f_T = zeros((len(self.comp), len(self.comp)))
        f_rho = zeros((len(self.comp), len(self.comp)))
        dFT_ik = zeros((len(self.comp), len(self.comp)))
        dFT_ki = zeros((len(self.comp), len(self.comp)))
        dFrho_ik = zeros((len(self.comp), len(self.comp)))
        dFrho_ki = zeros((len(self.comp), len(self.comp)))
        for i, x_i in enumerate(xy):
            for j, x_j in enumerate(xy):
                f_T[i, j] = x_i * x_j * (x_i + x_j) / (bt[self.id[i]][self.id[j]] ** 2 * x_i + x_j)
                f_rho[i, j] = x_i * x_j * (x_i + x_j) / (bv[self.id[i]][self.id[j]] ** 2 * x_i + x_j)
                dFT_ik[i, j] = x_j * (x_j + x_i) / (bt[self.id[i]][self.id[j]] ** 2 * x_i + x_j) + \
                               x_j * x_i / (bt[self.id[i]][self.id[j]] ** 2 * x_i + x_j) * \
                               (1 - bt[self.id[i]][self.id[j]] ** 2 * (x_j + x_i) / (bt[self.id[i]][self.id[j]] ** 2 * \
                                                                                     x_i + x_j))
                dFrho_ik[i, j] = x_j * (x_j + x_i) / (bv[self.id[i]][self.id[j]] ** 2 * x_i + x_j) + \
                                 x_j * x_i / (bv[self.id[i]][self.id[j]] ** 2 * x_i + x_j) * \
                                 (1 - bv[self.id[i]][self.id[j]] ** 2 * (x_j + x_i) / (bv[self.id[i]][self.id[j]] ** 2 * \
                                                                                       x_i + x_j))
                dFT_ki[j, i] = x_j * (x_j + x_i) / (bt[self.id[j]][self.id[i]] ** 2 * x_j + x_i) + x_j * x_i / \
                               (bt[self.id[j]][self.id[i]] ** 2 * x_j + x_i) * \
                               (1 - (x_j + x_i) / (bt[self.id[j]][self.id[i]] ** 2 * x_j + x_i))

                dFrho_ki[j, i] = x_j * (x_j + x_i) / (bv[self.id[j]][self.id[i]] ** 2 * x_j + x_i) + x_j * x_i / \
                                 (bv[self.id[j]][self.id[i]] ** 2 * x_j + x_i) * \
                                 (1 - (x_j + x_i) / (bv[self.id[j]][self.id[i]] ** 2 * x_j + x_i))

        sumai_v = 0
        sumaij_v = 0
        sumai_T = 0
        sumaij_T = 0
        m = 0
        for i, componentei in enumerate(self.comp):
            sumai_v += xy[i] ** 2 / componentei.rhocL
            sumai_T += xy[i] ** 2 * componentei.Tc
            m += xy[i] * componentei.M
            for j, componentej in enumerate(self.comp):
                if j > i:
                    sumaij_v += c_rho[i, j] * f_rho[i, j]
                    sumaij_T += c_T[i, j] * f_T[i, j]
        rhor = 1. / (sumai_v + sumaij_v) * m
        Tr = sumai_T + sumaij_T

        Tcxixy = [None] * len(self.comp)
        rhocxixy = [None] * len(self.comp)
        Tcxi = [None] * len(self.comp)
        rhocxi = [None] * len(self.comp)
        for i, componentei in enumerate(self.comp):
            sumav1 = 0
            sumat1 = 0
            for k in range(i):
                sumav1 += c_rho[k, i] * dFrho_ki[k, i]
                sumat1 += c_T[k, i] * dFT_ki[k, i]
            sumav2 = 0
            sumat2 = 0
            for k in range(i + 1, len(xy)):
                sumav2 += c_rho[i, k] * dFrho_ik[i, k]
                sumat2 += c_T[i, k] * dFT_ik[i, k]
            Tcxi[i] = 2 * xy[i] * componentei.Tc + sumat1 + sumat2
            rhocxi[i] = 2 * xy[i] / componentei.rhocL + sumav1 + sumav2
        for i in range(len(Tcxi)):
            Tcxixy[i] = Tcxi[i]
            rhocxixy[i] = rhocxi[i] / m
        return rhor, Tr, Tcxixy, rhocxixy, m