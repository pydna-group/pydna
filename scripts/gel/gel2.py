#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib

import numpy as _np
import matplotlib.ticker as _mtick
from matplotlib import pyplot as _plt
from matplotlib import cm as _cm
from matplotlib.ticker import FixedLocator as _FixedLocator

# from mpldatacursor       import datacursor   as _datacursor

from scipy.interpolate import griddata as _griddata
from scipy.optimize import leastsq as _leastsq
from scipy.optimize import fsolve as _fsolve
from scipy import stats as _stats
from io import BytesIO as _BytesIO
from pydna.dseqrecord import Dseqrecord as _Dseqrecord


horizontal_data_table = """
Table 1.Electrophoresis conditions and best-fit parameters
+------+-----+-------+-------+-----------+-----------+----------+--------+-------+--------+
|  E   | %T  | L_max | L_min | mu_S*10^8 | mu_L*10^8 |   gamma  |  chi^2 | L'max | chi^2' |
|(V/cm)|     | (bp)  | (bp)  |(m^2/V.sec)|(m^2/V.sec)|   (kbp)  |        | (bp)  |        |
+------+-----+-------+-------+-----------+-----------+----------+--------+-------+--------+
| 0.71 | 0.5 | 35000 | 1000  |    2.61   |    0.42   |   29.70  | 0.9999 |       |        |
| 0.71 | 1.0 | 20000 |  300  |    2.54   |  1.64E-07 | 2.81E+07 | 0.9981 |  5000 | 0.9992 |
| 0.71 | 1.3 | 20000 |  200  |    2.49   |  2.60E-09 | 1.06E+09 | 0.9647 |  2600 | 0.9990 |
| 1.00 | 1.0 | 48502 |  400  |    2.37   |    0.23   |   29.99  | 0.9997 |       |        |
| 1.00 | 1.1 | 30000 |  200  |    2.21   |    0.05   |   74.08  | 0.9987 | 25000 | 0.9991 |
| 1.00 | 1.2 | 27500 |  200  |    2.14   |    0.01   |  144.40  | 0.9964 |  1300 | 0.9990 |
| 1.23 | 0.7 | 47500 | 1000  |    2.60   |    0.66   |   19.86  | 0.9997 |       |        |
| 1.23 | 0.8 | 35000 |  200  |    2.55   |    0.23   |   23.84  | 0.9988 | 27500 | 0.9990 |
| 1.23 | 0.9 | 25000 |  200  |    2.44   |    0.08   |   49.58  | 0.9971 |  3000 | 0.9992 |
| 1.31 | 0.5 | 35000 |  400  |    2.54   |    0.39   |   19.69  | 0.9999 |       |        |
| 1.31 | 1.0 | 25000 |  200  |    2.52   |    0.15   |   25.75  | 0.9986 | 17500 | 0.9992 |
| 1.31 | 1.3 | 48502 |  200  |    2.48   |    0.05   |   46.75  | 0.9967 | 14000 | 0.9991 |
| 1.51 | 1.0 | 35000 |  200  |    2.24   |    0.29   |   17.86  | 0.9997 |       |        |
| 1.51 | 1.1 | 32500 |  300  |    2.26   |    0.14   |   23.78  | 0.9980 | 20000 | 0.9989 |
| 1.51 | 1.2 | 48502 |  200  |    2.17   |    0.04   |   64.31  | 0.9974 | 15000 | 0.9991 |
| 1.80 | 0.5 | 47500 |  300  |    2.48   |    0.51   |   17.23  | 0.9998 |       |        |
| 1.80 | 1.3 | 30000 |  200  |    2.39   |    0.22   |   11.93  | 0.9980 | 13000 | 0.9991 |
| 2.00 | 1.0 | 48502 |  200  |    2.22   |    0.23   |   13.96  | 0.9982 | 20000 | 0.9992 |
| 2.00 | 1.1 | 48502 |  200  |    2.24   |    0.16   |   13.06  | 0.9971 | 10000 | 0.9991 |
| 2.00 | 1.2 | 17500 |  200  |    2.15   |    0.07   |   30.99  | 0.9979 | 10000 | 0.9990 |
| 2.22 | 0.5 | 35000 |  300  |    2.68   |    0.65   |   12.90  | 0.9995 |       |        |
| 2.22 | 1.0 | 35000 |  200  |    2.64   |    0.43   |    9.56  | 0.9994 |       |        |
| 2.22 | 1.3 | 32500 |  200  |    2.56   |    0.29   |    8.42  | 0.9981 | 17500 | 0.9991 |
| 2.47 | 0.7 | 32500 |  200  |    2.59   |    0.51   |    9.50  | 0.9997 |       |        |
| 2.47 | 0.8 | 35000 |  500  |    2.64   |    0.47   |    8.77  | 0.9998 |       |        |
| 2.47 | 0.9 | 48502 |  200  |    2.57   |    0.39   |    8.98  | 0.9992 |       |        |
| 2.51 | 0.7 | 35000 |  200  |    2.44   |    0.45   |   14.19  | 0.9999 |       |        |
| 2.51 | 0.8 | 35000 |  200  |    2.44   |    0.39   |   12.41  | 0.9996 |       |        |
| 2.51 | 0.9 | 48502 |  200  |    2.43   |    0.27   |   11.03  | 0.9986 | 15000 | 0.9991 |
| 2.78 | 0.5 | 10000 |  200  |    3.32   |    0.92   |    9.30  | 0.9999 |       |        |
| 2.78 | 0.6 | 10000 |  200  |    3.25   |    0.80   |    8.47  | 1.0000 |       |        |
| 2.78 | 0.7 | 10000 |  200  |    3.18   |    0.66   |    9.80  | 0.9998 |       |        |
| 2.78 | 0.8 | 10000 |  200  |    3.06   |    0.55   |    9.84  | 0.9997 |       |        |
| 2.78 | 0.9 | 10000 |  200  |    2.93   |    0.46   |    9.64  | 0.9997 |       |        |
| 2.78 | 1.0 | 10000 |  200  |    2.93   |    0.39   |   10.40  | 0.9997 |       |        |
| 2.78 | 1.1 | 10000 |  200  |    2.93   |    0.33   |   10.35  | 0.9993 |       |        |
| 2.78 | 1.2 | 10000 |  200  |    2.85   |    0.30   |   10.82  | 0.9995 |       |        |
| 2.78 | 1.3 | 10000 |  200  |    2.84   |    0.25   |   11.14  | 0.9996 |       |        |
| 2.78 | 1.4 | 10000 |  200  |    2.79   |    0.22   |   11.52  | 0.9994 |       |        |
| 2.78 | 1.5 | 10000 |  200  |    2.70   |    0.19   |   11.63  | 0.9993 |       |        |
| 3.10 | 1.0 | 32500 |  200  |    2.65   |    0.58   |    8.61  | 0.9999 |       |        |
| 3.10 | 1.1 | 30000 |  200  |    2.62   |    0.55   |    7.99  | 0.9998 |       |        |
| 3.10 | 1.2 | 27500 |  200  |    2.56   |    0.33   |    6.06  | 0.9980 |  5000 | 0.9994 |
| 3.51 | 0.5 | 35000 |  900  |    2.58   |    0.75   |   13.34  | 0.9993 |       |        |
| 3.51 | 1.0 | 48502 |  200  |    2.58   |    0.52   |    6.82  | 0.9997 |       |        |
| 3.51 | 1.3 | 25000 |  200  |    2.51   |    0.38   |    6.05  | 0.9989 |  8000 | 0.9990 |
| 5.00 | 0.5 | 22500 | 1000  |    2.53   |    0.75   |   18.15  | 0.9996 |       |        |
| 5.00 | 1.0 | 20000 |  400  |    2.55   |    0.67   |    5.78  | 0.9995 |       |        |
| 5.00 | 1.3 | 17500 |  200  |    2.51   |    0.50   |    4.22  | 0.9995 |       |        |
+------+-----+-------+-------+-----------+-----------+----------+--------+-------+--------+
a) The maximum and minimum lengths of DNA fragments that could be enumerated.
b) Best-fit parameters obtained from Eq. (1a) using full range of observed lengths
c) The x2 values for fits to Eq. (1a) using the parameters listed
d) The maximum length DNA fragment which could be included in the fitting procedure to yield
   the neww20values (all$0.999).

Rill, R.L., Beheshti, A., Van Winkle, D.H., 2002. DNA electrophoresis in agarose gels:
effects of field and gel concentration on the exponential dependence of reciprocal mobility
on DNA length. Electrophoresis 23, 2710–2719.
"""

vertical_data_table = """
+------+-----+--------+-------+-----------+-----------+----------+--------+-------+--------+
|  E   | %T  | L_max  | L_min | mu_S*10^8 | mu_L*10^8 |   gamma  |  chi^2 | L'max | chi^2' |
|(V/cm)|     |  (bp)  | (bp)  |(m^2/V.sec)|(m^2/V.sec)|   (kbp)  |        | (bp)  |        |
+------+-----+--------+-------+-----------+-----------+----------+--------+-------+--------+
| 0.62 | 0.4 |  48502 | 1000  |    2.67   |    0.21   |  52.799  | 0.9993 |       |        |
| 0.93 | 0.4 | 194000 | 1000  |    2.54   |    0.34   |  25.304  | 0.9992 |       |        |
| 1.24 | 0.4 |  48502 |  200  |    3.67   |    0.81   |  14.922  | 0.9997 |       |        |
| 1.55 | 0.4 | 194000 |  200  |    3.18   |    0.77   |  14.488  | 0.9992 |       |        |
| 1.86 | 0.4 | 194000 |  200  |    3.12   |    0.78   |  12.603  | 0.9990 |       |        |
| 2.17 | 0.4 | 194000 |  200  |    3.07   |    0.89   |  11.911  | 0.9971 | 17500 | 0.9992 |
| 2.48 | 0.4 | 194000 |  200  |    3.13   |    0.90   |  12.648  | 0.9979 | 15000 | 0.9991 |
| 3.10 | 0.4 | 194000 |  200  |    3.14   |    0.98   |   9.423  | 0.9964 | 12000 | 0.9991 |
| 4.97 | 0.4 | 194000 |  200  |    3.34   |    1.22   |   9.630  | 0.9983 |  8000 | 0.9990 |
| 6.21 | 0.4 | 194000 |  200  |    3.11   |    1.20   |   9.214  | 0.9959 |  4400 | 0.9993 |
| 0.62 | 0.7 |  48502 |  200  |    2.86   |  8.36E-06 | 7.32E+05 | 0.9973 |  2700 | 0.9990 |
| 0.93 | 0.7 | 194000 |  200  |    2.38   |    0.09   |  54.264  | 0.9968 | 25000 | 0.9990 |
| 1.24 | 0.7 |  48502 |  200  |    3.71   |    0.49   |  15.255  | 0.9993 |       |        |
| 1.55 | 0.7 | 194000 |  200  |    3.61   |    0.53   |  13.125  | 0.9994 |       |        |
| 1.86 | 0.7 |  47500 |  300  |    3.61   |    0.62   |   9.834  | 0.9997 |       |        |
| 2.17 | 0.7 | 194000 |  200  |    3.07   |    0.55   |  10.151  | 0.9998 |       |        |
| 2.48 | 0.7 | 194000 |  200  |    3.42   |    0.72   |   9.036  | 0.9997 |       |        |
| 3.10 | 0.7 | 194000 |  200  |    3.14   |    0.69   |   7.713  | 0.9996 |       |        |
| 4.97 | 0.7 | 194000 |  200  |    3.08   |    0.93   |   6.188  | 0.9984 | 15000 | 0.9990 |
| 6.21 | 0.7 | 194000 |  200  |    3.24   |    1.03   |   6.371  | 0.9987 | 10000 | 0.9991 |
| 0.62 | 1.0 |  48502 |  200  |    2.48   |  1.82E-14 | 1.83E+14 | 0.9938 |  2900 | 0.9990 |
| 0.93 | 1.0 | 194000 |  300  |    2.41   |  1.37E-13 | 2.33E+13 | 0.9947 |  2900 | 0.9990 |
| 1.24 | 1.0 |  48502 |  200  |    3.34   |    0.17   |  25.855  | 0.9974 | 14000 | 0.9990 |
| 1.55 | 1.0 |  47500 |  200  |    3.25   |    0.33   |  12.773  | 0.9977 | 14000 | 0.9990 |
| 1.86 | 1.0 |  47500 |  400  |    3.69   |    0.47   |   8.502  | 0.9984 | 10000 | 0.9991 |
| 2.17 | 1.0 | 194000 |  200  |    2.87   |    0.38   |   9.372  | 0.9989 | 15000 | 0.9990 |
| 2.48 | 1.0 | 194000 |  200  |    3.39   |    0.52   |   7.209  | 0.9990 |       |        |
| 3.10 | 1.0 |  19400 |  300  |    3.36   |    0.61   |   5.852  | 0.9994 |       |        |
| 4.97 | 1.0 |  48502 |  200  |    3.04   |    0.77   |   4.749  | 0.9989 | 48502 | 0.9991 |
| 6.21 | 1.0 | 194000 |  200  |    3.13   |    0.88   |   4.417  | 0.9975 | 12500 | 0.9991 |
| 0.62 | 1.3 |  48502 |  200  |    2.12   |  1.87E-10 | 1.11E+10 | 0.9904 |  1800 | 0.9991 |
| 0.93 | 1.3 | 194000 |  300  |    1.96   |  3.80E-06 | 6.10E+05 | 0.9928 |  2600 | 0.9992 |
| 1.24 | 1.3 |  48502 |  300  |    3.32   |    0.02   | 179.189  | 0.9967 |  2600 | 0.9991 |
| 1.55 | 1.3 | 194000 |  200  |    3.47   |    0.27   |  11.920  | 0.9958 |  9000 | 0.9990 |
| 1.86 | 1.3 |  47500 |  400  |    3.97   |    0.36   |   8.594  | 0.9970 |  7000 | 0.9991 |
| 2.17 | 1.3 | 194000 |  200  |    2.76   |    0.21   |  12.412  | 0.9966 | 10000 | 0.9990 |
| 2.48 | 1.3 | 194000 |  200  |    3.26   |    0.45   |   7.464  | 0.9985 | 11000 | 0.9990 |
| 3.10 | 1.3 |  47500 |  200  |    3.37   |    0.50   |   5.159  | 0.9984 |  6000 | 0.9992 |
| 4.97 | 1.3 |  48502 |  200  |    2.68   |    0.49   |   4.048  | 0.9983 |  4000 | 0.9993 |
| 6.21 | 1.3 | 194000 |  200  |    2.70   |    0.59   |   3.555  | 0.9991 |       |        |
+------+-----+--------+-------+-----------+-----------+----------+--------+-------+--------+
"""

# Load data into numpy arrays

dset = _np.genfromtxt(
    _BytesIO(vertical_data_table.encode()),
    delimiter="|",
    dtype=None,
    skip_header=5,
    skip_footer=1,
    usecols=(1, 2, 5, 6, 7),
    names=("E", "T", "muS", "muL", "gamma"),
)


def size_to_mobility(dna_fragment_length, efield, gelconcentration):
    muS = dset["muS"] / (100 * 100)  # m**2 to cm**2/V/s
    muL = dset["muL"] / (100 * 100)  # m**2 to cm**2/V/s
    gamma = dset["gamma"] * 1000  # kbp  to bp
    alpha = 1 / muL - 1 / muS
    beta = 1 / muS

    mobility = _griddata(
        (dset["E"], dset["T"]),
        1 / (beta + alpha * (1 - _np.exp(-dna_fragment_length / gamma))),
        (efield, gelconcentration),
    )
    return mobility.item()


# Constants
kB = 1.3806488e-23  # m**2 * kg / s**2 / K Boltzmann constant
lp = 50  # persistence length of dsDNA (nm)
l = 2 * lp  # Kuhn length (nm)
b = 0.34  # nm/bp length dsDNA (nm/bp)
e = 1.602176487e-19  # elementary charge (1.6-19 A.s)
qeff = 2.2888235528571428e-20  # effective charge per dsDNA base pair (A.s/bp)

# Intrinsic band broadening as function of the diffusion coefficient and time
radius_gyration = lambda L, lp: (lp * L / 3 * (1 - lp / L + lp / L * _np.exp(-L / lp))) ** (1 / 2)

# Relations between basepairs (Nbp), Kuhn segments (NKuhn) and blobs (N)
Nbp_to_N = lambda Nbp, a, b, l: Nbp * (b / l) * (l / a) ** 2  # number of occupied pores

free_solution = lambda kB, T, eta, Rh: kB * T / (6 * _np.pi * eta * Rh)
Zimm_g = lambda Nbp, DRouse, qeff, mu0, kB, T: DRouse * Nbp * qeff / (mu0 * kB * T)
Ogston_Zimm = lambda D0, g: D0 * g
Ogston_Rouse = lambda Nbp, kB, T, a, eta, b, l: kB * T * a**3 / (eta * b**2 * l**2 * Nbp**2)

# Gaussian function
Gaussian = lambda x, hgt, ctr, dev: hgt * _np.exp(-((x - ctr) ** 2) / (2 * dev**2))

Gauss_dev = lambda FWHM: FWHM / (2 * _np.sqrt(2 * _np.log(2)))
Gauss_FWHM = lambda FWTM: FWTM * _np.sqrt(2 * _np.log(2)) / _np.sqrt(2 * _np.log(10))


class Gel:
    def __init__(
        self,
        samples,  # list of lists of linear Dseqrecord objects
        gel_concentration=1,  # % (w/v)
        gel_length=8,  # cm
        wellx=0.7,  # cm
        welly=0.2,  # cm
        wellsep=0.2,
    ):  # cm
        self.samples = samples
        self.gel_concentration = gel_concentration
        self.gel_length = gel_length
        self.wellx = wellx
        self.welly = welly
        self.wellsep = wellsep

    def run(
        self,
        field=5.0,  # V/cm
        temperature=295.15,  # K
        runtime=None,  # seconds
        exposure=0.5,  # [0-1]
        plot=True,
        res=200,  # px/cm
        cursor_ovr={"hover": False},
        back_col=0.3,
        band_col=1,
        well_col=0.05,
        noise=0.015,
        detectlim=0.04,
        interpol="linear",  # 'cubic','nearest'
        dset_name="vertical",  # 'horizontal'
        replNANs=True,
    ):  # replace NANs by 'nearest' interpolation
        max_mob = 0

        for i, lane in enumerate(self.samples):
            for j, frag in enumerate(lane):
                mob = size_to_mobility(len(frag), field, self.gel_concentration)  # cm/s
                frag.mobility = mob
                self.samples[i][j] = frag
                max_mob = max((max_mob, mob))

        # vWBR eq. parameters muL, muS, gamma

        mu = _np.zeros(100)
        for i, Li in enumerate(_np.linspace(100, 50000, 100)):
            mu[i] = size_to_mobility(Li, field, self.gel_concentration)

        muS0 = 3.5e-4  # cm^2/(V.sec)  ############################################
        muL0 = 1.0e-4  # cm^2/(V.sec)  ############################################
        gamma0 = 8000  # bp            ############################################

        vWBR = lambda L, muS, muL, gamma: (1 / muS + (1 / muL - 1 / muS) * (1 - _np.exp(-L / gamma))) ** -1

        def residuals(pars, L, mu):
            return mu - vWBR(L, *pars)

        pars, cov, infodict, mesg, ier = _leastsq(
            residuals,
            [muS0, muL0, gamma0],
            args=(_np.linspace(100, 50000, 100), mu),
            full_output=True,
        )
        muS, muL, gamma = pars

        time = runtime or (0.9 * self.gel_length) / max_mob  # sec

        # Free solution mobility estimate

        space = _np.logspace(_np.log10(100), _np.log10(3000), 10)
        DNAspace_for_mu0 = _np.array([round(val, 0) for val in space])

        Tvals = _np.unique(dset["T"])  # (g/(100 mL))*100

        # Mobility dependence on size (mu(L)) for each agarose percentage (Ti)
        ln_mu_LxT = []
        for Lj in DNAspace_for_mu0:
            ln_mu_T = []
            for Ti in Tvals:
                mu = size_to_mobility(Lj, field, Ti)
                ln_mu_T.append(_np.log(mu))
            ln_mu_LxT.append(ln_mu_T)
        ln_mu_LxT = _np.array(ln_mu_LxT)
        # Linear regression for each DNA size
        lregr_stats = []
        exclude = []
        for l in range(len(DNAspace_for_mu0)):
            not_nan = _np.logical_not(_np.isnan(ln_mu_LxT[l]))
            if sum(not_nan) > 1:
                # (enough points for linear regression)
                gradient, intercept, r_value, p_value, std_err = _stats.linregress(
                    Tvals[not_nan], ln_mu_LxT[l][not_nan]
                )
                lregr_stats.append((gradient, intercept, r_value, p_value, std_err))
                exclude.append(False)
            else:
                exclude.append(True)
        exclude = _np.array(exclude)
        ln_mu_LxT = ln_mu_LxT[~exclude]
        if len(lregr_stats) > 0:
            # Free solution mobility determination
            ln_mu0 = _np.mean([row[1] for row in lregr_stats])  # mean of intercepts
            mu0 = _np.exp(ln_mu0)  # cm^2/(V.seg)
        else:
            mu0 = None

        self.freesol_mob = mu0

        eta = 2.414e-5 * 10 ** (
            247.8 * (temperature - 140)
        )  # (Pa.s)=(kg/(m.s)) accurate to within 2.5% from 0 °C to 370 °C

        pore_size = lambda gamma, muL, mu0, lp, b: (gamma * muL * lp * b / mu0) ** (1 / 2)

        pore_size_fit = lambda C: 143 * C ** (-0.59)

        a = pore_size(gamma, muL, mu0, lp, b)
        a_fit = pore_size_fit(self.gel_concentration)  # ##################################
        a_fit = a_fit.to("m")  # ##################################
        self.poresize = a
        self.poresize_fit = a_fit
        reduced_field = lambda eta, a, mu0, E, kB, T: eta * a**2 * mu0 * E / (kB * T)
        epsilon = reduced_field(eta, a, mu0, field * 100, kB, temperature)
        # Diffusion coefficient of a blob
        Dblob = lambda kB, T, eta, a: kB * T / (eta * a)
        Db = Dblob(kB, temperature, eta, a)

        # Diffusion regime frontiers (in number of occupied pores)

        def diff_Zimm_Rouse(Nbp, args):
            kB, T, qeff, eta, mu0, a, b, l, lp = args
            Nbp = Nbp[0]
            L = Nbp * b
            Rg = radius_gyration(L, lp)
            D0 = free_solution(kB, T, eta, Rg)
            DRouse = Ogston_Rouse(Nbp, kB, T, a, eta, b, l)
            g = Zimm_g(Nbp, DRouse, qeff, mu0, kB, T)
            g = g.to_base_units()
            DZimm = Ogston_Zimm(D0, g)
            return DZimm - DRouse

        Zimm_Rouse = lambda x0, args: (Nbp_to_N(_fsolve(diff_Zimm_Rouse, x0, args)[0], args[5], args[6], args[7]))

        equil_accel = lambda epsilon: epsilon ** (-2 / 3)
        accel_plateau = lambda epsilon: epsilon ** (-1)

        N_lim1 = accel_plateau(epsilon)  # #################################
        N_lim2 = equil_accel(epsilon)  # # ***   Major problem    ***   ##
        N_lim3 = Zimm_Rouse(
            2000,  # #################################
            [kB, temperature, qeff, eta, mu0, a, b, l, lp],
        )

        N_to_Nbp = lambda N, a, b, l: N * (l / b) * (a / l) ** 2  # number of base pairs (bp)
        self.accel_to_plateau = N_to_Nbp(N_lim1, a, b, l)
        self.equil_to_accel = N_to_Nbp(N_lim2, a, b, l)
        self.Zimm_to_Rouse = N_to_Nbp(N_lim3, a, b, l)

        for i, lane in enumerate(self.samples):
            for j, frag in enumerate(lane):
                Nbp = len(frag)
                N = Nbp_to_N(Nbp, a, b, l)
                if N < N_lim3:  # Ogston-Zimm
                    L = Nbp * b  # (m)
                    Rg = radius_gyration(L, lp)  # (m)
                    D0 = free_solution(kB, temperature, eta, Rg)  # (m^2/s)
                    DRouse = Ogston_Rouse(Nbp, kB, temperature, a, eta, b, l)  # (m^2/s)
                    g = Zimm_g(Nbp, DRouse, qeff, mu0, kB, temperature)  # base
                    D = Ogston_Zimm(D0, g)  # unit
                elif N < N_lim2:  # Rouse/Reptation-equilibrium
                    D = Db / N**2
                elif N > N_lim1:  # Reptation-plateau (reptation with orientation)
                    D = Db * epsilon ** (3 / 2)
                else:  # Accelerated-reptation
                    D = Db * epsilon * N ** (-1 / 2)
                frag.band_width = _np.sqrt(2 * D * time) + self.welly
                self.samples[i][j] = frag

        # Total bandwidths
        time0 = self.welly / (mu0 * field)

        bandwidths0 = [mobs * time0[l] * field for l, mobs in enumerate(mobilities)]

        # Max intensities
        raw_Is = []
        maxI = Q_(-_np.inf, "ng/cm")
        minI = Q_(_np.inf, "ng/cm")
        for i, lane in enumerate(self.samples):
            lane_I = []
            for j, frag in enumerate(lane):
                frag_Qty = quantities[i][j]
                frag_Wth = bandwidths[i][j]  # w=FWHM or w=FWTM ???
                if False:
                    FWHM = Gauss_FWHM(frag_Wth)
                else:
                    FWHM = frag_Wth
                std_dev = Gauss_dev(FWHM)
                auc = frag_Qty  # area under curve proportional to DNA quantity
                Gauss_hgt = lambda auc, dev: auc / (dev * _np.sqrt(2 * _np.pi))
                frag_I = Gauss_hgt(auc, std_dev)  # peak height
                if frag_I > maxI:
                    maxI = frag_I
                if frag_I < minI:
                    minI = frag_I
                lane_I.append(frag_I)
            raw_Is.append(lane_I)

        # max intensity normalization
        satI = maxI + exposure * (minI - maxI)
        intensities = [[(1 - back_col) / satI * fragI for fragI in lane] for lane in raw_Is]
        self.intensities = intensities

        # Plot gel
        if plot:
            # Title
            mins, secs = divmod(time.to("s").magnitude, 60)  # time is in secs
            hours, mins = divmod(mins, 60)
            hours = Q_(hours, "hr")
            mins = Q_(mins, "min")
            secs = Q_(secs, "s")
            title = (
                "E = %.2f V/cm\n"
                "C = %.1f %%\n"
                "T = %.2f K\n"
                "t = %d h %02d m\n"
                "expo = %.1f" % (field, self.gel_concentration, temperature, hours, mins, exposure)
            )

            gel_width = len(samples) * (self.wellx + self.wellsep) + self.wellsep  # cm

            pxl_x = int(gel_width * res)
            pxl_y = int(gel_len * res)
            lane_centers = [(l + 1) * self.wellsep + sum(wellx[:l]) + 0.5 * wellx[l] for l in range(nlanes)]
            rgb_arr = _np.zeros(shape=(pxl_y, pxl_x, 3), dtype=_np.float32)
            bandlengths = wellx
            bands_pxlXYmid = []
            # Paint the bands
            for i in range(nlanes):
                distXmid = lane_centers[i]
                pxlXmid = int(round((distXmid * res).magnitude))
                bandlength = bandlengths[i]
                from_x = int(round(((distXmid - bandlength / 2.0) * res).magnitude))
                to_x = int(round(((distXmid + bandlength / 2.0) * res).magnitude))
                bands_pxlXYmid.append([])
                for j in range(len(lanes[i])):
                    distYmid = distances[i][j]
                    pxlYmid = int(round((distYmid * res).magnitude))
                    bands_pxlXYmid[i].append((pxlXmid, pxlYmid))
                    bandwidth = bandwidths[i][j]  # w=FWHM or w=FWTM ???
                    if False:
                        FWHM = Gauss_FWHM(bandwidth)
                    else:
                        FWHM = bandwidth
                    std_dev = Gauss_dev(FWHM)
                    maxI = intensities[i][j]
                    midI = Gaussian(distYmid, maxI, distYmid, std_dev)
                    if pxlYmid < len(rgb_arr):  # band within gel frontiers
                        rgb_arr[pxlYmid, from_x:to_x] += midI
                    bckwdYstop = False if pxlYmid > 0 else True
                    forwdYstop = False if pxlYmid < len(rgb_arr) - 1 else True
                    pxlYbck = pxlYmid - 1
                    pxlYfor = pxlYmid + 1
                    while not bckwdYstop or not forwdYstop:
                        if not bckwdYstop:
                            distYbck = Q_(pxlYbck, "px") / res
                            bckYI = Gaussian(distYbck, maxI, distYmid, std_dev)
                            if pxlYbck < len(rgb_arr):
                                rgb_arr[pxlYbck, from_x:to_x] += bckYI
                            pxlYbck -= 1
                            if bckYI <= 1e-5 or pxlYbck == -1:
                                bckwdYstop = True
                        if not forwdYstop:
                            distYfor = Q_(pxlYfor, "px") / res
                            forYI = Gaussian(distYfor, maxI, distYmid, std_dev)
                            rgb_arr[pxlYfor, from_x:to_x] += forYI
                            pxlYfor += 1
                            if forYI <= 1e-5 or pxlYfor == pxl_y:
                                forwdYstop = True

            # Background color
            if noise is None or noise <= 0:
                rgb_arr += back_col
            else:
                bckg = _np.random.normal(back_col, noise, (len(rgb_arr), len(rgb_arr[0])))
                rgb_arr += bckg[:, :, _np.newaxis]
            # Saturation
            rgb_arr[rgb_arr > 1] = 1
            rgb_arr[rgb_arr < 0] = 0
            # bands_arr = _np.ma.masked_where(rgb_arr == back_col, rgb_arr)  ###########
            bands_arr = rgb_arr

            # Plot

            fig = _plt.figure()

            ax1 = fig.add_subplot(111, facecolor=str(back_col))
            ax1.xaxis.tick_top()
            ax1.yaxis.set_ticks_position("left")
            ax1.spines["left"].set_position(("outward", 8))
            ax1.spines["left"].set_bounds(0, gel_len)
            ax1.spines["right"].set_visible(False)
            ax1.spines["bottom"].set_visible(False)
            ax1.spines["top"].set_visible(False)
            ax1.spines["right"].set_color(str(back_col))
            ax1.spines["bottom"].set_color(str(back_col))
            ax1.xaxis.set_label_position("top")

            # _plt.xticks(lane_centers, names)
            majorLocator = _FixedLocator(list(range(int(gel_len + 1))))
            minorLocator = _FixedLocator(
                [j / 10.0 for k in range(0, int(gel_len + 1) * 10, 10) for j in range(1 + k, 10 + k, 1)]
            )
            ax1.yaxis.set_major_locator(majorLocator)
            ax1.yaxis.set_minor_locator(minorLocator)
            ax1.tick_params(axis="x", which="both", top="off")

            # Gel image
            bands_plt = ax1.imshow(bands_arr, extent=[0, gel_width, gel_len, 0], interpolation="none")
            # Draw wells

            for i in range(nlanes):
                ctr = lane_centers[i]
                wx = wellx[i]
                wy = welly[i]
                ax1.fill_between(
                    x=[ctr - wx / 2, ctr + wx / 2],
                    y1=[0, 0],
                    y2=[-wy, -wy],
                    color=str(well_col),
                )
            # Invisible rectangles overlapping the bands for datacursor to detect
            bands = []
            for i in range(nlanes):
                bandlength = bandlengths[i]
                center = lane_centers[i]
                x = center - bandlength / 2.0
                for j in range(len(lanes[i])):
                    dna_frag = lanes[i][j]
                    bandwidth = bandwidths[i][j]
                    dist = distances[i][j].magnitude
                    y = dist - bandwidth / 2.0
                    pxlX, pxlY = bands_pxlXYmid[i][j]
                    band_midI = bands_arr[pxlY, pxlX][0]
                    alpha = 0 if abs(band_midI - back_col) >= detectlim else 0.4
                    band = _plt.Rectangle(
                        (x, y),
                        bandlength,
                        bandwidth,
                        fc="none",
                        ec="w",
                        ls=":",
                        alpha=alpha,
                        label="{} bp".format(len(dna_frag)),
                    )
                    _plt.gca().add_patch(band)
                    bands.append(band)
            _plt.ylim(gel_len, -max(welly))
            xlim = sum(wellx) + (nlanes + 1) * self.wellsep
            _plt.xlim(0, xlim)
            _plt.ylabel("Distance (cm)")
            _plt.xlabel("Lanes")
            bbox_args = dict(boxstyle="round,pad=0.6", fc="none")
            an1 = _plt.annotate(
                title,
                xy=(0, 0),
                xytext=(xlim + 0.4, (gel_len + max(welly)) / 2.0),
                va="center",
                bbox=bbox_args,
            )
            an1.draggable()
            _plt.gca().set_aspect("equal", adjustable="box")
            cursor_args = dict(
                display="multiple",
                draggable=True,
                hover=False,
                bbox=dict(fc="white"),
                arrowprops=dict(arrowstyle="simple", fc="white", alpha=0.5),
                xytext=(15, -15),
                formatter="{label}".format,
            )
            if cursor_ovr:
                for key in cursor_ovr:
                    cursor_args[key] = cursor_ovr[key]
            if cursor_args["hover"] == True:
                cursor_args["display"] = "single"
            # _datacursor(bands, **cursor_args)
            return _plt
        return None


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
