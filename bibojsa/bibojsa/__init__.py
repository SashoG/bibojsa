"""bibojsa: computational core for GVM-engineered SPDC in biaxial crystals.

Working title context: "Group-velocity-matched photon-pair generation in the
monoclinic crystal BiB3O6: joint spectral engineering and heralded-photon
purity."

Modules
-------
sellmeier    : Sellmeier dispersion, group index k'(omega), GVD k''(omega).
biaxial      : Fresnel equation of wave normals, eigenindices, walk-off.
phasematch   : collinear SPDC phase matching (Delta k, angle solver, QPM).
deff         : effective nonlinearity, LPR2010 Table 1 (Kleinman, point group 2).
gvm          : group-velocity-matching classification (theta_PM loci).
jsa          : joint spectral amplitude, Schmidt decomposition, purity.
materials    : hard-coded Sellmeier coefficients for BIBO, BBO, KDP, KTP.

References
----------
[LPR2010] V. Petrov, ... "Application of the nonlinear crystal BiB3O6 for
    mid-infrared femtosecond pulse generation" (and related surveys),
    Laser & Photonics Reviews 4, 53 (2010). Table 1 (o/e classification,
    deff), Table 2 (Sellmeier fit), Table 3 (index values), Fig. 10b (BIBO
    GVM survey). Cited throughout as LPR2010.
"""

__version__ = "1.0.0"
