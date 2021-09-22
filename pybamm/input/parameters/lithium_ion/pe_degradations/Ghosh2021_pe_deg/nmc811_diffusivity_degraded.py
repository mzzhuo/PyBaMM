from pybamm import exp, constants


def nmc811_diffusivity_degraded(sto, T):
    """
    NMC811 diffusivity as a function of stochiometry.

    References
    ----------

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """

    D_ref = 1.0E-15

    aEne = 9977.16

    arrhenius = exp(aEne / constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius
