from pybamm import exp, constants


def initial_oxygen_concentration(x):
    """
    NMC811 diffusivity as a function of stochiometry.

    References
    ----------

    Parameters
    ----------
    x: :class:`pybamm.Symbol`
        Dimensionless position

    Returns
    -------
    :class:`pybamm.Symbol`
        initial oxygen concentration
    """

    # c_o_ini_ref = 15219.321 
    c_o_ini_ref = 0 

    return c_o_ini_ref * (1 - x ** 2)