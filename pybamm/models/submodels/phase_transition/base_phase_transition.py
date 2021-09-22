#
# Base class for cathode phase transition models.
#
# Developed by Mingzhao Zhuo @Imperial July 2021
#
import numpy as np
import pybamm

# spatial variables for DFN many particle model
eta = pybamm.SpatialVariable(
    "eta", 
    domain=["positive core"],
    auxiliary_domains={
        "secondary": "positive electrode",
        "tertiary": "current collector",
    },
    coord_sys="cartesian"
)
chi = pybamm.SpatialVariable(
    "chi", 
    domain=["positive shell"],
    auxiliary_domains={
        "secondary": "positive electrode",
        "tertiary": "current collector",
    },
    coord_sys="cartesian"
)
psi = pybamm.SpatialVariable(
    "psi", 
    domain=["positive shell oxygen"],
    auxiliary_domains={
        "secondary": "positive electrode",
        "tertiary": "current collector",
    },
    coord_sys="cartesian"
)
# spatial variables for SPM single particle model
eta_xav = pybamm.SpatialVariable(
    "eta", 
    domain=["positive core"],
    auxiliary_domains={
        "secondary": "current collector",
    },
    coord_sys="cartesian"
)
chi_xav = pybamm.SpatialVariable(
    "chi", 
    domain=["positive shell"],
    auxiliary_domains={
        "secondary": "current collector",
    },
    coord_sys="cartesian"
)
psi_xav = pybamm.SpatialVariable(
    "psi", 
    domain=["positive shell oxygen"],
    auxiliary_domains={
        "secondary": "current collector",
    },
    coord_sys="cartesian"
)

class BaseTransition(pybamm.BaseSubModel):
    """
    Base class for positive electrode phase transition models. 

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : dict, optional
        Dictionary of either the electrode for "Positive" or "Nagative"

    References
    ----------
    .. [1] Ghosh, Abir, et al. "A Shrinking-Core Model for the Degradation of
           High-Nickel Cathodes (NMC811) in Li-Ion Batteries: Passivation Layer
           Growth and Oxygen Evolution." Journal of The Electrochemical Society
           168.2 (2021): 020509.

    **Extends:** :class:`pybamm.BaseSubModel`
    """

    def __init__(self, param, domain):
        super().__init__(param, domain)

        pybamm.citations.register("Ghosh2021")

        tau_diffusion_c = self.param.R_p_typ ** 2 / self.D_c_dimensional(
            pybamm.Scalar(1), self.param.T_ref
        )
        tau_diffusion_s = self.param.R_p_typ ** 2 / self.D_s_dimensional(
            pybamm.Scalar(1), self.param.T_ref
        )
        tau_diffusion_o = self.param.R_p_typ ** 2 / self.D_o_dimensional(
            pybamm.Scalar(1), self.param.T_ref
        )

        self.C_c = tau_diffusion_c / self.param.timescale
        self.C_s = tau_diffusion_s / self.param.timescale
        self.C_o = tau_diffusion_o / self.param.timescale

        self.phi = pybamm.Parameter(
            "Ratio of maximum lithium concentration between core and shell")

        c_o_core_dim = pybamm.Parameter(
            "Constant oxygen concentration in particle core [mol.m-3]")
        c_p_thrd_dim = pybamm.Parameter(
            "Threshold concentration for phase transition [mol.m-3]")

        self.c_o_core = c_o_core_dim / self.param.c_p_max
        self.c_p_thrd = c_p_thrd_dim / self.param.c_p_max

        # dimensionless coefficients relating to chemical reactions
        self.kappa_1 = (
            self.k_1_dimensional(self.param.T_ref) 
            * self.param.timescale
            / self.param.R_p_typ
        )
        self.kappa_2 = (
            self.k_2_dimensional(self.param.T_ref) 
            * self.param.timescale
            / self.param.R_p_typ
            * self.param.c_p_max
        )

    def _get_standard_concentration_variables(
        self, c_c, c_s, c_o, s,
        c_c_xav=None, c_s_xav=None, c_o_xav=None, s_xav=None
    ):
        # c_c, c_c_xav=None, # c_c_rav=None, c_c_av=None, c_c_surf=None,
        # c_s, c_s_xav=None, # c_s_rav=None, c_s_av=None, c_s_surf=None,
        # c_o, c_o_xav=None, # c_o_rav=None, c_o_av=None, c_o_surf=None,
        # s, s_xav=None
        """
        All phase transition submodels must provide the core concentration, shell 
        concentration, oxygen concentration in shell, and phase boundary location
        as arguments.
        """

        # Get surface and center concentration if not provided as fundamental
        # variable to solve for
        c_c_surf = pybamm.surf(c_c)
        c_c_surf_av = pybamm.x_average(c_c_surf)

        c_s_surf = pybamm.surf(c_s)
        c_s_surf_av = pybamm.x_average(c_s_surf)
        c_s_cent = pybamm.boundary_value(c_s, "left")
        c_s_cent_av = pybamm.x_average(c_s_cent)

        c_o_surf = pybamm.surf(c_o)
        c_o_surf_av = pybamm.x_average(c_o_surf)
        c_o_cent = pybamm.boundary_value(c_o, "left")
        c_o_cent_av = pybamm.x_average(c_o_cent)

        # maximum concentration as the ref.
        c_scale = self.param.c_p_max
        R_scale = self.param.R_p_typ

        # s = s
        s_xav = s_xav or pybamm.x_average(s)

        # Get average concentration(s) if not provided as fundamental variable to
        # solve for
        c_c_xav = c_c_xav or pybamm.x_average(c_c)
        c_c_rav = self._pe_r_average(c_c, s)
        # c_c_rav = pybamm.r_average(c_c)
        c_c_av = pybamm.x_average(c_c_rav)

        c_s_xav = c_s_xav or pybamm.x_average(c_s)
        c_s_rav = self._pe_r_average(c_s, s)
        c_s_av = pybamm.x_average(c_s_rav)

        c_o_xav = c_o_xav or pybamm.x_average(c_o)
        c_o_rav = self._pe_r_average(c_o, s)
        c_o_av = pybamm.x_average(c_o_rav)

        # boundary cell value
        # surface (rightmost) cell for core c_c and
        c_c_N = pybamm.boundary_cell_value(c_c, "right")
        c_c_N_av = pybamm.x_average(c_c_N)
        # compare with
        # c_c_N_av = pybamm.boundary_cell_value(c_c_xav, "right")

        # center (leftmost) cell for shell c_s
        c_s_1 = pybamm.boundary_cell_value(c_s, "left")
        c_s_1_av = pybamm.x_average(c_s_1)
        # center (leftmost) cell for shell c_o
        c_o_1 = pybamm.boundary_cell_value(c_o, "left")
        c_o_1_av = pybamm.x_average(c_o_1)

        # boundary cell length: node to the edge
        dx_cp = pybamm.boundary_cell_length(c_c, "right")
        dx_cp_av = pybamm.x_average(dx_cp)

        dx_cs = pybamm.boundary_cell_length(c_s, "left")
        dx_cs_av = pybamm.x_average(dx_cs)

        dx_co = pybamm.boundary_cell_length(c_o, "left")
        dx_co_av = pybamm.x_average(dx_co)

        variables = {
            # for core concentration c_c
            "Positive core concentration": c_c,
            "Positive core concentration [mol.m-3]": c_c * c_scale,
            "X-averaged positive core concentration": c_c_xav,
            "X-averaged positive core concentration [mol.m-3]": c_c_xav * c_scale,
            "R-averaged positive core concentration": c_c_rav,
            "R-averaged positive core concentration [mol.m-3]": c_c_rav * c_scale,
            "Average positive core concentration": c_c_av,
            "Average positive core concentration [mol.m-3]": c_c_av * c_scale,
            "Positive core surface concentration": c_c_surf,
            "Positive core surface concentration [mol.m-3]": c_scale * c_c_surf,
            "X-averaged positive core surface concentration": c_c_surf_av,
            "X-averaged positive core surface concentration [mol.m-3]": c_scale * c_c_surf_av,
            "Minimum positive core concentration": pybamm.min(c_c),
            "Maximum positive core concentration": pybamm.max(c_c),
            "Minimum positive core concentration [mol.m-3]": pybamm.min(c_c) * c_scale,
            "Maximum positive core concentration [mol.m-3]": pybamm.max(c_c) * c_scale,
            "Minimum positive core surface concentration": pybamm.min(c_c_surf),
            "Maximum positive core surface concentration": pybamm.max(c_c_surf),
            "Minimum positive core surface concentration [mol.m-3]": pybamm.min(c_c_surf)
            * c_scale,
            "Maximum positive core surface concentration [mol.m-3]": pybamm.max(c_c_surf)
            * c_scale,
            # for shell concentration c_s
            "Positive shell concentration": c_s,
            "Positive shell concentration [mol.m-3]": c_s * c_scale,
            "X-averaged positive shell concentration": c_s_xav,
            "X-averaged positive shell concentration [mol.m-3]": c_s_xav * c_scale,
            "R-averaged positive shell concentration": c_s_rav,
            "R-averaged positive shell concentration [mol.m-3]": c_s_rav * c_scale,
            "Positive shell surface concentration": c_s_surf,
            "Positive shell surface concentration [mol.m-3]": c_scale * c_s_surf,
            "X-averaged positive shell surface concentration": c_s_surf_av,
            "X-averaged positive shell surface concentration [mol.m-3]": c_scale * c_s_surf_av,
            #
            # for shell concentration of oxygen c_o
            "Positive shell concentration of oxygen": c_o,
            "Positive shell concentration of oxygen [mol.m-3]": c_o * c_scale,
            "X-averaged positive shell concentration of oxygen": c_o_xav,
            "X-averaged positive shell concentration of oxygen [mol.m-3]": c_o_xav * c_scale,
            "R-averaged positive shell concentration of oxygen": c_o_rav,
            "R-averaged positive shell concentration of oxygen [mol.m-3]": c_o_rav * c_scale,
            "Positive shell center concentration of oxygen": c_o_cent,
            "Positive shell center concentration of oxygen [mol.m-3]": c_scale * c_o_cent,
            "X-averaged positive shell center concentration of oxygen": c_o_cent_av,
            "X-averaged positive shell center concentration of oxygen [mol.m-3]": c_scale * c_o_cent_av,
            # for moving phase boundary
            "Moving phase boundary location": s,
            "Moving phase boundary location [m]": s * R_scale,
            "X-averaged moving phase boundary location": s_xav,
            "X-averaged moving phase boundary location [m]": s_xav * R_scale,
            # for boundary cell value and length (geometry)
            "Positive core surface cell concentration": c_c_N,
            "Positive core surface cell concentration [mol.m-3]": c_scale * c_c_N,
            "Positive shell center cell concentration": c_s_1,
            "Positive shell center cell concentration [mol.m-3]": c_s_1 * c_scale,
            "Positive shell center cell concentration of oxygen": c_o_1,
            "Positive shell center cell concentration of oxygen [mol.m-3]": c_o_1 * c_scale,
            "X-averaged positive core surface cell concentration": c_c_N_av,
            "X-averaged positive core surface cell concentration [mol.m-3]": c_scale * c_c_N_av,
            "X-averaged positive shell center cell concentration": c_s_1_av,
            "X-averaged positive shell center cell concentration [mol.m-3]": c_s_1_av * c_scale,
            "X-averaged positive shell center cell concentration of oxygen": c_o_1_av,
            "X-averaged positive shell center cell concentration of oxygen [mol.m-3]": c_o_1_av * c_scale,
            # half cell length --- center node to edge
            "Positive core surface cell length": dx_cp,
            "Positive shell center cell length": dx_cs,
            "Positive shell center cell length of oxygen": dx_co,
            "X-averaged positive core surface cell length": dx_cp_av,
            "X-averaged positive shell center cell length": dx_cs_av,
            "X-averaged positive shell center cell length of oxygen": dx_co_av,
            # ---------------------------------------------------------------------
            # we need to supply a variable to
            # submodel 'positive interface' for voltage calculation
            "Positive particle surface concentration": c_s_surf,
            "Positive particle surface concentration [mol.m-3]": c_scale * c_s_surf,
            # ---------------------------------------------------------------------
        }

        return variables

    def _pe_r_average(self, symbol, s):

        if symbol.domain == ["positive core"]:
            eta = pybamm.SpatialVariable("eta", symbol.domain, symbol.auxiliary_domains)
            v = pybamm.FullBroadcast(
                pybamm.Scalar(1), symbol.domain, symbol.auxiliary_domains
            )
            coeff = 4 * np.pi * (eta * s) ** 2 * s
            return pybamm.Integral(coeff * symbol, eta) / pybamm.Integral(coeff * v, eta)
        elif symbol.domain in [["positive shell"], ["positive shell oxygen"]]:
            chi = pybamm.SpatialVariable("chi", symbol.domain, symbol.auxiliary_domains)
            v = pybamm.FullBroadcast(
                pybamm.Scalar(1), symbol.domain, symbol.auxiliary_domains
            )
            coeff = 4 * np.pi * ((1 - s) * chi + s) ** 2 * (1 - s)
            return pybamm.Integral(coeff * symbol, chi) / pybamm.Integral(coeff * v, chi)
        else:
            raise pybamm.DomainError(
                "domain must be positive core or shell (oxygen)."
                # "symbol '{}' domain {} is not 'positive core or shell (oxygen)'.".format(
                #     symbol.name, symbol.domain
                # )
            )

    def _get_total_concentration_variables(self, variables):

        s = variables["Moving phase boundary location"]

        # c_c = variables["Positive core concentration"]
        # c_s = variables["Positive shell concentration"]
        # # calculate r_average
        # c_c_rav = self._pe_r_average(c_c)
        # c_s_rav = self._pe_r_average(c_s, s)

        c_c_rav = variables["R-averaged positive core concentration"]
        c_s_rav = variables["R-averaged positive shell concentration"]

        eps_p = variables["Positive electrode active material volume fraction"]
        eps_p_av = pybamm.x_average(eps_p)

        # the question is to multiply phi or not self.phi * 
        c_c_vol_av = (
            pybamm.x_average(eps_p * (c_c_rav * s ** 3 + self.phi * c_s_rav * (1 - s ** 3))) 
            # pybamm.x_average(eps_p * (c_c_rav)) 
            / eps_p_av
        )
        c_scale = self.param.c_p_max

        # Positive electrode thickness [m]
        L = self.param.L_p
        # Area of current collector
        A = self.param.A_cc

        variables.update(
            {
                # "Positive electrode SOC": c_c_vol_av,
                "Positive electrode volume-averaged concentration": c_c_vol_av,
                "Positive electrode volume-averaged concentration [mol.m-3]": c_c_vol_av * c_scale,
                "Total lithium in positive electrode [mol]": pybamm.yz_average(c_c_vol_av  * eps_p_av)
                    * c_scale 
                    * L 
                    * A,
            }
        )
        return variables

    # def set_events(self, variables):

    #     c_s_surf = variables["Positive shell surface concentration"]
    #     tol = 1e-4

    #     self.events.append(
    #         pybamm.Event(
    #             "Minimum positive shell surface concentration",
    #             pybamm.min(c_s_surf) - tol,
    #             pybamm.EventType.TERMINATION,
    #         )
    #     )

    #     self.events.append(
    #         pybamm.Event(
    #             "Maximum positive shell surface concentration",
    #             (1 - tol) - pybamm.max(c_s_surf),
    #             pybamm.EventType.TERMINATION,
    #         )
    #     )

    # --------------------------------------------------------------------
    # define parameters that are exclusive to PE phase transition model
    # positive core will use already defined positive particle
    # --------------------------------------------------------------------

    def D_c_dimensional(self, sto, T):
        """Dimensional diffusivity in positive core. Note this is defined as a
        function of stochiometry"""
        inputs = {"Positive core stoichiometry": sto, "Temperature [K]": T}
        return pybamm.FunctionParameter("Positive core diffusivity [m2.s-1]", inputs)

    def D_c(self, c_c, T):
        """Dimensionless positive core diffusivity"""
        sto = c_c
        T_dim = self.param.Delta_T * T + self.param.T_ref
        return self.D_c_dimensional(sto, T_dim) / self.D_c_dimensional(
            pybamm.Scalar(1), self.param.T_ref
        )

    def D_s_dimensional(self, sto, T):
        """Dimensional diffusivity in positive shell. Note this is defined as a
        function of stochiometry"""
        inputs = {"Positive shell stoichiometry": sto, "Temperature [K]": T}
        return pybamm.FunctionParameter("Positive shell diffusivity [m2.s-1]", inputs)

    def D_s(self, c_s, T):
        """Dimensionless positive shell diffusivity"""
        sto = c_s
        T_dim = self.param.Delta_T * T + self.param.T_ref
        return self.D_s_dimensional(sto, T_dim) / self.D_s_dimensional(
            pybamm.Scalar(1), self.param.T_ref
        )

    def D_o_dimensional(self, sto, T):
        """Dimensional oxygen diffusivity in positive shell. Note this is defined as a
        function of stochiometry"""
        inputs = {"Positive shell oxygen stoichiometry": sto,
                  "Temperature [K]": T}
        return pybamm.FunctionParameter("Positive shell oxygen diffusivity [m2.s-1]", inputs)

    def D_o(self, c_o, T):
        """Dimensionless oxygen diffusivity in positive shell"""
        sto = c_o
        T_dim = self.param.Delta_T * T + self.param.T_ref
        return self.D_o_dimensional(sto, T_dim) / self.D_o_dimensional(
            pybamm.Scalar(1), self.param.T_ref
        )

    def k_1_dimensional(self, T):
        """Dimensional forward chemical reaction coefficient"""
        inputs = {"Temperature [K]": T}
        return pybamm.FunctionParameter("Forward chemical reaction coefficient [m.s-1]", inputs)

    def k_1(self, T):
        """Dimensionless forward chemical reaction coefficient"""
        T_dim = self.param.Delta_T * T + self.param.T_ref
        return self.k_1_dimensional(T_dim) / self.k_1_dimensional(self.param.T_ref)

    def k_2_dimensional(self, T):
        """Dimensional reverse chemical reaction coefficient"""
        inputs = {"Temperature [K]": T}
        return pybamm.FunctionParameter("Reverse chemical reaction coefficient [m4.mol-1.s-1]", inputs)

    def k_2(self, T):
        """Dimensionless reverse chemical reaction coefficient"""
        T_dim = self.param.Delta_T * T + self.param.T_ref
        return self.k_2_dimensional(T_dim) / self.k_2_dimensional(self.param.T_ref)

    def c_c_init_dimensional(self, x):
        """Initial concentration as a function of dimensionless position x"""
        inputs = {"Dimensionless through-cell position (x_p)": x}
        return pybamm.FunctionParameter(
            "Initial concentration in positive core [mol.m-3]", inputs
        )

    def c_c_init(self, x):
        """
        Dimensionless initial concentration as a function of dimensionless position x
        """
        return self.c_c_init_dimensional(x) / self.param.c_p_max

    def c_s_init_dimensional(self, x):
        """Initial concentration as a function of dimensionless position x"""
        inputs = {"Dimensionless through-cell position (x_p)": x}
        return pybamm.FunctionParameter(
            "Initial concentration in positive shell [mol.m-3]", inputs
        )

    def c_s_init(self, x):
        """
        Dimensionless initial concentration as a function of dimensionless position x
        """
        return self.c_s_init_dimensional(x) / self.param.c_p_max

    def c_o_init_dimensional(self, psi):
        """Initial oxygen concentration as a function of dimensionless position x"""
        inputs = {"Dimensionless transformed position shell oxygen (psi)": psi}
        return pybamm.FunctionParameter(
            "Initial oxygen concentration in positive shell [mol.m-3]", inputs
        )

    def c_o_init(self, psi):
        """
        Dimensionless initial oxygen concentration as a function of dimensionless position x
        """
        return self.c_o_init_dimensional(psi) / self.param.c_p_max

    def s_init_dimensional(self, x):
        """Initial phase boundary location as a function of dimensionless position x"""
        inputs = {"Dimensionless through-cell position (x_p)": x}
        return pybamm.FunctionParameter(
            "Initial phase boundary location [m]", inputs
        )

    def s_init(self, x):
        """
        Dimensionless phase boundary location as a function of dimensionless position x
        """
        return self.s_init_dimensional(x) / self.param.R_p_typ
