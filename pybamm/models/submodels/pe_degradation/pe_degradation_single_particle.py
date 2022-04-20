#
# Class for a single particle with cathode degradation.
#
import pybamm

from .base_pe_degradation import BasePeDegradation

class PeDegradationSingleParticle(BasePeDegradation):
    """
    Class for phase transition in a single x-averaged particle of Positive Electrode for degradation.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str
        The domain of the model must be 'Positive'

    **Extends:** :class:`pybamm.phase_transition.BaseTransition`
    """

    def __init__(self, param, domain="Positive"):

        # check whether the domain is correct
        if domain != "Positive":
            raise ValueError(
                "Value of domain must be 'Positive' for phase transition degradation"
            )
        super().__init__(param, domain)

    def get_fundamental_variables(self):

        # concentration in particle core
        c_c_xav = pybamm.Variable(
            "X-averaged positive core concentration",
            domain="positive core",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, 1),
        )
        c_c = pybamm.SecondaryBroadcast(c_c_xav, ["positive electrode"])

        # oxygen concentration in degraded passivation layer/shell
        c_o_xav = pybamm.Variable(
            "X-averaged positive shell concentration of oxygen",
            domain="positive shell oxygen",
            auxiliary_domains={"secondary": "current collector"},
        )
        c_o = pybamm.SecondaryBroadcast(c_o_xav, ["positive electrode"])

        # location of core-shell boundary
        s_xav = pybamm.Variable(
            "X-averaged moving phase boundary location",
            domain="current collector",
            bounds=(0, 1),
        )
        s = pybamm.PrimaryBroadcast(s_xav, ["positive electrode"])

        variables = self._get_standard_concentration_variables(
            c_c=c_c, c_o=c_o, s=s, 
            c_c_xav=c_c_xav, c_o_xav=c_o_xav, s_xav=s_xav
        )

        return variables

    def get_coupled_variables(self, variables):

        T_xav = variables["X-averaged positive electrode temperature"]

        # get s_dot
        c_o_cent_av = variables[
            "X-averaged positive shell center concentration of oxygen"
        ]
        c_c_surf_av = variables["X-averaged positive core surface concentration"]

        # EqualHeaviside evaluates whether left <= right
        s_dot = ( 
            -(self.kappa_1 * self.k_1(T_xav) - self.kappa_2 * self.k_2(T_xav) * c_o_cent_av)
            * pybamm.EqualHeaviside(c_c_surf_av, self.c_p_thrd)
            # * pybamm.sigmoid(c_c_surf_av, self.c_p_thrd, 10)
        ) 
        # s_dot = pybamm.Scalar(-0.07)

        # get c_c and c_o at boundary
        c_c_xav = variables["X-averaged positive core concentration"]
        c_o_xav = variables["X-averaged positive shell concentration of oxygen"]
        s_xav = variables["X-averaged moving phase boundary location"]

        T_xav_c = pybamm.PrimaryBroadcast(T_xav, ["positive core"])
        T_xav_o = pybamm.PrimaryBroadcast(T_xav, ["positive shell oxygen"])

        # defined in base_phase_transition
        D_c = pybamm.surf(self.D_c(c_c_xav, T_xav_c))
        D_o = pybamm.boundary_value(self.D_o(c_o_xav, T_xav_o), "left")

        c_c_N_av = variables["X-averaged positive core surface cell concentration"]
        c_o_1_av = variables[
            "X-averaged positive shell center cell concentration of oxygen"
        ]
        dx_cp_av = variables["X-averaged positive core surface cell length"]
        dx_co_av = variables[
            "X-averaged positive shell center cell length of oxygen"
        ]

        j_xav = variables[
            "X-averaged positive electrode interfacial current density"
        ]

        # the interface boundary value is calculated from applied boundary condition
        # not extrapolated afterwards
        c_c_b_xav = (
            (c_c_N_av - dx_cp_av / D_c * (
              self.C_c / self.param.C_d * s_xav / self.param.a_R_p / self.param.gamma_p * j_xav / (s_xav ** 2)
              - self.C_c * s_xav * s_dot * self.c_s_trap)
            )
            / (1 + dx_cp_av / D_c * self.C_c * s_xav * s_dot)
        )
        c_c_b = pybamm.PrimaryBroadcast(c_c_b_xav, ["positive electrode"])

        # boundary oxygen concentration from applied bc
        # not extrapolated afterwards
        c_o_b_xav = (
            (c_o_1_av - dx_co_av / D_o * self.C_o * s_dot * (1 - s_xav) * self.c_o_core)
            / (1 - dx_co_av / D_o * self.C_o * s_dot * (1 - s_xav))
        )
        c_o_b = pybamm.PrimaryBroadcast(c_o_b_xav, ["positive electrode"])

        variables.update(
            {
                "X-averaged time derivative of moving phase boundary location": s_dot,
                "X-averaged lithium concentration at core-shell interface": c_c_b_xav,
                "X-averaged lithium concentration at core-shell interface [mol.m-3]": 
                    c_c_b_xav * self.param.c_p_max,
                "Lithium concentration at core-shell interface": c_c_b,
                "Lithium concentration at core-shell interface [mol.m-3]": 
                    c_c_b * self.param.c_p_max,
                "X-averaged oxygen concentration at core-shell interface": c_o_b_xav,
                "X-averaged oxygen concentration at core-shell interface [mol.m-3]": 
                    c_o_b_xav * self.param.c_p_max,
                "Oxygen concentration at core-shell interface": c_o_b,
                "Oxygen concentration at core-shell interface [mol.m-3]": 
                    c_o_b * self.param.c_p_max,
            }
        )

        variables.update(self._get_total_concentration_variables(variables))

        return variables

    def set_rhs(self, variables):

        c_c_xav = variables["X-averaged positive core concentration"]
        c_o_xav = variables["X-averaged positive shell concentration of oxygen"]
        s_xav = variables["X-averaged moving phase boundary location"]

        T_xav = variables["X-averaged positive electrode temperature"]
        T_xav_c = pybamm.PrimaryBroadcast(T_xav, ["positive core"])
        T_xav_o = pybamm.PrimaryBroadcast(T_xav, ["positive shell oxygen"])

        # use FunctionParameter "Positive electrode diffusivity [m2.s-1]"
        # in lithium_ion_parameters for the one of core phase
        # defined in base_phase_transition
        D_c = self.D_c(c_c_xav, T_xav_c)
        D_o = self.D_o(c_o_xav, T_xav_o)

        eta = pybamm.pe_degradation.eta_xav
        psi = pybamm.pe_degradation.psi_xav

        s_dot = variables["X-averaged time derivative of moving phase boundary location"]

        # cartesian coordinate
        # self.rhs[c_c_xav] = (
        #     pybamm.inner(eta * s_dot / s_xav, pybamm.grad(c_c_xav)) 
        #     + 1 / self.C_c / ((eta * s_xav) ** 2) 
        #     * pybamm.div(eta ** 2 * D_c * pybamm.grad(c_c_xav))
        # )
        # spherical polar
        self.rhs[c_c_xav] = (
            pybamm.inner(eta * s_dot / s_xav, pybamm.grad(c_c_xav)) 
            + 1 / self.C_c / (s_xav ** 2) 
            * pybamm.div(D_c * pybamm.grad(c_c_xav))
        )
        # test 
        # self.rhs[c_c_xav] = pybamm.div(pybamm.grad(c_c_xav))

        self.rhs[c_o_xav] = (
            pybamm.inner((1 - psi) * s_dot / (1 - s_xav), pybamm.grad(c_o_xav))
            + 1 / self.C_o / ((1 - s_xav) ** 2 * (psi * (1 - s_xav) + s_xav) ** 2)
            * pybamm.div((psi * (1 - s_xav) + s_xav) ** 2 * D_o * pybamm.grad(c_o_xav))
        )
        self.rhs[s_xav] = s_dot 

    def set_boundary_conditions(self, variables):

        c_c_xav = variables["X-averaged positive core concentration"]
        c_o_xav = variables["X-averaged positive shell concentration of oxygen"]
        s_xav = variables["X-averaged moving phase boundary location"]

        T_xav = variables["X-averaged positive electrode temperature"]

        c_c_N_av = variables["X-averaged positive core surface cell concentration"]
        c_o_1_av = variables["X-averaged positive shell center cell concentration of oxygen"]
        dx_cp_av = variables["X-averaged positive core surface cell length"] 
        dx_co_av = variables["X-averaged positive shell center cell length of oxygen"]

        c_c_b = variables["X-averaged lithium concentration at core-shell interface"]
        c_o_b   = variables["X-averaged oxygen concentration at core-shell interface"]

        rbc_cc = (c_c_b - c_c_N_av) / dx_cp_av
        lbc_co = (c_o_1_av - c_o_b) / dx_co_av

        self.boundary_conditions[c_c_xav] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (rbc_cc, "Neumann")
            # "right": (rbc_cs, "Neumann")
        }
        self.boundary_conditions[c_o_xav] = {
            "left":  (lbc_co, "Neumann"),
            "right": (pybamm.Scalar(0), "Dirichlet")
        }

    def set_initial_conditions(self, variables):
        """
        For single particle models, initial conditions can't depend on x so we
        arbitrarily set the initial values of the single particles to be given
        by the values at x=1 in the positive electrode. 
        Typically, supplied initial conditions are uniform x.
        """
        c_c_xav = variables["X-averaged positive core concentration"]
        c_o_xav = variables["X-averaged positive shell concentration of oxygen"]
        s_xav = variables["X-averaged moving phase boundary location"]

        eta = pybamm.pe_degradation.eta_xav
        psi = pybamm.pe_degradation.psi_xav

        self.initial_conditions[c_c_xav] = self.c_c_init(eta)
        self.initial_conditions[c_o_xav] = self.c_o_init(psi)
        self.initial_conditions[s_xav] = self.s_init(1)