#
# Class for a single particle with cathode phase transition.
#
import pybamm

from .base_phase_transition import BaseTransition

class PhaseTransitionSingleParticle(BaseTransition):
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

        # concentration in degraded passivation layer/shell
        c_s_xav = pybamm.Variable(
            "X-averaged positive shell concentration",
            domain="positive shell",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, 1),
        )
        c_s = pybamm.SecondaryBroadcast(c_s_xav, ["positive electrode"])

        # oxygen concentration in degraded passivation layer/shell
        c_o_xav = pybamm.Variable(
            "X-averaged positive shell concentration of oxygen",
            domain="positive shell oxygen",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, 1),
        )
        c_o = pybamm.SecondaryBroadcast(c_o_xav, ["positive electrode"])

        # location of two-phase boundary
        s_xav = pybamm.Variable(
            "X-averaged moving phase boundary location",
            domain="current collector",
            bounds=(0, 1),
        )
        s = pybamm.PrimaryBroadcast(s_xav, ["positive electrode"])

        variables = self._get_standard_concentration_variables(
            c_c=c_c, c_s=c_s, c_o=c_o, s=s, 
            c_c_xav=c_c_xav, c_s_xav=c_s_xav, c_o_xav=c_o_xav, s_xav=s_xav
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
        # s_dot = ( 
        #     -(self.kappa_1 * self.k_1(T_xav) - self.kappa_2 * self.k_2(T_xav) * c_o_cent_av)
        #     * pybamm.EqualHeaviside(c_c_surf_av, self.c_p_thrd)
        # ) 
        s_dot = pybamm.Scalar(0)

        # get c_share and c_o at boundary
        c_c_xav = variables["X-averaged positive core concentration"]
        c_s_xav = variables["X-averaged positive shell concentration"]
        c_o_xav = variables["X-averaged positive shell concentration of oxygen"]
        s_xav = variables["X-averaged moving phase boundary location"]

        T_xav_c = pybamm.PrimaryBroadcast(T_xav, ["positive core"])
        T_xav_s = pybamm.PrimaryBroadcast(T_xav, ["positive shell"])
        T_xav_o = pybamm.PrimaryBroadcast(T_xav, ["positive shell oxygen"])

        # defined in base_phase_transition
        D_c = pybamm.surf(self.D_c(c_c_xav, T_xav_c))
        D_s = pybamm.boundary_value(self.D_s(c_s_xav, T_xav_s), "left")
        D_o = pybamm.boundary_value(self.D_o(c_o_xav, T_xav_o), "left")

        c_c_N_av = variables["X-averaged positive core surface cell concentration"]
        c_s_1_av = variables["X-averaged positive shell center cell concentration"]
        c_o_1_av = variables[
            "X-averaged positive shell center cell concentration of oxygen"
        ]
        dx_cp_av = variables["X-averaged positive core surface cell length"]
        dx_cs_av = variables["X-averaged positive shell center cell length"]
        dx_co_av = variables[
            "X-averaged positive shell center cell length of oxygen"
        ]

        # the interface boundary value is calculated from applied boundary condition
        # not extrapolated afterwards
        c_share_xav = (
            (D_c * (1 - s_xav) / self.C_c / dx_cp_av * c_c_N_av 
             + D_s * s_xav / self.C_s / dx_cs_av * c_s_1_av
             + s_dot * s_xav * (1 - s_xav) * self.c_o_core)
            / (D_c * (1 - s_xav) / self.C_c / dx_cp_av 
               + D_s * s_xav / self.C_s / dx_cs_av
               + s_dot * s_xav * (1 - s_xav) * (1 - self.phi))
        )
        # no lithium eaten during phase transition
        # c_share_xav = (
        #     (D_c * (1 - s_xav) / self.C_c / dx_cp_av * c_c_N_av 
        #      + D_s * s_xav / self.C_s / dx_cs_av * c_s_1_av)
        #     / (D_c * (1 - s_xav) / self.C_c / dx_cp_av 
        #        + D_s * s_xav / self.C_s / dx_cs_av)
        # )
        # NO C_oc chemical reaction effect
        # c_share_xav = (
        #     (D_c * (1 - s_xav) / self.C_c / dx_cp_av * c_c_N_av 
        #      + D_s * s_xav / self.C_s / dx_cs_av * c_s_1_av)
        #     / (D_c * (1 - s_xav) / self.C_c / dx_cp_av 
        #        + D_s * s_xav / self.C_s / dx_cs_av
        #        + s_dot * s_xav * (1 - s_xav) * (1 - self.phi))
        # )
        c_share = pybamm.PrimaryBroadcast(c_share_xav, ["positive electrode"])

        # boundary oxygen concentration from applied bc
        # not extrapolated afterwards
        c_o_b_xav = (
            ((1 - s_xav) * s_dot * dx_co_av * self.c_o_core - D_o / self.C_o * c_o_1_av)
            / ((1 - s_xav) * s_dot * dx_co_av - D_o / self.C_o)
        )
        c_o_b = pybamm.PrimaryBroadcast(c_o_b_xav, ["positive electrode"])

        variables.update(
            {
                "X-averaged time derivative of moving phase boundary location": s_dot,
                "X-averaged shared concentration at core-shell interface": c_share_xav,
                "X-averaged shared concentration at core-shell interface [mol.m-3]": 
                    c_share_xav * self.param.c_p_max,
                "Shared concentration at core-shell interface": c_share,
                "Shared concentration at core-shell interface [mol.m-3]": 
                    c_share * self.param.c_p_max,
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
        c_s_xav = variables["X-averaged positive shell concentration"]
        c_o_xav = variables["X-averaged positive shell concentration of oxygen"]
        s_xav = variables["X-averaged moving phase boundary location"]

        T_xav = variables["X-averaged positive electrode temperature"]
        T_xav_c = pybamm.PrimaryBroadcast(T_xav, ["positive core"])
        T_xav_s = pybamm.PrimaryBroadcast(T_xav, ["positive shell"])
        T_xav_o = pybamm.PrimaryBroadcast(T_xav, ["positive shell oxygen"])

        # use FunctionParameter "Positive electrode diffusivity [m2.s-1]"
        # in lithium_ion_parameters for the one of core phase
        # defined in base_phase_transition
        D_c = self.D_c(c_c_xav, T_xav_c)
        D_s = self.D_s(c_s_xav, T_xav_s)
        D_o = self.D_o(c_o_xav, T_xav_o)

        eta = pybamm.phase_transition.eta_xav
        chi = pybamm.phase_transition.chi_xav
        psi = pybamm.phase_transition.psi_xav

        s_dot = variables["X-averaged time derivative of moving phase boundary location"]

        self.rhs[c_c_xav] = (
            pybamm.inner(eta * s_dot / s_xav, pybamm.grad(c_c_xav)) 
            + 1 / self.C_c / ((eta * s_xav) ** 2) 
            * pybamm.div(eta ** 2 * D_c * pybamm.grad(c_c_xav))
        )
        # self.rhs[c_c_xav] = (
        #     1 / self.C_c / (eta ** 2) 
        #     * pybamm.div(eta ** 2 * D_c * pybamm.grad(c_c_xav))
        # )

        self.rhs[c_s_xav] = (
            pybamm.inner((1 - chi) * s_dot / (1 - s_xav), pybamm.grad(c_s_xav))
            + 1 / self.C_s / (self.phi * (1 - s_xav) ** 2 * (chi * (1 - s_xav) + s_xav) ** 2) 
            * pybamm.div((chi * (1 - s_xav) + s_xav) ** 2 * D_s * pybamm.grad(c_s_xav))
        )
        # self.rhs[c_s_xav] = (
        #     1 / self.C_s / (chi ** 2) 
        #     * pybamm.div(chi ** 2 * D_s * pybamm.grad(c_s_xav))
        # )

        self.rhs[c_o_xav] = (
            pybamm.inner((1 - psi) * s_dot / (1 - s_xav), pybamm.grad(c_o_xav))
            + 1 / self.C_o / ((1 - s_xav) ** 2 * (psi * (1 - s_xav) + s_xav) ** 2)
            * pybamm.div((psi * (1 - s_xav) + s_xav) ** 2 * D_o * pybamm.grad(c_o_xav))
        )
        self.rhs[s_xav] = s_dot 

    def set_boundary_conditions(self, variables):

        c_c_xav = variables["X-averaged positive core concentration"]
        c_s_xav = variables["X-averaged positive shell concentration"]
        c_o_xav = variables["X-averaged positive shell concentration of oxygen"]
        s_xav = variables["X-averaged moving phase boundary location"]

        T_xav = variables["X-averaged positive electrode temperature"]
        T_xav_s = pybamm.PrimaryBroadcast(T_xav, ["positive shell"])

        # # D_s and D_o defined in base_phase_transition
        # D_s = pybamm.boundary_value(self.D_s(c_s_xav, T_xav_s), "right")

        c_c_N_av = variables["X-averaged positive core surface cell concentration"]
        c_s_1_av = variables["X-averaged positive shell center cell concentration"]
        c_o_1_av = variables["X-averaged positive shell center cell concentration of oxygen"]
        dx_cp_av = variables["X-averaged positive core surface cell length"]
        dx_cs_av = variables["X-averaged positive shell center cell length"]
        dx_co_av = variables["X-averaged positive shell center cell length of oxygen"]

        j_xav = variables[
            "X-averaged positive electrode interfacial current density"
        ]

        c_share = variables["X-averaged shared concentration at core-shell interface"]
        c_o_b   = variables["X-averaged oxygen concentration at core-shell interface"]

        rbc_cc = (c_share - c_c_N_av) / dx_cp_av
        lbc_cs = (c_s_1_av - c_share) / dx_cs_av

        rbc_cs = (
            -self.C_s
            * (1 - s_xav)
            * j_xav
            / self.param.a_R_p
            / self.param.gamma_p
            / pybamm.surf(self.D_s(c_s_xav, T_xav_s))
        )

        lbc_co = (c_o_1_av - c_o_b) / dx_co_av

        self.boundary_conditions[c_c_xav] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (rbc_cc, "Neumann")
            # "right": (rbc_cs, "Neumann")
        }
        self.boundary_conditions[c_s_xav] = {
            "left":  (lbc_cs, "Neumann"),
            # "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (rbc_cs, "Neumann")
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
        c_s_xav = variables["X-averaged positive shell concentration"]
        c_o_xav = variables["X-averaged positive shell concentration of oxygen"]
        s_xav = variables["X-averaged moving phase boundary location"]

        self.initial_conditions[c_c_xav] = self.c_c_init(1)
        self.initial_conditions[c_s_xav] = self.c_s_init(1)
        self.initial_conditions[c_o_xav] = self.c_o_init(1)
        self.initial_conditions[s_xav] = self.s_init(1)