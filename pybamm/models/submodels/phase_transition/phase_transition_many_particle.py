#
# Class for many particles with cathode phase transition.
#
import pybamm
from .base_phase_transition import BaseTransition

class PhaseTransitionManyParticle(BaseTransition):
    """
    Class for phase transition in many particles of Positive Electrode for degradation.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str
        The domain of the model either 'Negative' or 'Positive'

    **Extends:** :class:`pybamm.particle.BaseParticle`
    """

    def __init__(self, param, domain="Positive"):
        # check whether the domain is correct
        if domain != "Positive":
            raise ValueError(
                "domain must be 'Positive' for phase transition degradation"
            )
        super().__init__(param, domain)

    def get_fundamental_variables(self):

        # concentration in particle core
        c_c = pybamm.Variable(
            "Positive core concentration",
            domain="positive core",
            auxiliary_domains={
                "secondary": "positive electrode",
                "tertiary": "current collector",
            },
            bounds=(0, 1),
        )

        # concentration in degraded passivation layer/shell
        c_s = pybamm.Variable(
            "Positive shell concentration",
            domain="positive shell",
            auxiliary_domains={
                "secondary": "positive electrode",
                "tertiary": "current collector",
            },
            bounds=(0, 1),
        )

        # oxygen concentration in degraded passivation layer/shell
        c_o = pybamm.Variable(
            "Positive shell concentration of oxygen",
            domain="positive shell oxygen", # use different 
            auxiliary_domains={
                "secondary": "positive electrode",
                "tertiary": "current collector",
            },
            bounds=(0, 1),
        )

        # location of two-phase boundary
        s = pybamm.Variable(
            "Moving phase boundary location",
            domain="positive electrode",
            auxiliary_domains={
                "secondary": "current collector",
            },
            bounds=(0, 1),
        )

        variables = self._get_standard_concentration_variables(
            c_c=c_c, c_s=c_s, c_o=c_o, s=s
        )

        return variables

    def get_coupled_variables(self, variables):

        T = variables["Positive electrode temperature"]

        # get s_dot

        c_o_cent = variables[
            "Positive shell center concentration of oxygen"
        ]
        c_c_surf = variables["Positive core surface concentration"]

        # s_dot = ( 
        #     -(self.kappa_1 * self.k_1(T) - self.kappa_2 * self.k_2(T) * c_o_cent)
        #     * pybamm.EqualHeaviside(c_c_surf, self.c_p_thrd)
        # )
        s_dot = pybamm.Scalar(0)

        # # get c_share and c_o at boundary
        c_c = variables["Positive core concentration"]
        c_s = variables["Positive shell concentration"]
        c_o = variables["Positive shell concentration of oxygen"]
        s = variables["Moving phase boundary location"]

        T_c = pybamm.PrimaryBroadcast(T, ["positive core"])
        T_s = pybamm.PrimaryBroadcast(T, ["positive shell"])
        T_o = pybamm.PrimaryBroadcast(T, ["positive shell oxygen"])

        # defined in base_phase_transition
        D_c = pybamm.surf(self.D_c(c_c, T_c))
        D_s = pybamm.boundary_value(self.D_s(c_s, T_s), "left")
        D_o = pybamm.boundary_value(self.D_o(c_o, T_o), "left")

        c_c_N = variables["Positive core surface cell concentration"]
        c_s_1 = variables["Positive shell center cell concentration"]
        c_o_1 = variables["Positive shell center cell concentration of oxygen"]
        dx_cp = variables["Positive core surface cell length"]
        dx_cs = variables["Positive shell center cell length"]
        dx_co = variables["Positive shell center cell length of oxygen"]

        # the interface boundary value is calculated from applied boundary condition
        # not extrapolated afterwards
        c_share = (
            (D_c * (1 - s) / self.C_c / dx_cp * c_c_N 
             + D_s * s / self.C_s / dx_cs * c_s_1
             + s_dot * s * (1 - s) * self.c_o_core)
            / (D_c * (1 - s)/ self.C_c / dx_cp 
               + D_s * s / self.C_s / dx_cs
               + s_dot * s * (1 - s) * (1 - self.phi))
        )
        # no lithium loss
        # c_share = (
        #     (D_c * (1 - s) / self.C_c / dx_cp * c_c_N 
        #      + D_s * s / self.C_s / dx_cs * c_s_1)
        #     / (D_c * (1 - s)/ self.C_c / dx_cp 
        #        + D_s * s / self.C_s / dx_cs
        #        + s_dot * s * (1 - s) * (1 - self.phi))
        # )
        
        # c_share = (
        #     (D_c * (1 - s) / self.C_c / dx_cp * c_c_N 
        #      + D_s * s / self.C_s / dx_cs * c_s_1)
        #     / (D_c * (1 - s)/ self.C_c / dx_cp 
        #        + D_s * s / self.C_s / dx_cs)
        # )
        c_share_xav = pybamm.x_average(c_share)

        # boundary oxygen concentration from applied bc
        # not extrapolated afterwards
        c_o_b = (
            ((1 - s) * s_dot * dx_co * self.c_o_core - D_o * c_o_1)
            / ((1 - s) * s_dot * dx_co - D_o)
        )
        # c_o_b = (
        #     ((1 - s) * s_dot * dx_co * self.c_o_core - D_o / self.C_o * c_o_1)
        #     / ((1 - s) * s_dot * dx_co - D_o / self.C_o)
        # )
        c_o_b_xav = pybamm.x_average(c_o_b)

        variables.update(
            {
                "Time derivative of moving phase boundary location": s_dot,
                "Shared concentration at core-shell interface": c_share,
                "Shared concentration at core-shell interface [mol.m-3]": 
                    c_share * self.param.c_p_max,
                "X-averaged shared concentration at core-shell interface": c_share_xav,
                "X-averaged shared concentration at core-shell interface [mol.m-3]": 
                    c_share_xav * self.param.c_p_max,
                "Oxygen concentration at core-shell interface": c_o_b,
                "Oxygen concentration at core-shell interface [mol.m-3]": 
                    c_o_b * self.param.c_p_max,
                "X-averaged oxygen concentration at core-shell interface": c_o_b_xav,
                "X-averaged oxygen concentration at core-shell interface [mol.m-3]": 
                    c_o_b_xav * self.param.c_p_max,
            }
        )

        variables.update(self._get_total_concentration_variables(variables))

        return variables

    def set_rhs(self, variables):

        c_c = variables["Positive core concentration"]
        c_s = variables["Positive shell concentration"]
        c_o = variables["Positive shell concentration of oxygen"]
        s = variables["Moving phase boundary location"]
        R = variables["Positive particle radius"]

        T = variables["Positive electrode temperature"]
        T_c = pybamm.PrimaryBroadcast(T, ["positive core"])
        T_s = pybamm.PrimaryBroadcast(T, ["positive shell"])
        T_o = pybamm.PrimaryBroadcast(T, ["positive shell oxygen"])

        # use FunctionParameter "Positive electrode diffusivity [m2.s-1]"
        # in lithium_ion_parameters for the one of core phase
        # defined in base_phase_transition
        D_c = self.D_c(c_c, T_c)
        D_s = self.D_s(c_s, T_s)
        D_o = self.D_o(c_o, T_o)

        eta = pybamm.phase_transition.eta
        chi = pybamm.phase_transition.chi
        psi = pybamm.phase_transition.psi

        s_dot = variables["Time derivative of moving phase boundary location"]

        self.rhs[c_c] = (
            pybamm.inner(eta * s_dot / s, pybamm.grad(c_c)) 
            + 1 / (R ** 2 * self.C_c) / ((eta * s) ** 2) 
            * pybamm.div(eta ** 2 * D_c * pybamm.grad(c_c))
        )
        self.rhs[c_s] = (
            pybamm.inner((1 - chi) * s_dot / (1 - s), pybamm.grad(c_s))
            + 1 / (R ** 2 * self.C_s) / (self.phi * (1 - s) ** 2 * (chi * (1 - s) + s) ** 2) 
            * pybamm.div((chi * (1 - s) + s) ** 2 * D_s * pybamm.grad(c_s))
        )
        self.rhs[c_o] = (
            pybamm.inner((1 - psi) * s_dot / (1 - s), pybamm.grad(c_o))
            + 1 / (R ** 2 * self.C_o) / ((1 - s) ** 2 * (psi * (1 - s) + s) ** 2)
            * pybamm.div((psi * (1 - s) + s) ** 2 * D_o * pybamm.grad(c_o))
        )
        self.rhs[s] = s_dot


    def set_boundary_conditions(self, variables):

        c_c = variables["Positive core concentration"]
        c_s = variables["Positive shell concentration"]
        c_o = variables["Positive shell concentration of oxygen"]
        s = variables["Moving phase boundary location"]

        T = variables["Positive electrode temperature"]
        T_s = pybamm.PrimaryBroadcast(T, ["positive shell"])

        # # D_s and D_o defined in base_phase_transition
        # D_s = pybamm.boundary_value(self.D_s(c_s, T_s), "right")

        c_c_N = variables["Positive core surface cell concentration"]
        c_s_1 = variables["Positive shell center cell concentration"]
        c_o_1 = variables["Positive shell center cell concentration of oxygen"]
        dx_cp = variables["Positive core surface cell length"]
        dx_cs = variables["Positive shell center cell length"]
        dx_co = variables["Positive shell center cell length of oxygen"]

        j = variables["Positive electrode interfacial current density"]
        R = variables["Positive particle radius"]

        c_share = variables["Shared concentration at core-shell interface"]
        c_o_b   = variables["Oxygen concentration at core-shell interface"]

        rbc_cc = (c_share - c_c_N) / dx_cp
        lbc_cs = (c_s_1 - c_share) / dx_cs

        rbc_cs = (
            -self.C_s
            * (1 - s)
            * j
            * R
            / self.param.a_R_p
            / self.param.gamma_p
            / pybamm.surf(self.D_s(c_s, T_s))
        )

        lbc_co = (c_o_1 - c_o_b) / dx_co

        self.boundary_conditions[c_c] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (rbc_cc, "Neumann")
        }
        self.boundary_conditions[c_s] = {
            "left":  (lbc_cs, "Neumann"),
            "right": (rbc_cs, "Neumann")
        }
        self.boundary_conditions[c_o] = {
            "left":  (lbc_co, "Neumann"),
            "right": (pybamm.Scalar(0), "Dirichlet")
        }


    def set_initial_conditions(self, variables):

        c_c = variables["Positive core concentration"]
        c_s = variables["Positive shell concentration"]
        c_o = variables["Positive shell concentration of oxygen"]
        s = variables["Moving phase boundary location"]

        x_p = pybamm.standard_spatial_vars.x_p
        # x_p_p = pybamm.PrimaryBroadcast(x_p, "positive core")
        # x_p_s = pybamm.PrimaryBroadcast(x_p, "positive shell")
        # x_p_o = pybamm.PrimaryBroadcast(x_p, "positive shell oxygen")

        eta = pybamm.phase_transition.eta
        chi = pybamm.phase_transition.chi
        psi = pybamm.phase_transition.psi

        self.initial_conditions[c_c] = self.c_c_init(eta)
        self.initial_conditions[c_s] = self.c_s_init(chi)
        self.initial_conditions[c_o] = self.c_o_init(psi)
        self.initial_conditions[s] = self.s_init(x_p)
