#
# Basic Positive Electrode Degradation Model embedded in SPM
#

import numpy as np
import pybamm

# utility

class HelpParameters:
    def __init__(self):

        self.n_cells = 1

        # self.I_typ = pybamm.Parameter("Typical current [A]")
        self.I_typ = 1

        self.Q = pybamm.Parameter("Nominal cell capacity [A.h]")

        # self.n_electrodes_parallel = pybamm.Parameter(
            # "Number of electrodes connected in parallel to make a cell"
        # )
        self.n_electrodes_parallel = 1

        # self.timescale = pybamm.Parameter("Typical timescale [s]")
        # self.dimensional_current_with_time = pybamm.FunctionParameter(
        #     "Current function [A]", {"Time[s]": pybamm.t * self.timescale}
        # )
        self.dimensional_current_with_time = pybamm.Parameter("Current function [A]")

        self.current_with_time = (
            self.dimensional_current_with_time
            / self.I_typ
            * pybamm.Function(np.sign, self.I_typ)
        )

        # self.L_y = pybamm.Parameter("Electrode width [m]")
        # self.L_z = pybamm.Parameter("Electrode height [m]")
        self.A_cc = pybamm.Parameter("Positive electroactive surface area in total [m2]")
        # self.A_cc = pybamm.Parameter("Negative electroactive surface area in total [m2]")


class PeDegradationModel (pybamm.BaseModel):
    """NMC cathode degradation model embedded in Single Particle Model (SPM) model of
    a lithium-ion battery, from [1]_.

    Parameters
    ----------
    name : str, optional
        The name of the model.

    References
    ----------
    .. [1] Ghosh, Abir, et al. "A Shrinking-Core Model for the Degradation of
           High-Nickel Cathodes (NMC811) in Li-Ion Batteries: Passivation Layer
           Growth and Oxygen Evolution." Journal of The Electrochemical Society
           168.2 (2021): 020509.

    **Extends:** :class:`pybamm.lithium_ion.BaseModel`
    """

    def __init__(self, name="PE Degradation Model"):
        super().__init__(name)

        # pybamm.citations.register("Ghosh2021")

        ######################
        # parameters
        ######################

        # param = pybamm.PeDegradationParameters()

        # gas constant R = 8.314 [J.K-1.mol-1]
        R = pybamm.constants.R
        # %
        # Faraday constant F = 96485.33212 [C·mol-1]
        F = pybamm.constants.F

        # physical parameters
        c_p_max = pybamm.Parameter(
            "Maximum lithium concentration in positive electrode [mol.m-3]")
        c_n_max = pybamm.Parameter(
            "Maximum lithium concentration in negative electrode [mol.m-3]")

        # the constant coefficients in the temperature dependent diffusivities
        D_p_ref = pybamm.Parameter(
            "Reference lithium diffusivity in positive particle core [m2.s-1]")
        D_s_ref = pybamm.Parameter(
            "Reference lithium diffusivity in positive particle shell [m2.s-1]")
        D_o_ref = pybamm.Parameter(
            "Reference oxygen diffusivity in positive particle shell [m2.s-1]")
        D_n_ref = pybamm.Parameter(
            "Reference lithium diffusivity in negative electrode [m2.s-1]")

        # geometry parameters
        A_p = pybamm.Parameter(
            "Positive electroactive surface area in total [m2]")
        A_n = pybamm.Parameter(
            "Negative electroactive surface area in total [m2]")
        # radii of positive and negative particles
        R_p = pybamm.Parameter("Positive electrode particle radius [m]")
        R_n = pybamm.Parameter("Negative electrode particle radius [m]")

        # chracteristic time scale as reference
        t_ref = R_p ** 2 / D_p_ref

        # ds / dt = k_1 - k_2 *c_o
        k_1_ref = pybamm.Parameter(
            "Forward chemical reaction (degradation) coefficient [m.s-1]")
        k_2_ref = pybamm.Parameter(
            "Reverse chemical reaction (degradation) coefficient [m4.mol-1.s-1]")

        # electrochemical reaction parameters
        alpha = pybamm.Parameter("Charge transfer coefficient [-]")
        k_p_ref = pybamm.Parameter(
            "Reference exchange-current density in positive particle [m2.5.mol-0.5.s-1]")
        k_n_ref = pybamm.Parameter(
            "Reference exchange-current density in negative particle [m2.5.mol-0.5.s-1]")

        #
        aEne_diff_p = pybamm.Parameter(
            "Activation energy for lithium diffusivity in positive particle core [J.mol-1]")
        aEne_diff_s = pybamm.Parameter(
            "Activation energy for lithium diffusivity in positive particle shell [J.mol-1]")
        aEne_diff_o = pybamm.Parameter(
            "Activation energy for oxygen diffusivity in positive particle shell [J.mol-1]")
        aEne_diff_n = pybamm.Parameter(
            "Activation energy for lithium diffusivity in negative electrode [J.mol-1]")
        # #
        aEne_chem_reac_for = pybamm.Parameter(
            "Activation energy for forward degradation reaction [J.mol-1]")
        aEne_chem_reac_rev = pybamm.Parameter(
            "Activation energy for reverse degradation reaction [J.mol-1]")
        #
        aEne_surf_reac_p = pybamm.Parameter(
            "Activation energy for exchange-current density in positive particle [J.mol-1]")
        aEne_surf_reac_n = pybamm.Parameter(
            "Activation energy for exchange-current density in negative particle [J.mol-1]")

        # non-changing parameters
        c_e = pybamm.Parameter("Electrolyte concentration [mol.m-3]")
        theta_1 = pybamm.Parameter(
            "Electrolyte resistance coefficient 1 [Ohm]")
        theta_2 = pybamm.Parameter(
            "Electrolyte resistance coefficient 2 [Ohm.K-1]")

        # cell parameters
        rho = pybamm.Parameter("Cell density [kg.m-3]")
        vol = pybamm.Parameter("Cell volume [m3]")
        c_rho = pybamm.Parameter(
            "Specific heat capacity of the cell [J.kg-1.K-1]")
        hA = pybamm.Parameter(
            "Heat transfer coefficient x cell surface area [J.s-1.K-1]")

        T_ref = pybamm.Parameter("Reference temperature [K]")
        T_amb = pybamm.Parameter("Ambient temperature [K]")

        # simulation related parameters
        # discharge positive agreeing with pybamm experiment class definition
        I_app = pybamm.Parameter("Current function [A]")

        # upper cut-off voltage
        voltage_high_cut = pybamm.Parameter("Upper cut-off voltage [V]")

        # passivation shell layer Parameters (Positive Electrode)

        phi = pybamm.Parameter(
            "Ratio of maximum lithium concentration between core and shell [-]")

        c_o_core = pybamm.Parameter(
            "Constant oxygen concentration in particle core [mol.m-3]")
        c_o_core_scaled = c_o_core / c_p_max

        # initial conditions
        c_p_0 = pybamm.Parameter(
            "Initial normalized lithium concentration in positive particle core [-]")
        c_s_0 = pybamm.Parameter(
            "Initial normalized lithium concentration in positive particle shell [-]")
        c_o_0 = pybamm.Parameter(
            "Initial normalized oxygen concentration in positive particle shell [-]")
        c_n_0 = pybamm.Parameter(
            "Initial normalized lithium concentration in negative electrode [-]")
        s_0 = pybamm.Parameter(
            "Initial normalized phase boundary location [-]")

        # normalized critical lithium concentration in core at which phase transition may occur
        c_p_thrd = pybamm.Parameter(
            "Normalized threshold lithium concentration in positive particle core [-]")

        ############################################
        # Variables
        ############################################
        c_p = pybamm.Variable(
            "Normalized lithium concentration in positive particle core [-]", domain="positive particle core")
        c_s = pybamm.Variable(
            "Normalized lithium concentration in positive particle shell [-]", domain="positive particle shell")
        c_o = pybamm.Variable(
            "Normalized oxygen concentration in positive particle shell [-]", domain="positive particle shell")
        # normalized by c_n_max
        c_n = pybamm.Variable(
            "Normalized lithium concentration in negative particle [-]", domain="negative particle")

        # normalized moving phase boundary location by R_p
        s = pybamm.Variable("Normalized moving phase boundary location [-]")

        # normalized temp. by T_ref
        T = pybamm.Variable("Normalized temperature [-]")

        # Variables that depend on time only are created without a domain
        Q = pybamm.Variable("Discharge capacity [A.h]")

        # surface value of concentration variables via extrapolation
        c_p_surf = pybamm.boundary_value(c_p, "right")
        c_s_cent = pybamm.boundary_value(c_s, "left")
        c_s_surf = pybamm.boundary_value(c_s, "right")
        c_o_cent = pybamm.boundary_value(c_o, "left")
        c_n_surf = pybamm.boundary_value(c_n, "right")

        c_p_N = pybamm.boundary_cell_value(c_p, "right")
        c_s_1 = pybamm.boundary_cell_value(c_s, "left")
        c_o_1 = pybamm.boundary_cell_value(c_o, "left")

        # boundary cell length: node to the edge
        dx_cp = pybamm.boundary_cell_length(c_p, "right")
        dx_cs = pybamm.boundary_cell_length(c_s, "left")
        dx_co = pybamm.boundary_cell_length(c_o, "left")

        ############################################
        # temp dependent diffusivity
        ############################################
        # arrhenius takes non-scaled temp. with unit [K]
        D_p = D_p_ref * pybamm.arrhenius(aEne_diff_p, T_ref, T * T_ref)
        D_s = D_s_ref * pybamm.arrhenius(aEne_diff_s, T_ref, T * T_ref)
        D_o = D_o_ref * pybamm.arrhenius(aEne_diff_o, T_ref, T * T_ref)
        D_n = D_n_ref * pybamm.arrhenius(aEne_diff_n, T_ref, T * T_ref)

        # we use D_p_ref to normalize other diffusivities
        D_p_scaled = D_p / D_p_ref
        D_s_scaled = D_s / D_p_ref
        D_o_scaled = D_o / D_p_ref
        D_n_scaled = D_n / D_p_ref

        k_1 = k_1_ref * pybamm.arrhenius(aEne_chem_reac_for, T_ref, T * T_ref)
        k_2 = k_2_ref * pybamm.arrhenius(aEne_chem_reac_rev, T_ref, T * T_ref)

        k_1_scaled = k_1 * t_ref / R_p
        k_2_scaled = k_2 * t_ref / R_p * c_p_max

        # time derivative of phase boundary location
        # / c_o_core_scaled is used in line with Abir matlab code, to pick same param values
        # not in the model 
        s_dot = -(k_1_scaled - k_2_scaled * c_o_cent) / c_o_core_scaled * \
            pybamm.EqualHeaviside(c_p_surf, c_p_thrd)
        # s_dot = pybamm.Scalar(0)

        k_p = k_p_ref * pybamm.arrhenius(aEne_surf_reac_p, T_ref, T * T_ref)
        k_n = k_n_ref * pybamm.arrhenius(aEne_surf_reac_n, T_ref, T * T_ref)

        # boundary flux
        # molar flux of lithium at the boundary [mol.m-2.s-1]
        # discharge positive in agreement with pybamm experiment class definition
        J_p = - I_app / (F * A_p)
        J_n = I_app / (F * A_n)

        # scaled flux by nondimensionalization and variable change
        J_p_scaled = (1 - s) * J_p * R_p / (D_p_ref * c_p_max)
        J_n_scaled = J_n * R_n / (D_p_ref * c_n_max)

        ############################################
        # spatial variables
        ############################################
        eta = pybamm.SpatialVariable(
            "eta", domain=["positive particle core"], coord_sys="cartesian")
        chi = pybamm.SpatialVariable(
            "chi", domain=["positive particle shell"], coord_sys="cartesian")
        r_n = pybamm.SpatialVariable(
            "r_n", domain=["negative particle"], coord_sys="spherical polar")
        # r_n = pybamm.SpatialVariable("r_n", domain=["negative particle"], coord_sys="cartesian")

        ############################################
        # set timescale and length scales
        ############################################

        self.timescale = t_ref
        self.length_scales = {
            # "positive particle core": R_p,
            # "positive particle shell": R_p,
            # "negative particle": R_n,

            # default length scale is 1 m, change to 1 um
            # to fool pybamm countercancel spatial_factor = 1e6 
            # default plot unit
            # so that x [um] is actually dimensionless
            "positive particle core": pybamm.Scalar(1e-6),
            "positive particle shell": pybamm.Scalar(1e-6),
            "negative particle": pybamm.Scalar(1e-6)
        }

        ############################################
        # rhs of equations
        ############################################

        self.rhs[c_p] = (
            pybamm.inner(eta * s_dot / s, pybamm.grad(c_p)) + D_p_scaled / (
                (eta * s) ** 2) * pybamm.div(eta ** 2 * pybamm.grad(c_p))
        )

        self.rhs[c_s] = (
            pybamm.inner((1 - chi) * s_dot / (1 - s), pybamm.grad(c_s))
            + D_s_scaled / (phi * (1 - s) ** 2 * (chi * (1 - s) + s) ** 2) *
            pybamm.div((chi * (1 - s) + s) ** 2 * pybamm.grad(c_s))
        )

        self.rhs[c_o] = (
            pybamm.inner((1 - chi) * s_dot / (1 - s), pybamm.grad(c_o))
            + D_o_scaled / ((1 - s) ** 2 * (chi * (1 - s) + s) ** 2) *
            pybamm.div((chi * (1 - s) + s) ** 2 * pybamm.grad(c_o))
        )

        self.rhs[c_n] = (
            (R_p / R_n) ** 2 * D_n_scaled *
            pybamm.div(pybamm.grad(c_n))
        )

        self.rhs[s] = s_dot

        self.rhs[Q] = I_app * t_ref / 3600

        ############################################
        # cell voltage and temp rhs
        ############################################
        inputs_p = {"Positive particle stoichiometry": c_s_surf}
        inputs_n = {"Negative particle stoichiometry": c_n_surf}

        # dUdT (derivative of OCV wit respect to Temp), entropy at reference temp.
        dUdT_p = 1e-3 * \
            pybamm.FunctionParameter(
                "Positive electrode OCP entropy coefficient [V.T-1]", inputs_p)
        dUdT_n = 1e-3 * \
            pybamm.FunctionParameter(
                "Negative electrode OCP entropy coefficient [V.T-1]", inputs_n)
        # OCP (open-circuit potential) as a function of concentration at reference temp.
        U_p_ref = pybamm.FunctionParameter(
            "Positive electrode OCP at reference temp [V]", inputs_p)
        U_n_ref = pybamm.FunctionParameter(
            "Negative electrode OCP at reference temp [V]", inputs_n)

        # electrode open-circuit potential as function of stochiometry and temp
        U_p = U_p_ref + dUdT_p * T_ref * (T - 1)
        U_n = U_n_ref + dUdT_n * T_ref * (T - 1)

        m_p = -I_app / (F * A_p * k_p * c_e ** 0.5 * c_p_max *
                        (1 - c_s_surf) ** 0.5 * c_s_surf ** 0.5)
        m_n = -I_app / (F * A_n * k_n * c_e ** 0.5 * c_n_max *
                        (1 - c_n_surf) ** 0.5 * c_n_surf ** 0.5)
        # overpotential [V], temp back to dimensional one
        eta_p = 2 * R * (T * T_ref) / F * \
            pybamm.log(((m_p ** 2 + 4) ** 0.5 + m_p) / 2)
        eta_n = 2 * R * (T * T_ref) / F * \
            pybamm.log(((m_n ** 2 + 4) ** 0.5 + m_n) / 2)
        # temp back to dimensional one
        R_cell = theta_1 + theta_2 * (T * T_ref - T_amb)

        V_cell = U_p - U_n + eta_p + eta_n - I_app * R_cell

        # self.rhs[T] = pybamm.Scalar(0)
        self.rhs[T] = (
            R_p ** 2 / (rho * vol * c_rho * D_p_ref)
            * (-I_app * (dUdT_p - dUdT_n) * T - I_app / T_ref * (V_cell - U_p + U_n) - hA * (T - T_amb / T_ref))
        )

        ############################################
        # boundary conditions
        ############################################

        # c_shared = (
        #     (D_p_scaled * (1 - s) / dx_cp * c_p_N + D_s_scaled * s / dx_cs * c_s_1
        #         + s_dot * s * (1 - s) * c_o_core_scaled)
        #     / (D_p_scaled * (1 - s) / dx_cp + D_s_scaled * s / dx_cs
        #         + s_dot * s * (1 - s) * (1 - phi))
        # )

        # no lithium eaten during phase transition
        c_shared = (
            (D_p_scaled * (1 - s) / dx_cp * c_p_N + D_s_scaled * s / dx_cs * c_s_1)
            / (D_p_scaled * (1 - s) / dx_cp + D_s_scaled * s / dx_cs)
        )

        # only C_oc chemical reaction effect
        # c_shared = (
        #     (D_p_scaled * (1 - s) / dx_cp * c_p_N + D_s_scaled * s / dx_cs * c_s_1
        #         + s_dot * s * (1 - s) * c_o_core_scaled)
        #     / (D_p_scaled * (1 - s) / dx_cp + D_s_scaled * s / dx_cs)
        # )

        # check that what eats the lithium
        # with negative c_o_core, lithium increased
        # c_shared = (
        #     (D_p_scaled * (1 - s) / dx_cp * c_p_N + D_s_scaled * s / dx_cs * c_s_1 
        #         + s_dot * s * (1 - s) * (-c_o_core_scaled/5))
        #     / (D_p_scaled * (1 - s) / dx_cp + D_s_scaled * s / dx_cs)
        # )

        rbc_cp = (c_shared - c_p_N) / dx_cp
        lbc_cs = (c_s_1 - c_shared) / dx_cs

        self.boundary_conditions[c_p] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (rbc_cp, "Neumann")
        }

        self.boundary_conditions[c_s] = {
            "left":  (lbc_cs, "Neumann"),
            "right": (-J_p_scaled / D_s_scaled, "Neumann")
            # "right": (pybamm.Scalar(0), "Neumann")
        }

        c_o_b = ((1 - s) * s_dot * dx_co * c_o_core_scaled -
                 D_o_scaled * c_o_1) / ((1 - s) * s_dot * dx_co - D_o_scaled)
        lbc_co = (c_o_1 - c_o_b) / dx_co
        # or
        # lbc_co = (1 - s) * s_dot * (c_o_1 - c_o_core_scaled) / (dx_co * (1 - s) * s_dot - D_o_scaled)

        self.boundary_conditions[c_o] = {
            # "left":  (s_dot * (1 - s) * (c_o_core_scaled - c_o_cent) / D_o_scaled, "Neumann"),
            "left":  (lbc_co, "Neumann"),
            "right": (pybamm.Scalar(0), "Dirichlet")
        }

        self.boundary_conditions[c_n] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (-J_n_scaled / D_n_scaled, "Neumann")
        }

        ############################################
        # initial conditions
        ############################################

        self.initial_conditions[c_p] = c_p_0
        self.initial_conditions[c_s] = c_s_0
        self.initial_conditions[c_o] = c_o_0 * (1 - chi ** 2)
        self.initial_conditions[c_n] = c_n_0
        self.initial_conditions[s] = s_0
        self.initial_conditions[T] = T_amb / T_ref
        self.initial_conditions[Q] = pybamm.Scalar(0)

        ############################################
        # Events
        ############################################
        # voltage_high_cut = 4.2
        self.events += [
            pybamm.Event(
                "Maximum voltage",
                V_cell - voltage_high_cut,
                pybamm.EventType.TERMINATION
            ),
            # pybamm.Event(
            #     "Minimum negative particle surface concentration",
            #     pybamm.min(c_n_surf) - 0.01,
            # ),
            # pybamm.Event(
            #     "Maximum negative particle surface concentration",
            #     (1 - 0.01) - pybamm.max(c_n_surf),
            # ),
            # pybamm.Event(
            #     "Minimum positive shell surface concentration",
            #     pybamm.min(c_s_surf) - 0.015,
            # )
        ]

        # self.param = pybamm.LithiumIonParameters()
        self.param = HelpParameters()
        param = self.param


        ############################################
        # The `variables` dictionary contains all variables that might be useful for
        # visualising the solution of the model
        ############################################
        self.variables = {
            "Normalized lithium concentration in positive particle core [-]": c_p,
            "Normalized lithium concentration in positive core-shell interface - core [-]":  c_p_surf,
            "Normalized lithium concentration in positive particle shell [-]": c_s,
            "Normalized lithium concentration in positive core-shell interface - shell [-]": c_s_cent,
            "Normalized lithium concentration in positive shell surface [-]": c_s_surf,
            "Normalized oxygen concentration in positive particle shell [-]": c_o,
            "Normalized oxygen concentration in positive core-shell interface [-]": c_o_cent,
            "Normalized lithium concentration in negative particle [-]": c_n,
            "Normalized lithium concentration in negative particle surface [-]": c_n_surf,
            "Normalized moving phase boundary location [-]": s,
            "Normalized temperature [-]": T,
            "Temperature [Kelvin]": T * T_ref,
            "Normalized interface concentration [-]": c_shared,
            "Terminal voltage [V]": V_cell,
            "J_p_scaled": J_p_scaled,
            "J_n_scaled": J_n_scaled,
            "Total current density": I_app / param.I_typ,
            "Total current density [A.m-2]": I_app / A_p,
            "Current [A]":  I_app,
            "C-rate": I_app / param.Q,
            "Discharge capacity [A.h]": Q
        }


    @property
    def t_end(self):
        t_end = pybamm.Parameter("Total simulation time (half C) [s]")
        return self.default_parameter_values.evaluate(t_end)

    ############################################
    # default set up
    ############################################
    @property
    def eta(self):
        return pybamm.SpatialVariable(
            "eta", domain=["positive particle core"], coord_sys="cartesian"
        )

    @property
    def chi(self):
        return pybamm.SpatialVariable(
            "chi", domain=["positive particle shell"], coord_sys="cartesian"
        )

    @property
    def r_n(self):
        return pybamm.SpatialVariable(
            "r_n", domain=["negative particle"], coord_sys="spherical polar"
        )

    @property
    def default_parameter_values(self):
        return pybamm.ParameterValues(
            "pybamm/myDevelop/data/parameters_abir.csv"
        )

    @property
    def default_var_pts(self):
        return {self.eta: 20, self.chi: 5, self.r_n: 20}

    @property
    def default_geometry(self):
        return {
            "negative particle": {self.r_n: {"min": pybamm.Scalar(0.0), "max": pybamm.Scalar(1)}},
            "positive particle core": {self.eta: {"min": pybamm.Scalar(0.0), "max": pybamm.Scalar(1)}},
            "positive particle shell": {self.chi: {"min": pybamm.Scalar(0.0), "max": pybamm.Scalar(1)}}
        }

    @property
    def default_submesh_types(self):
        return {
            "positive particle core": pybamm.MeshGenerator(pybamm.Uniform1DSubMesh),
            "positive particle shell": pybamm.MeshGenerator(pybamm.Uniform1DSubMesh),
            "negative particle": pybamm.MeshGenerator(pybamm.Uniform1DSubMesh)
        }

    @property
    def default_spatial_methods(self):

        return {
            "positive particle core": pybamm.FiniteVolume(),
            "positive particle shell": pybamm.FiniteVolume(),
            "negative particle": pybamm.FiniteVolume()
        }

    @property
    def default_solver(self):
        """Return default solver based on whether model is ODE model or DAE model."""
        return pybamm.CasadiSolver(mode="safe")
        # return pybamm.ScipySolver()

    def new_empty_copy(self):
        return pybamm.BaseModel.new_empty_copy(self)
