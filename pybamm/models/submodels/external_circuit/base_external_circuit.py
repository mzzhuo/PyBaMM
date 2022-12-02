#
# Base model for the external circuit
#
import pybamm


class BaseModel(pybamm.BaseSubModel):
    """Model to represent the behaviour of the external circuit."""

    def __init__(self, param):
        super().__init__(param)

    def get_fundamental_variables(self):
        Q = pybamm.standard_variables.Q
        # Q_thrpt = pybamm.standard_variables.Q_thrpt
        variables = {
            "Discharge capacity [A.h]": Q,
            "Charge capacity [A.h]": -Q,
            # "Total charge throughput [A.h]": Q_thrpt
            # divided by "Nominal cell capacity [A.h]"
            "SoC_chg": -Q / self.param.Q,
            "SoC_dis": 1 - Q / self.param.Q,
        }
        return variables

    def set_initial_conditions(self, variables):
        Q = variables["Discharge capacity [A.h]"]
        self.initial_conditions[Q] = pybamm.Scalar(0)
        # Q_thrpt = variables["Total charge throughput [A.h]"]
        # self.initial_conditions[Q_thrpt] = pybamm.Scalar(0)

    def set_rhs(self, variables):
        # ODE for discharge capacity
        Q = variables["Discharge capacity [A.h]"]
        # Q_thrpt = variables["Total charge throughput [A.h]"]
        I = variables["Current [A]"]
        self.rhs[Q] = I * self.param.timescale / 3600
        # self.rhs[Q_thrpt] = pybamm.abs(I) * self.param.timescale / 3600


class LeadingOrderBaseModel(BaseModel):
    """Model to represent the behaviour of the external circuit, at leading order."""

    def __init__(self, param):
        super().__init__(param)

    def get_fundamental_variables(self):
        Q = pybamm.Variable("Leading-order discharge capacity [A.h]")
        variables = {"Discharge capacity [A.h]": Q}
        return variables
