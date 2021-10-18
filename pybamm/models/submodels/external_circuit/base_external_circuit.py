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
        # Q_tct = pybamm.standard_variables.Q_tct
        variables = {
            "Discharge capacity [A.h]": Q,
            # "Total charge throughput [A.h]": Q_tct
        }
        return variables

    def set_initial_conditions(self, variables):
        Q = variables["Discharge capacity [A.h]"]
        self.initial_conditions[Q] = pybamm.Scalar(0)
        # Q_tct = variables["Total charge throughput [A.h]"]
        # self.initial_conditions[Q_tct] = pybamm.Scalar(0)

    def set_rhs(self, variables):
        # ODE for discharge capacity
        Q = variables["Discharge capacity [A.h]"]
        # Q_tct = variables["Total charge throughput [A.h]"]
        I = variables["Current [A]"]
        self.rhs[Q] = I * self.param.timescale / 3600
        # self.rhs[Q_tct] = pybamm.abs(I) * self.param.timescale / 3600


class LeadingOrderBaseModel(BaseModel):
    """Model to represent the behaviour of the external circuit, at leading order."""

    def __init__(self, param):
        super().__init__(param)

    def get_fundamental_variables(self):
        Q = pybamm.Variable("Leading-order discharge capacity [A.h]")
        variables = {"Discharge capacity [A.h]": Q}
        return variables
