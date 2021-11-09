    # --------------------------------------------------------------------
    # voltage correction due to the shell layer
    # --------------------------------------------------------------------

    def _get_corrected_terminal_voltage(self, variables):

        pot_scale = self.param.potential_scale 
        R = pybamm.Parameter("Shell resistance [ohm]")
        I = variables["Current [A]"]

        V = variables["Terminal voltage"]
        V_dim = variables["Terminal voltage [V]"]

        # print(I, R)

        variables.update(
            {
                "Terminal voltage": V - I * R / pot_scale,
                "Terminal voltage [V]": V_dim - I * R,
            }
        )
        return variables 