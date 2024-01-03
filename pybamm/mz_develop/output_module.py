#%%
output_variables_spm = [
    "Current [A]",
    "Terminal voltage [V]",
    # "X-averaged negative particle concentration",
    #---------------------------------------------------------
    # "X-averaged positive core concentration",
    #---------------------------------------------------------
    #---------------------------------------------------------
    # "X-averaged positive shell concentration of oxygen [mol.m-3]",
    # "Positive shell concentration of oxygen [mol.m-3]",
    #---------------------------------------------------------
    "X-averaged moving phase boundary location",
    # #---------------------------------------------------------
    [
         # "X-averaged lithium concentration at core-shell interface [mol.m-3]",
          # "X-averaged positive core surface concentration [mol.m-3]",
          # "X-averaged negative particle surface concentration [mol.m-3]",
          "X-averaged positive core surface concentration",
          "X-averaged negative particle surface concentration",
    ],
    #---------------------------------------------------------
    # "X-averaged positive electrode temperature [K]",
    # "X-averaged time derivative of moving phase boundary location",
    # "Loss of active material in positive electrode (MZ)",
    "X-averaged loss of active material in positive electrode (MZ)",
    # "Loss of active material in positive electrode [%]",
    # "Discharge capacity [A.h]",
    "SoC_dis",
    # "Total charge throughput [A.h]"
]
#%
output_variables_spm.extend(
    [
        #---------------------------------------------------------
        # [
        #     "Total lithium in positive electrode [mol]",
        #     "Total lithium in negative electrode [mol]",
        #     "Total lithium in particles [mol]",
        # ],
        [
            "Total cyclable lithium in positive electrode [mol]",
            "Total cyclable lithium in negative electrode [mol]",
            "Total cyclable lithium in particles [mol]",
        ],
        # "Loss of lithium inventory",
        # "LLI",
        "LLI_cyc",
        # "Positive electrode interfacial current density [A.m-2]",
        # "Negative electrode interfacial current density [A.m-2]",
        # "X-averaged PE shell layer overpotential [V]",
        # "PE shell layer overpotential [V]"
        # "Electrolyte concentration [mol.m-3]",
        # "X-averaged negative electrode potential [V]",
        # "X-averaged positive electrode potential [V]",
        # "Electrolyte potential [V]",
        # "X-averaged positive electrode reaction overpotential [V]",
        # "X-averaged negative electrode reaction overpotential [V]",
        # "X-averaged positive electrode interfacial current density [A.m-2]",
        # "X-averaged negative electrode interfacial current density [A.m-2]",
        # "X-averaged positive electrode exchange current density [A.m-2]",
        # "X-averaged negative electrode exchange current density [A.m-2]",
        # [
        #     "X-averaged positive electrode open circuit potential [V]",
        #     "X-averaged negative electrode open circuit potential [V]",
        # ],
    ]
)


#%%
which_average = "XR"

if which_average == "X":
    output_variables_dfn = [
        "X-averaged negative particle concentration [mol.m-3]",
        "X-averaged positive core concentration [mol.m-3]",
        "X-averaged positive shell concentration of oxygen [mol.m-3]",
    ]
elif which_average == "R":
    output_variables_dfn = [
        "R-averaged negative particle concentration",
        "R-averaged positive core concentration [mol.m-3]",
        "R-averaged positive shell concentration of oxygen [mol.m-3]",
    ]
elif which_average == "XR":
    output_variables_dfn = [
        # "Negative particle concentration [mol.m-3]",
        # "Positive core concentration [mol.m-3]",
        "Positive shell concentration of oxygen [mol.m-3]"
    ] 
else:
    output_variables_dfn = []
    
# output_variables_dfn.extend(
#     [
#         "Current [A]",
#         "Terminal voltage [V]",
#         "Moving phase boundary location",
#         # "X-averaged moving phase boundary location",
#         #---------------------------------------------------------
#         # "Lithium concentration at core-shell interface [mol.m-3]",     
#         # "X-averaged lithium concentration at core-shell interface [mol.m-3]",
#         #---------------------------------------------------------
#         # [
#         #       "X-averaged positive core surface concentration",
#         #       "X-averaged negative particle surface concentration",
#         # ],
#         "Positive core surface concentration [mol.m-3]",
#         "Negative particle surface concentration [mol.m-3]",
#         # "X-averaged positive electrode temperature [K]",
#         # "Time derivative of moving phase boundary location",
#         # "X-averaged time derivative of moving phase boundary location",
#         # "Loss of active material in positive electrode (MZ)",
#         "X-averaged loss of active material in positive electrode (MZ)",
#         # [
#         # "Total lithium in positive electrode [mol]",
#         # "Total lithium in negative electrode [mol]",
#         # "Total lithium in particles [mol]",
#         # ],
#         # [
#         #     "Total cyclable lithium in positive electrode [mol]",
#         #     "Total cyclable lithium in negative electrode [mol]",
#         #     "Total cyclable lithium in particles [mol]",
#         # ],
#         # "LLI",
#         "LLI_cyc",
#         # "Discharge capacity [A.h]",
#         "SoC",
#         # "X-averaged PE shell layer overpotential [V]",
#         # "PE shell layer overpotential [V]",
#         "Positive electrode interfacial current density [A.m-2]",
#         "Negative electrode interfacial current density [A.m-2]",
#         # "Electrolyte concentration [mol.m-3]",
#         # "Negative electrode potential [V]",
#         # "Electrolyte potential [V]",
#         # "Positive electrode potential [V]",
#     ]
# )

#%%
output_variables_spm_par = [
    # "X-averaged negative particle concentration [mol.m-3]",
    # "X-averaged positive particle concentration [mol.m-3]",
    # "X-averaged negative particle concentration",
    # "X-averaged positive particle concentration",
    # "X-averaged negative particle surface concentration",
    "X-averaged positive particle surface concentration",
    #---------------------------------------------------------
    "Current [A]",
    "Terminal voltage [V]",
    # "X-averaged positive electrode temperature [K]",
    # "X-averaged positive electrode interfacial current density [A.m-2]",
    # "Positive electrode interfacial current density [A.m-2]",
    # "Negative electrode interfacial current density [A.m-2]",
    "Discharge capacity [A.h]",
    # "SoC",
    # "Total charge throughput [A.h]"
]
#%
output_variables_spm_par.extend(
    [
        #---------------------------------------------------------
        # [
        #     "Total lithium in positive electrode [mol]",
        #     "Total lithium in negative electrode [mol]",
        #     "Total lithium [mol]",
        # ],
        # # "Loss of lithium inventory",
        # "LLI",
        # "Electrolyte concentration [mol.m-3]",
        # "X-averaged negative electrode open circuit potential [V]",
        # "Electrolyte potential [V]",
        # "X-averaged positive electrode open circuit potential [V]",
        "X-averaged positive electrode reaction overpotential [V]",
        "X-averaged negative electrode reaction overpotential [V]",
        # "X-averaged positive electrode exchange current density [A.m-2]",
        # "X-averaged positive electrode interfacial current density [A.m-2]",
        # "X-averaged negative electrode interfacial current density [A.m-2]",
    ]
)

#%%
output_variables_dfn_par = [
    "R-averaged negative particle concentration [mol.m-3]",
    "R-averaged positive particle concentration [mol.m-3]",
    "Negative particle surface concentration [mol.m-3]",
    "Positive particle surface concentration [mol.m-3]",
    # "X-averaged negative particle concentration [mol.m-3]",
    # "X-averaged positive particle concentration [mol.m-3]",
    #---------------------------------------------------------
    "Current [A]",
    "Terminal voltage [V]",
    # "X-averaged positive electrode temperature [K]",
    # "X-averaged positive electrode interfacial current density [A.m-2]",
    "Charge capacity [A.h]",
    # "Total charge throughput [A.h]"
]
#%
output_variables_dfn_par.extend(
    [
        #---------------------------------------------------------
        # [
        #     "Total lithium in positive electrode [mol]",
        #     "Total lithium in negative electrode [mol]",
        #     "Total lithium [mol]",
        # ],
        # # "Loss of lithium inventory",
        # "LLI",
        # "Electrolyte concentration [mol.m-3]",
        # "Negative electrode potential [V]",
        # "Electrolyte potential [V]",
        # "Positive electrode potential [V]",
        # "Positive electrode reaction overpotential [V]",
        # "Negative electrode reaction overpotential [V]",
        # "Positive electrode surface potential difference [V]",
        # "Negative electrode surface potential difference [V]",
        "Positive electrode interfacial current density [A.m-2]",
        "Negative electrode interfacial current density [A.m-2]",
        # "X-averaged positive electrode interfacial current density [A.m-2]",
        # "X-averaged negative electrode interfacial current density [A.m-2]",
        # "Positive electrode open circuit potential [V]",
        # "Negative electrode open circuit potential [V]",
    ]
)


#%%
# output_variables.extend(
#     [
#         "Moving phase boundary location",
#         "X-averaged moving phase boundary location",
#         #---------------------------------------------------------
#         # "Lithium concentration at core-shell interface [mol.m-3]",     
#         # "X-averaged lithium concentration at core-shell interface [mol.m-3]",
#         "Loss of lithium inventory [%]",
#         #---------------------------------------------------------
#         # "Electrolyte concentration [mol.m-3]",
#         "Current [A]",
#         "Terminal voltage [V]",
#         # "X-averaged positive electrode temperature [K]",
#         # "Time derivative of moving phase boundary location",
#         # "X-averaged time derivative of moving phase boundary location",
#         # [
#         # "Total lithium in positive electrode [mol]",
#         # "Total lithium in negative electrode [mol]",
#         # "Total lithium [mol]",
#         # ],
#         "Negative electrode interfacial current density [A.m-2]",
#         "Positive electrode interfacial current density [A.m-2]",   
#         # "Electrolyte concentration [mol.m-3]",
#         # "Negative electrode potential [V]",
#         # "Electrolyte potential [V]",
#         # "Positive electrode potential [V]",
#     ]
# )



#%%
# output_variables = [
#     # "Negative particle surface concentration [mol.m-3]",
#     "Negative electrode volume-averaged concentration [mol.m-3]",
#     # "Electrolyte concentration [mol.m-3]",
#     # "Positive particle surface concentration [mol.m-3]",
#     "Positive electrode volume-averaged concentration [mol.m-3]",
#     "Current [A]",
#     "Negative electrode surface potential difference [V]",
#     # "Electrolyte potential [V]",
#     "Positive electrode surface potential difference [V]",
#     "Terminal voltage [V]",
#     "X-averaged negative electrode interfacial current density [A.m-2]",
#     "X-averaged positive electrode interfacial current density [A.m-2]",   
#     # "Exchange current density",
#     "Total current density [A.m-2]"
# ]
# 

# output_variables = [
#         "Current [A]",
#         "Terminal voltage [V]",
#         # "X-averaged moving phase boundary location",
#         # "X-averaged lithium concentration at core-shell interface [mol.m-3]",
#         # "Moving phase boundary location",
#         # "Lithium concentration at core-shell interface [mol.m-3]",     
#         "Loss of lithium inventory [%]",
#         [
#         "Total lithium in positive electrode [mol]",
#         "Total lithium in negative electrode [mol]",
#         "Total lithium [mol]",
#         ],
#         "X-averaged moving phase boundary location",
#         #---------------------------------------------------------
#         # "Lithium concentration at core-shell interface [mol.m-3]",     
#         # "X-averaged lithium concentration at core-shell interface [mol.m-3]",
#         "X-averaged positive shell surface concentration [mol.m-3]",
#         # "Negative particle concentration [mol.m-3]",
#         # "Positive shell concentration of oxygen [mol.m-3]",
# ]