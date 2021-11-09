#%%
output_variables_spm = [
    # "X-averaged negative particle concentration [mol.m-3]",
    # "X-averaged positive particle concentration [mol.m-3]",
    #---------------------------------------------------------
    "X-averaged positive core concentration [mol.m-3]",
    #---------------------------------------------------------
    # "X-averaged positive shell concentration [mol.m-3]",
    # #---------------------------------------------------------
    # "X-averaged positive shell concentration of oxygen [mol.m-3]",
    #---------------------------------------------------------
    "X-averaged moving phase boundary location",
    # #---------------------------------------------------------
    "X-averaged shared concentration at core-shell interface [mol.m-3]",
    #---------------------------------------------------------
    "Current [A]",
    "Terminal voltage [V]",
    # "X-averaged positive electrode temperature [K]",
    # "X-averaged time derivative of moving phase boundary location",
    # "Loss of active material in positive electrode (MZ)",
    "X-averaged loss of active material in positive electrode (MZ)",
    # "X-averaged positive electrode interfacial current density [A.m-2]",
    # "Loss of active material in positive electrode [%]",
    # "X-averaged positive core surface concentration [mol.m-3]",
    "Discharge capacity [A.h]",
    # "Total charge throughput [A.h]"
]
#%
output_variables_spm.extend(
    [
        #---------------------------------------------------------
        [
            "Total lithium in positive electrode [mol]",
            "Total lithium in negative electrode [mol]",
            "Total lithium [mol]",
        ],
        "Loss of lithium inventory [%]",
        # "Positive electrode interfacial current density [A.m-2]",
        # "Negative electrode interfacial current density [A.m-2]",
        # "X-averaged PE shell layer overpotential [V]",
        "PE shell layer overpotential [V]"
        # "Electrolyte concentration [mol.m-3]",
        # "Negative electrode potential [V]",
        # "Electrolyte potential [V]",
        # "Positive electrode potential [V]",
    ]
)


#%%
which_average = "X"

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
else:
    output_variables_dfn = [
        "Negative particle concentration [mol.m-3]",
        "Positive core concentration [mol.m-3]",
        "Positive shell concentration of oxygen [mol.m-3]"
    ] 
    
    
output_variables_dfn.extend(
    [
        "Moving phase boundary location",
        # "X-averaged moving phase boundary location",
        #---------------------------------------------------------
        "Shared concentration at core-shell interface [mol.m-3]",     
        # "X-averaged shared concentration at core-shell interface [mol.m-3]",
        #---------------------------------------------------------
        # "Current [A]",
        "Terminal voltage [V]",
        # "X-averaged positive electrode temperature [K]",
        # "Time derivative of moving phase boundary location",
        # "X-averaged time derivative of moving phase boundary location",
        "Loss of active material in positive electrode (MZ)",
        # "X-averaged loss of active material in positive electrode (MZ)",
        [
        "Total lithium in positive electrode [mol]",
        "Total lithium in negative electrode [mol]",
        "Total lithium [mol]",
        ],
        "LLI [%]",
        "X-averaged PE shell layer overpotential [V]",
        "PE shell layer overpotential [V]",
        "Positive electrode interfacial current density [A.m-2]",
        # "Electrolyte concentration [mol.m-3]",
        # "Negative electrode potential [V]",
        # "Electrolyte potential [V]",
        # "Positive electrode potential [V]",
    ]
)

    
# output_variables.extend(
#     [
#         "Moving phase boundary location",
#         "X-averaged moving phase boundary location",
#         #---------------------------------------------------------
#         # "Shared concentration at core-shell interface [mol.m-3]",     
#         # "X-averaged shared concentration at core-shell interface [mol.m-3]",
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
#         # "X-averaged shared concentration at core-shell interface [mol.m-3]",
#         # "Moving phase boundary location",
#         # "Shared concentration at core-shell interface [mol.m-3]",     
#         "Loss of lithium inventory [%]",
#         [
#         "Total lithium in positive electrode [mol]",
#         "Total lithium in negative electrode [mol]",
#         "Total lithium [mol]",
#         ],
#         "X-averaged moving phase boundary location",
#         #---------------------------------------------------------
#         # "Shared concentration at core-shell interface [mol.m-3]",     
#         # "X-averaged shared concentration at core-shell interface [mol.m-3]",
#         "X-averaged positive shell surface concentration [mol.m-3]",
#         # "Negative particle concentration [mol.m-3]",
#         # "Positive shell concentration of oxygen [mol.m-3]",
# ]