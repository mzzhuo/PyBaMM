#%%
bol = 1
eol = 20

solu_bol = solution.cycles[bol-1]
solu_eol = solution.cycles[eol-1]

#%%
V_t_bol = solu_bol["Terminal voltage [V]"].entries
V_t_eol = solu_eol["Terminal voltage [V]"].entries
#%%
time_bol_3      = solu_bol.steps[3]["Time"].entries
soc_bol_3       = solu_bol.steps[3]["SoC"].entries
c_p_bol_3       = solu_bol.steps[3]["X-averaged positive core surface concentration"].entries
V_p_bol_3       = solu_bol.steps[3]["X-averaged positive electrode open circuit potential [V]"].entries
c_n_bol_3       = solu_bol.steps[3]["X-averaged negative particle surface concentration"].entries
V_n_bol_3       = solu_bol.steps[3]["X-averaged negative electrode open circuit potential [V]"].entries
V_t_bol_3       = solu_bol.steps[3]["Terminal voltage [V]"].entries
eta_shell_bol_3 = solu_bol.steps[3]["X-averaged PE shell layer overpotential [V]"].entries

soc_eol_3       = solu_eol.steps[3]["SoC"].entries
c_p_eol_3       = solu_eol.steps[3]["X-averaged positive core surface concentration"].entries
V_p_eol_3       = solu_eol.steps[3]["X-averaged positive electrode open circuit potential [V]"].entries
c_n_eol_3       = solu_eol.steps[3]["X-averaged negative particle surface concentration"].entries
V_n_eol_3       = solu_eol.steps[3]["X-averaged negative electrode open circuit potential [V]"].entries
V_t_eol_3       = solu_eol.steps[3]["Terminal voltage [V]"].entries
eta_shell_eol_3 = solu_eol.steps[3]["X-averaged PE shell layer overpotential [V]"].entries

#%%
time_bol_2      = solu_bol.steps[2]["Time"].entries
soc_bol_2       = solu_bol.steps[2]["SoC"].entries
c_p_bol_2       = solu_bol.steps[2]["X-averaged positive core surface concentration"].entries
V_p_bol_2       = solu_bol.steps[2]["X-averaged positive electrode open circuit potential [V]"].entries
c_n_bol_2       = solu_bol.steps[2]["X-averaged negative particle surface concentration"].entries
V_n_bol_2       = solu_bol.steps[2]["X-averaged negative electrode open circuit potential [V]"].entries
V_t_bol_2       = solu_bol.steps[2]["Terminal voltage [V]"].entries
eta_shell_bol_2 = solu_bol.steps[2]["X-averaged PE shell layer overpotential [V]"].entries

soc_eol_2       = solu_eol.steps[2]["SoC"].entries
c_p_eol_2       = solu_eol.steps[2]["X-averaged positive core surface concentration"].entries
V_p_eol_2       = solu_eol.steps[2]["X-averaged positive electrode open circuit potential [V]"].entries
c_n_eol_2       = solu_eol.steps[2]["X-averaged negative particle surface concentration"].entries
V_n_eol_2       = solu_eol.steps[2]["X-averaged negative electrode open circuit potential [V]"].entries
V_t_eol_2       = solu_eol.steps[2]["Terminal voltage [V]"].entries
eta_shell_eol_2 = solu_eol.steps[2]["X-averaged PE shell layer overpotential [V]"].entries

#%%
V_t_bol_3[0] = V_t_bol_2[-1]
V_t_eol_3[0] = V_t_eol_2[-1]

#%%
plt.figure(figsize=(8, 6))
plt.plot(soc_bol, V_p_bol, '-', mfc='none', label=str(bol))
plt.plot(soc_bol, V_n_bol, '-', mfc='none')
plt.plot(soc_bol, V_t_bol, '-', mfc='none')
plt.plot(soc_bol, eta_shell_bol, '-', mfc='none')

plt.plot(soc_eol, V_p_eol, '--', mfc='none', label=str(eol))
plt.plot(soc_eol, V_n_eol, '--', mfc='none')
plt.plot(soc_eol, V_t_eol, '--', mfc='none')
plt.plot(soc_eol, eta_shell_eol, '--', mfc='none')

plt.axis([1,0, 0,4.5])
plt.legend()

#%%
f = open("pybamm\mz_develop\output\sce1_stochishift_bol.txt", "w")
np.savetxt(f, np.c_[soc_bol_3, V_p_bol_3, V_n_bol_3, V_t_bol_3, eta_shell_bol_3, c_p_bol_3, c_n_bol_3], fmt='%1.6e', delimiter=', ')
f.close()
#%%
f = open("pybamm\mz_develop\output\sce1_stochishift_eol.txt", "w")
np.savetxt(f, np.c_[soc_eol_3, V_p_eol_3, V_n_eol_3, V_t_eol_3, eta_shell_eol_3, c_p_eol_3, c_n_eol_3], fmt='%1.6e', delimiter=', ')
f.close()

#%%
eol = 20
soc_eol = solution.cycles[eol-1].steps[1]["SoC"].entries

#%%
# time_bol_2 = solution.cycles[bol-1].steps[2]["Time"].entries

# V_t_bol_2 = solution.cycles[bol-1].steps[2]["Terminal voltage [V]"].entries

# time = solution.cycles[bol-1]["Time"].entries
# current = solution.cycles[bol-1]["Current [A]"].entries
# voltage = solution.cycles[bol-1]["Terminal voltage [V]"].entries


