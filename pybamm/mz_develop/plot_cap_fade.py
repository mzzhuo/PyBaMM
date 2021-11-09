#%%
timescale = solution.timescale_eval
time_in_sec = solution.cycles[0].t * timescale

V_cell = solution.cycles[0]["Terminal voltage [V]"].entries
I = solution.cycles[0]["Current [A]"].entries
Q = solution.cycles[0]["Discharge capacity [A.h]"].entries

#%%
plt.figure(figsize=(8, 6))
plt.plot(time_in_sec, V_cell_out_value, 'b-')


#%% loop to calculate discharge capacity of each cycle
import numpy as np
total_cycles = len(solution.cycles)
Q_dis_cycles = []
cycle_numbers = np.arange(1, total_cycles+1)
for i in range(total_cycles):
    Q_dis_cyc = solution.cycles[i].steps[3]["Discharge capacity [A.h]"].entries
    Q_dis = Q_dis_cyc[-1] - Q_dis_cyc[0]
    Q_dis_cycles.append(Q_dis)
    
#%%
plt.figure(figsize=(8, 6))
# markerfacecolor
plt.plot(cycle_numbers, Q_dis_cycles, 'o', mfc='none')

#%% calculate total charge throughput
from scipy import integrate

I = solution["Current [A]"].entries
timescale = solution.timescale_eval
time_in_sec = solution.t * timescale

Q_thrpt = integrate.cumtrapz(np.abs(I), time_in_sec, initial=0) / 3600

#%%
plt.figure(figsize=(8, 6))
# markerfacecolor
plt.plot(time_in_sec, Q_thrpt, '-', mfc='none')

#%%
plt.figure(figsize=(8, 6))
# markerfacecolor
plt.plot(Q_thrpt, Q_thrpt, '-', mfc='none')