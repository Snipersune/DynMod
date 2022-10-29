
#%%
import numpy as np
from matplotlib import pyplot as plt
alpha = 0.9
beta = 2
dt = 0.02
tsteps = 12960000
n_variants = 100

print("Start loading data")
q_vars = np.loadtxt("data/q_vars.txt")
Np_ts = np.loadtxt("data/Np_ts.txt")
end_occs = np.loadtxt("data/end_occs.txt")
print("Done")

#%%
codon_length = len(q_vars[0,:])

tsteps_transient = int(tsteps/10)
Np_t_noTransient = Np_ts[:,tsteps_transient:]

mean = np.mean(Np_t_noTransient, axis=1)
variance = np.var(Np_t_noTransient, axis=1)
ff = variance/mean

end_occs = end_occs[:, tsteps_transient:]
pL = np.mean(end_occs, axis=1)

print("Fano factor:", ff)

J_all = beta*pL
J_orig = J_all[0]


idx_sorted_J = np.argsort(J_all)
idx_top25 = idx_sorted_J[int((n_variants+1)*3/4):]
idx_bot25 = idx_sorted_J[0:int((n_variants+1)*3/4)]

top25_q = q_vars[idx_top25,:]
bot25_q = q_vars[idx_bot25,:]

top25_qi_mean = np.mean(top25_q, axis=0)
bot25_qi_mean = np.mean(bot25_q, axis=0)

q_mean = np.mean(q_vars, axis=1) 

idx_sorted_q = np.argsort(q_mean)
q_mean_sorted = q_mean[idx_sorted_q]
J_all_sorted_q = J_all[idx_sorted_q]

#%%
print("Begin plotting")
plt.figure()
plt.plot(np.arange(0,tsteps*dt, dt), Np_ts[0,:], label="orig")
plt.plot(np.arange(0,tsteps*dt, dt), Np_ts[1,:])
plt.plot(np.arange(0,tsteps*dt, dt), Np_ts[2,:])
plt.plot(np.arange(0,tsteps*dt, dt), Np_ts[3,:])
plt.plot(np.arange(0,tsteps*dt, dt), Np_ts[4,:])
plt.legend()

plt.figure()
plt.plot(range(1,codon_length+1), q_vars[0,:], label="orig")
plt.plot(range(1,codon_length+1), top25_qi_mean, label="top25")
plt.plot(range(1,codon_length+1), bot25_qi_mean, label="bot25")
plt.legend()

plt.figure()
plt.plot(q_mean_sorted, J_all_sorted_q)

plt.figure()
plt.hist([J_all[1:]])
plt.vlines(J_all[0], ymin=0, ymax=30, colors="red")

plt.show()
# %%
