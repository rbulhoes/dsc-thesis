############################## True values for the response matrix

###### Exporting Y
time = np.zeros((N, 1))
site = np.arange(1, N + 1).reshape(N, 1)
for t in range(1, T + 1):
  now_time = t*np.ones((N, 1))
  now_site = np.arange(1, N + 1).reshape(N, 1)
  time = np.vstack([time, now_time])
  site = np.vstack([site, now_site])
export_Y = np.hstack([time, site, Y.reshape(N*(T + 1), q)])
np.savetxt("D:/Exportações/Y.txt",
           export_Y,
           delimiter = ",")

###### Exporting Y_i
time_i = np.zeros((N_i, 1))
site_i = np.arange(1, N_i + 1).reshape(N_i, 1)
for t in range(1, T + 1):
  now_time = t*np.ones((N_i, 1))
  now_site = np.arange(1, N_i + 1).reshape(N_i, 1)
  time_i = np.vstack([time_i, now_time])
  site_i = np.vstack([site_i, now_site])
export_Y_i = np.hstack([time_i, site_i, Y_i.reshape(N_i*(T + 1), q)])
np.savetxt("D:/Exportações/Y_i.txt",
           export_Y_i,
           delimiter = ",")

###### Exporting Y_p
time_p = np.ones((N, 1))
site_p = np.arange(1, N + 1).reshape(N, 1)
for t in range(2, T_p + 1):
  now_time = t*np.ones((N, 1))
  now_site = np.arange(1, N + 1).reshape(N, 1)
  time_p = np.vstack([time_p, now_time])
  site_p = np.vstack([site_p, now_site])
export_Y_p = np.hstack([time_p, site_p, Y_p.reshape(T_p*N, q)])
np.savetxt("D:/Exportações/Y_p.txt",
           export_Y_p,
           delimiter = ",")

###### Exporting Y_p_i
time_p_i = np.ones((N_i, 1))
site_p_i = np.arange(1, N_i + 1).reshape(N_i, 1)
for t in range(2, T_p + 1):
  now_time = t*np.ones((N_i, 1))
  now_site = np.arange(1, N_i + 1).reshape(N_i, 1)
  time_p_i = np.vstack([time_p_i, now_time])
  site_p_i = np.vstack([site_p_i, now_site])
export_Y_p_i = np.hstack([time_p_i, site_p_i, Y_p_i.reshape(T_p*N_i, q)])
np.savetxt("D:/Exportações/Y_p_i.txt",
           export_Y_p_i,
           delimiter = ",")

###### Exporting S
np.savetxt("D:/Exportações/S.txt",
           np.hstack([np.arange(1, N + 1).reshape(N, 1), S.T]),
           delimiter = ",")

###### Exporting S_i
np.savetxt("D:/Exportações/S_i.txt",
           np.hstack([np.arange(1, N_i + 1).reshape(N_i, 1), S_i.T]),
           delimiter = ",")

###### Exporting D (true matrix)
np.savetxt("D:/Exportações/D.txt",
           np.hstack([np.arange(1, N + 1).reshape(N, 1), D.T]),
           delimiter = ",")

###### Exporting D_i (true matrix)
np.savetxt("D:/Exportações/D_i.txt",
           np.hstack([np.arange(1, N_i + 1).reshape(N_i, 1), D_i.T]),
           delimiter = ",")

###### Exporting Beta (true matrix)
np.savetxt("D:/Exportações/Beta.txt",
           np.hstack([np.arange(0, T + 1).reshape(T + 1, 1),
                      Beta[:,0,0].reshape(T + 1, 1),
                      Beta[:,1,0].reshape(T + 1, 1),
                      Beta[:,0,1].reshape(T + 1, 1),
                      Beta[:,1,1].reshape(T + 1, 1)]),
           delimiter = ",")

###### Exporting Beta_p (true matrix)
np.savetxt("D:/Exportações/Beta_p.txt",
           np.hstack([np.arange(1, T_p + 1).reshape(T_p, 1),
                      Beta_p[:,0,0].reshape(T_p, 1),
                      Beta_p[:,1,0].reshape(T_p, 1),
                      Beta_p[:,0,1].reshape(T_p, 1),
                      Beta_p[:,1,1].reshape(T_p, 1)]),
           delimiter = ",")



############################## Estimated values 
############################## (depending on the model and/or the algorithm)

###### Exporting e_Y
# Use the quantities time and site previously defined
export_e_Y = np.hstack([1*np.ones((N*(T + 1), 1)),
                        time, site,
                        e_Y[0].reshape(N*(T + 1), q)])
for j in range(1, int(smp_size)):
  now = np.hstack([(j+1)*np.ones((N*(T + 1), 1)),
                   time, site,
                   e_Y[j].reshape(N*(T + 1), q)])
  export_e_Y = np.vstack([export_e_Y, now])
np.savetxt("D:/Exportações/e_Y.txt",
           export_e_Y,
           delimiter = ",")

###### Exporting e_Y_i
# Use the quantities time_i and site_i previously defined
export_e_Y_i = np.hstack([1*np.ones((N_i*(T + 1), 1)),
                          time_i, site_i,
                          e_Y_i[0].reshape(N_i*(T + 1), q)])
for j in range(1, int(smp_size)):
  now = np.hstack([(j+1)*np.ones((N_i*(T + 1), 1)),
                   time_i, site_i,
                   e_Y_i[j].reshape(N_i*(T + 1), q)])
  export_e_Y_i = np.vstack([export_e_Y_i, now])
np.savetxt("D:/Exportações/e_Y_i.txt",
           export_e_Y_i,
           delimiter = ",")

###### Exporting e_Y_p
export_e_Y_p = np.hstack([1*np.ones((N*T_p, 1)),
                          time_p, site_p,
                          e_Y_p[0].reshape(N*T_p, q)])
for j in range(1, int(smp_size)):
  now = np.hstack([(j+1)*np.ones((N*T_p, 1)),
                   time_p, site_p,
                   e_Y_p[j].reshape(N*T_p, q)])
  export_e_Y_p = np.vstack([export_e_Y_p, now])
np.savetxt("D:/Exportações/e_Y_p.txt",
           export_e_Y_p,
           delimiter = ",")

###### Exporting e_Y_p_i
export_e_Y_p_i = np.hstack([1*np.ones((N_i*T_p, 1)),
                            time_p_i, site_p_i,
                            e_Y_p_i[0].reshape(N_i*T_p, q)])
for j in range(1, int(smp_size)):
  now = np.hstack([(j+1)*np.ones((N_i*T_p, 1)),
                   time_p_i, site_p_i,
                   e_Y_p_i[j].reshape(N_i*T_p, q)])
  export_e_Y_p_i = np.vstack([export_e_Y_p_i, now])
np.savetxt("D:/Exportações/e_Y_p_i.txt",
           export_e_Y_p_i,
           delimiter = ",")

###### Exporting e_D
site = np.arange(1, N + 1).reshape(N, 1)
export_e_D = np.hstack([1*np.ones(N).reshape(N, 1), site,
                        (e_D[0][0,:].T).reshape(N, 1),
                        (e_D[0][1,:].T).reshape(N, 1)])
for j in range(1, int(smp_size)):
  now = np.hstack([(j+1)*np.ones(N).reshape(N, 1), site,
                   (e_D[j][0,:].T).reshape(N, 1),
                   (e_D[j][1,:].T).reshape(N, 1)])
  export_e_D = np.vstack([export_e_D, now])
np.savetxt("D:/Exportações/e_D.txt",
           export_e_D,
           delimiter = ",")

###### Exporting e_D_i
site_i = np.arange(1, N_i + 1).reshape(N_i, 1)
export_e_D_i = np.hstack([1*np.ones(N_i).reshape(N_i, 1), site_i,
                          (e_D_i[0][0,:].T).reshape(N_i, 1),
                          (e_D_i[0][1,:].T).reshape(N_i, 1)])
for j in range(1, int(smp_size)):
  now = np.hstack([(j+1)*np.ones(N_i).reshape(N_i, 1), site_i,
                   (e_D_i[j][0,:].T).reshape(N_i, 1),
                   (e_D_i[j][1,:].T).reshape(N_i, 1)])
  export_e_D_i = np.vstack([export_e_D_i, now])
np.savetxt("D:/Exportações/e_D_i.txt",
           export_e_D_i,
           delimiter = ",")

###### Exporting e_phi
export_e_phi = e_phi[0]
for j in range(1, int(smp_size)):
  now = e_phi[j]
  export_e_phi = np.concatenate((export_e_phi, now), axis = None)
np.savetxt("D:/Exportações/e_phi.txt",
           np.hstack([np.arange(1,
                                int(smp_size)+1).reshape(int(smp_size), 1),
                      export_e_phi.reshape(int(smp_size), 1)]),
           delimiter = ",")

###### Exporting e_V * e_Sigma
export_e_VSigma = np.hstack([e_V[0]*e_Sigma[0][0,0],
                             e_V[0]*e_Sigma[0][0,1],
                             e_V[0]*e_Sigma[0][1,1]])
for j in range(1, int(smp_size)):
  now = np.hstack([e_V[j]*e_Sigma[j][0,0],
                   e_V[j]*e_Sigma[j][0,1],
                   e_V[j]*e_Sigma[j][1,1]])
  export_e_VSigma = np.vstack([export_e_VSigma, now])
np.savetxt("D:/Exportações/e_VSigma.txt",
           np.hstack([np.arange(1,
                                int(smp_size)+1).reshape(int(smp_size), 1),
                      export_e_VSigma]),
           delimiter = ",")

###### Exporting e_Beta
export_e_Beta = np.hstack([1*np.ones(T + 1).reshape(T + 1, 1),
                           np.arange(0, T + 1).reshape(T + 1, 1),
                           e_Beta[0][:,0,0].reshape(T + 1, 1),
                           e_Beta[0][:,1,0].reshape(T + 1, 1),
                           e_Beta[0][:,0,1].reshape(T + 1, 1),
                           e_Beta[0][:,1,1].reshape(T + 1, 1)])
for j in range(1, int(smp_size)):
  now = np.hstack([(j+1)*np.ones(T + 1).reshape(T + 1, 1),
                   np.arange(0, T + 1).reshape(T + 1, 1),
                   e_Beta[j][:,0,0].reshape(T + 1, 1),
                   e_Beta[j][:,1,0].reshape(T + 1, 1),
                   e_Beta[j][:,0,1].reshape(T + 1, 1),
                   e_Beta[j][:,1,1].reshape(T + 1, 1)])
  export_e_Beta = np.vstack([export_e_Beta, now])
np.savetxt("D:/Exportações/e_Beta.txt",
           export_e_Beta,
           delimiter = ",")

###### Exporting e_Beta_p
export_e_Beta_p = np.hstack([1*np.ones(T_p).reshape(T_p, 1),
                             np.arange(1, T_p + 1).reshape(T_p, 1),
                             e_Beta_p[0][:,0,0].reshape(T_p, 1),
                             e_Beta_p[0][:,1,0].reshape(T_p, 1),
                             e_Beta_p[0][:,0,1].reshape(T_p, 1),
                             e_Beta_p[0][:,1,1].reshape(T_p, 1)])
for j in range(1, int(smp_size)):
  now = np.hstack([(j+1)*np.ones(T_p).reshape(T_p, 1),
                   np.arange(1, T_p + 1).reshape(T_p, 1),
                   e_Beta_p[j][:,0,0].reshape(T_p, 1),
                   e_Beta_p[j][:,1,0].reshape(T_p, 1),
                   e_Beta_p[j][:,0,1].reshape(T_p, 1),
                   e_Beta_p[j][:,1,1].reshape(T_p, 1)])
  export_e_Beta_p = np.vstack([export_e_Beta_p, now])
np.savetxt("D:/Exportações/e_Beta_p.txt",
           export_e_Beta_p,
           delimiter = ",")
