# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 09:37:54 2023

@author: rodri
"""

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

###### Exporting S
np.savetxt("D:/Exportações/S.txt",
           np.hstack([np.arange(1, N + 1).reshape(N, 1), S.T]),
           delimiter = ",")

###### Exporting D (true matrix)
np.savetxt("D:/Exportações/D.txt",
           np.hstack([np.arange(1, N + 1).reshape(N, 1), D.T]),
           delimiter = ",")

###### Exporting Beta (true matrix)
np.savetxt("D:/Exportações/Beta.txt",
           np.hstack([np.arange(0, T + 1).reshape(T + 1, 1),
                      Beta[:,0,0].reshape(T + 1, 1),
                      Beta[:,1,0].reshape(T + 1, 1),
                      Beta[:,0,1].reshape(T + 1, 1),
                      Beta[:,1,1].reshape(T + 1, 1)]),
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

