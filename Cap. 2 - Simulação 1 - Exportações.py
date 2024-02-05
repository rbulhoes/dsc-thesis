# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:51:23 2021

@author: rodri
"""

############################## True values for the response matrix

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

