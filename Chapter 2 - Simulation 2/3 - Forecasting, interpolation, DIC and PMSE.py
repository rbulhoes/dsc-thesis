########################### Interpolation and forecasting

def DIC(M_Y: np.ndarray, M_X: np.ndarray, e_Beta = e_Beta, e_V = e_V,
        e_phi = e_phi, e_D = e_D, e_Sigma = e_Sigma, df = 1) -> np.float64:
  K = np.shape(e_Beta)[0]
  T = np.shape(e_Beta)[1] - 1
  p = np.shape(e_Beta)[2]
  q = np.shape(e_Beta)[3]
  N = np.shape(e_D)[2]
  L = np.zeros((K, T))
  for k in range(0, K):
    leftcov_k = e_V[k] * nb_BR(matrix = e_D[k],
                               need_transp = True,
                               scalar = e_phi[k])
    for t in range(1, T + 1):
      L[k,t-1] = matrix_normal.logpdf(X = M_Y[t],
                                      mean = np.matmul(M_X[t],
                                                       e_Beta[k][t]),
                                      rowcov = leftcov_k,
                                      colcov = e_Sigma[k])
  Beta_sum = np.zeros((T + 1, p, q))
  V_sum = 0
  phi_sum = 0
  D_sum = np.zeros((2, N))
  Sigma_sum = np.zeros((q, q))
  for k in range(0, K):
    V_sum += e_V[k]
    phi_sum += e_phi[k]
    Sigma_sum += e_Sigma[k]
    D_sum += e_D[k]
    for t in range(0, T + 1):
      Beta_sum[t] += e_Beta[k][t]
  Beta_bar = Beta_sum / K
  V_bar = V_sum / K
  phi_bar = phi_sum / K
  D_bar = D_sum / K
  Sigma_bar = Sigma_sum / K
  leftcov_bar = V_bar * nb_BR(matrix = D_bar,
                              need_transp = True,
                              scalar = phi_bar)
  F = np.zeros(T)
  for t in range(1, T + 1):
    F[t-1] += matrix_normal.logpdf(X = M_Y[t],
                                   mean = np.matmul(M_X[t],
                                                    Beta_bar[t]),
                                   rowcov = leftcov_bar,
                                   colcov = Sigma_bar)
  return -(4/K)*np.sum(L) + 2*np.sum(F)
DIC_value = DIC(M_Y = Y, M_X = X)

e_Y = [np.zeros((T + 1, N, q)) for j in range(0, smp_size)]
for k in range(0, int(smp_size)):
  for t in range(1, T + 1):
    e_Y[k][t] = Y[t]

e_Beta_p = [np.zeros((T_p, p, q)) for k in range(0, smp_size)]
for k in range(0, int(smp_size)):
  for t in range(1, T_p + 1):
    aux_t = G[T + t]
    cov_t = W.copy()
    if t > 1:
      for t2 in reversed(range(2, t + 1)):
        aux_t = np.matmul(aux_t, G[T + t2])
        cov_t += np.matmul(aux_t.T, W @ aux_t)
    else:
      pass
    e_Beta_p[k][t - 1] = rMVN(n = 1,
                              avg = np.matmul(aux_t, e_Beta[k][T]),
                              left = e_V[k] * cov_t,
                              right = e_Sigma[k],
                              n_seed = seed_v + k*t)
    
e_Y_p = [np.zeros((T_p, N, q)) for k in range(0, smp_size)]
for k in range(0, int(smp_size)):
  for t in range(1, T_p + 1):
    e_Y_p[k][t - 1] = rMVN(n = 1,
                           avg = np.matmul(X[T + t], e_Beta_p[k][t - 1]),
                           left = e_V[k] * nb_BR(matrix = e_D[k],
                                                 need_transp = True,
                                                 scalar = e_phi[k]),
                           right = e_Sigma[k],
                           n_seed = seed_v + k*t)

e_D_i = [S_i for k in range(0, smp_size)]

e_Y_i = [np.zeros((T + 1, N_i, q)) for j in range(0, smp_size)]
e_Y_p_i = [np.zeros((T_p, N_i, q)) for j in range(0, smp_size)]
for k in range(0, int(smp_size)):
  B_k = nb_BR(matrix = e_D[k],
              need_transp = True,
              scalar = e_phi[k])
  inv_B_k = np.linalg.inv(B_k)
  dist_oi_k = np.zeros((N, N_i))
  for n in range(0, N):
    for m in range(0, N_i):
      dist_oi_k[n, m] = dist(e_D[k][:,n], e_D_i[k][:,m]) 
  B_oi_k = np.exp(-e_phi[k] * dist_oi_k)
  B_io_k = B_oi_k.T
  B_i_k = nb_BR(matrix = e_D_i[k],
                need_transp = True,
                scalar = e_phi[k])
  aux_B_k = np.matmul(B_io_k, inv_B_k)
  cov_B_k = B_i_k - np.matmul(aux_B_k, B_oi_k)
  for t in range(1, T + 1):
    avg_t_k = (X_i[t] @ e_Beta[k][t]) + (aux_B_k @ (Y[t]-(X[t] @ e_Beta[k][t])))
    e_Y_i[k][t] = rMVN(n = 1,
                       avg = avg_t_k,
                       left = e_V[k] * cov_B_k,
                       right = e_Sigma[k],
                       n_seed = seed_v + k*t)
  for t in range(T + 1, T_tot + 1):
    avg1 = (X_i[t] @ e_Beta_p[k][t - T - 1])
    avg2 = aux_B_k @ (e_Y_p[k][t - T - 1]-(X[t] @ e_Beta_p[k][t - T - 1]))
    avg_sum = avg1 + avg2
    e_Y_p_i[k][t - T - 1] = rMVN(n = 1,
                                 avg = avg_sum,
                                 left = e_V[k] * cov_B_k,
                                 right = e_Sigma[k],
                                 n_seed = seed_v + k*t)

def PMSE(M_Y_i: np.ndarray, e_Y_i = e_Y_i) -> np.float64:
  K = np.shape(e_Y_i)[0]
  T = np.shape(e_Y_i)[1] - 1
  N_i = np.shape(e_Y_i)[2]
  q = np.shape(e_Y_i)[3]
  square_dif_sum = 0
  for n in np.arange(0, N_i):
    for i in np.arange(0, q):
      for t in np.arange(1, T + 1):
        aux_sum = 0
        for k in np.arange(0, K):
          aux_sum += e_Y_i[k][t][n,i]
        aux_mean = aux_sum / K
        square_dif_sum += pow(M_Y_i[t][n,i] - aux_mean, 2)
  square_dif_mean = square_dif_sum / (N_i*q*T)
  return square_dif_mean
PMSE_value = PMSE(M_Y_i = Y_i)
