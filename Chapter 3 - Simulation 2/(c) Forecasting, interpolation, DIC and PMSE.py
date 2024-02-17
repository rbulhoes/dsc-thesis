########################### Interpolation and forecasting

e_Beta = [np.zeros((T + 1, p, q)) for k in range(0, smp_size)]
aux_kron = nb_kp(((np.eye(q).T).reshape((q*q, 1))).T, np.eye(p))
for k in range(0, int(smp_size)):
  for t in range(1, T + 1):
    e_Beta[k][t] = np.matmul(aux_kron, nb_kp(np.eye(q), e_beta[k][t]))

e_Y = [np.zeros((T + 1, N, q)) for k in range(0, smp_size)]
for k in range(0, int(smp_size)):
  for t in range(1, T + 1):
    e_Y[k][t] = (e_y[k][t].reshape(q, N)).T
 
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
    avg1 = (X_i[t] @ e_Beta[k][t])
    avg2 = (aux_B_k @ (e_Y[k][t]-(X[t] @ e_Beta[k][t])))
    avg_t_k = avg1 + avg2 
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
  aux_tot, square_dif_sum = 0, 0
  for n in np.arange(0, N_i):
    for i in np.arange(0, q):
      for t in np.arange(1, T + 1):
        if np.isnan(Y_i[t][n,i]):
          pass
        else:
          aux_tot += 1
          aux_sum = 0
          for k in np.arange(0, K):
            aux_sum += e_Y_i[k][t][n,i]
          aux_mean = aux_sum / K
          square_dif_sum += pow(M_Y_i[t][n,i] - aux_mean, 2)
  square_dif_mean = square_dif_sum / aux_tot
  return square_dif_mean
PMSE_value = PMSE(M_Y_i = Y_i)

def DIC(M_y: np.ndarray, M_x: np.ndarray, e_beta = e_beta, e_V = e_V,
        e_phi = e_phi, e_D = e_D, e_Sigma = e_Sigma, e_Beta = e_Beta,
        num_obs = N_o, M_P = P, Id_obs = Id_o, e_y = e_y) -> np.float64:
  K = np.shape(e_Beta)[0]
  T = np.shape(e_Beta)[1] - 1
  p = np.shape(e_Beta)[2]
  q = np.shape(e_Beta)[3]
  N = np.shape(e_D)[2]
  L = np.zeros((K, T))
  for k in range(0, K):
    cov_k = e_V[k] * nb_kp(e_Sigma[k],
                           nb_BR(matrix = e_D[k],
                                 need_transp = True,
                                 scalar = e_phi[k]))
    for t in range(1, T + 1):
      if num_obs[t] == N*q:
        mu_t = np.matmul(M_x[t], e_beta[k][t]).reshape(N*q)
        L[k,t-1] = multivariate_normal.logpdf(x = M_y[t].reshape(N*q),
                                              mean = mu_t,
                                              cov = cov_k)
      elif num_obs[t] > 0 and num_obs[t] < N*q:
        mu_t = M_P[t] @ (M_x[t] @ e_beta[k][t])
        mu_o_t = mu_t[0:num_obs[t]]
        mu_m_t = mu_t[num_obs[t]:(N*q)]
        Delta_t = (M_P[t] @ cov_k) @ M_P[t].T
        Delta_oo_t = Delta_t[0:num_obs[t],0:num_obs[t]]
        Delta_om_t = Delta_t[0:num_obs[t],num_obs[t]:(N*q)]
        Delta_mo_t = Delta_om_t.T
        Delta_mm_t = Delta_t[num_obs[t]:(N*q),
                             num_obs[t]:(N*q)]
        inv_Delta_mm_t = np.linalg.inv(Delta_mm_t)
        y_obs_t = (M_y[t][Id_obs[t]]).reshape(num_obs[t])
        y_mis_t = e_y[k][t][np.logical_not(Id_obs[t])]
        dif_m_t = y_mis_t - mu_m_t
        aux_mult = Delta_om_t @ inv_Delta_mm_t
        mu_cond_t = (mu_o_t + (aux_mult @ dif_m_t)).reshape(num_obs[t])
        Delta_cond_t = Delta_oo_t - (aux_mult @ Delta_mo_t)
        L[k,t-1] = multivariate_normal.logpdf(x = y_obs_t,
                                              mean = mu_cond_t,
                                              cov = Delta_cond_t)
      else:
        pass
  beta_sum = np.zeros((T + 1, p*q, 1))
  V_sum = 0
  phi_sum = 0
  D_sum = np.zeros((2, N))
  Sigma_sum = np.zeros((q, q))
  y_sum = np.zeros((T + 1, N*q, 1))
  for k in range(0, K):
    V_sum += e_V[k]
    phi_sum += e_phi[k]
    Sigma_sum += e_Sigma[k]
    D_sum += e_D[k]
    for t in range(0, T + 1):
      beta_sum[t] += e_beta[k][t]
      y_sum[t] += e_y[k][t]
  y_bar = y_sum / K
  beta_bar = beta_sum / K
  V_bar = V_sum / K
  phi_bar = phi_sum / K
  D_bar = D_sum / K
  Sigma_bar = Sigma_sum / K
  cov_bar = V_bar * nb_kp(Sigma_bar,
                          nb_BR(matrix = D_bar,
                                need_transp = True,
                                scalar = phi_bar))
  F = np.zeros(T)
  for t in range(1, T + 1):
    if num_obs[t] == N*q:
      mu_bar = (np.matmul(M_x[t], beta_bar[t])).reshape(N*q)
      F[t-1] = multivariate_normal.logpdf(x = M_y[t].reshape(N*q),
                                          mean = mu_bar,
                                          cov = cov_bar)
    elif num_obs[t] > 0 and num_obs[t] < N*q:
      mu_t = M_P[t] @ (M_x[t] @ beta_bar[t])
      mu_o_t = mu_t[0:num_obs[t]]
      mu_m_t = mu_t[num_obs[t]:(N*q)]
      Delta_t = (M_P[t] @ cov_bar) @ M_P[t].T
      Delta_oo_t = Delta_t[0:num_obs[t],0:num_obs[t]]
      Delta_om_t = Delta_t[0:num_obs[t],num_obs[t]:(N*q)]
      Delta_mo_t = Delta_om_t.T
      Delta_mm_t = Delta_t[num_obs[t]:(N*q),
                           num_obs[t]:(N*q)]
      inv_Delta_mm_t = np.linalg.inv(Delta_mm_t)
      y_obs_t = (M_y[t][Id_obs[t]]).reshape(num_obs[t])
      y_mis_t = y_bar[t][np.logical_not(Id_obs[t])]
      dif_m_t = y_mis_t - mu_m_t
      aux_mult = Delta_om_t @ inv_Delta_mm_t
      mu_cond_t = (mu_o_t + (aux_mult @ dif_m_t)).reshape(num_obs[t])
      Delta_cond_t = Delta_oo_t - (aux_mult @ Delta_mo_t)
      F[t-1] = multivariate_normal.logpdf(x = y_obs_t,
                                          mean = mu_cond_t,
                                          cov = Delta_cond_t)
    else:
      pass
  return -(4/K)*np.sum(L) + 2*np.sum(F)
DIC_value = DIC(M_y = y, M_x = x)
