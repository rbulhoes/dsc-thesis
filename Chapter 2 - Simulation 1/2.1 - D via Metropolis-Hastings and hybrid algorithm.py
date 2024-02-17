@njit
def updt_D(it: int64, S_obs: float64[:,:], M_X: float64[:,:,:],
           M_Y: float64[:,:,:], Beta_usage: float64[:,:,:], V_usage: float64,
           phi_usage: float64, sigma2d_usage: float64[:,:],
           D_usage: float64[:,:], R_usage: float64[:,:],
           Sigma_usage: float64[:,:], warm_up: int64, acc_now: int64,
           delta_usage: float64[:,:], size = 50, rate = 0.44,
           cte = 0.01) -> (float64[:,:], float64[:,:]):
  D_prev, delta_prev = D_usage.copy(), delta_usage.copy()
  N_obs = M_Y.shape[1]
  for i in range(0, 2):
    for j in range(2, N_obs):
      D_prop = D_prev.copy()
      delta_now = delta_prev.copy()
      np.random.seed(it + i + j)
      D_prop_ij = np.random.normal(loc = D_prev[i,j],
                                   scale = delta_now[i,j])
      D_prop[i,j] = D_prop_ij
      lcd_D_prop = log_krnl_D(par = D_prop,
                              S_obs = S_obs,
                              M_X = M_X,
                              M_Y = M_Y,
                              Sigma_usage = Sigma_usage,
                              Beta_usage = Beta_usage,
                              V_usage = V_usage,
                              phi_usage = phi_usage,
                              sigma2d_usage = sigma2d_usage,
                              R_usage = R_usage)
      lcd_D_prev = log_krnl_D(par = D_prev,
                              S_obs = S_obs,
                              M_X = M_X,
                              M_Y = M_Y,
                              Sigma_usage = Sigma_usage,
                              Beta_usage = Beta_usage,
                              V_usage = V_usage,
                              phi_usage = phi_usage,
                              sigma2d_usage = sigma2d_usage,
                              R_usage = R_usage)
      log_rho = lcd_D_prop - lcd_D_prev
      log_alpha = min([0, log_rho])
      log_u = log_rUnif(it + i + j)
      
      seq = [j for j in range(0, warm_up + 1, size)]
      if it in seq:
        frac = acc_now[i,j] / it
        if frac < rate:
          delta_now_ij = exp(log(delta_now[i,j]) - min(cte, 1/sqrt(it)))
          delta_now[i,j] = delta_now_ij
        else:
          delta_now_ij = exp(log(delta_now[i,j]) + min(cte, 1/sqrt(it)))
          delta_now[i,j] = delta_now_ij
        if isfinite(log_alpha) and log_u <= log_alpha:
          D_now = D_prop
        else:
          D_now = D_prev
      else:
        if isfinite(log_alpha) and log_u <= log_alpha:
          D_now = D_prop
        else:
          D_now = D_prev
      D_prev = D_now.copy()
      delta_prev = delta_now.copy()
        
  return D_now, delta_now



########################### Creating objects to store the MCMC samples
burn_in = 1000
num_iter = 10000 + burn_in
thin = 10
smp_size = int((num_iter - burn_in)/thin)

e_phi = [0 for k in range(0, smp_size)]
e_V = [0 for k in range(0, smp_size)]
e_Beta = [np.zeros((T + 1, p, q)) for k in range(0, smp_size)]
e_Sigma = [np.zeros((q, q)) for k in range(0, smp_size)]
e_D = [np.zeros((2, N)) for k in range(0, smp_size)]



########################### Initial values for MCMC estimation
phi_prev = 1.0
V_prev = 2.0
Sigma_prev = 2.0*np.identity(q)
Beta_prev = np.zeros((T + 1, p, q))
D_prev = np.copy(S)
delta_phi = 12.0
delta_D = np.hstack([np.zeros((2,2)), 0.05*np.ones((2, N-2))]) 



########################### Metropolis-within-Gibbs algorithm
# ind will vary between 0 and smp_size - 1 (i.e. k = 1, ..., K)
ind = 0 
acc_phi, acc_D = 0, np.zeros((2, N))

# Algorithm and processing time in seconds
start = timer()
for j in range(1, num_iter + 1):
   
#  Beta_curr = np.copy(Beta)
  Beta_curr = updt_Beta(it = seed_v + j,
                        init_M0 = M0,
                        init_C0 = C0,
                        M_Y = Y,
                        M_X = X,
                        M_G = G,
                        D_usage = D_prev,
                        W_usage = W,
                        phi_usage = phi_prev,
                        V_usage = V_prev,
                        Sigma_usage = Sigma_prev)

#  V_curr = np.copy(V).item()
  V_curr = updt_V(it = seed_v + j,
                  init_a_V = a_V,
                  init_b_V = b_V,
                  init_M0 = M0,
                  init_C0 = C0,
                  M_Y = Y,
                  M_X = X,
                  M_G = G,
                  D_usage = D_prev,
                  Beta_usage = Beta_curr,
                  phi_usage = phi_prev,
                  Sigma_usage = Sigma_prev,
                  W_usage = W)

#  phi_curr = np.copy(phi).item()
  phi_curr = updt_phi(it = seed_v + j,
                      init_a_phi = a_phi,
                      init_b_phi = b_phi,
                      D_usage = D_prev,
                      M_Y = Y,
                      M_X = X,
                      Beta_usage = Beta_curr,
                      V_usage = V_curr,
                      phi_usage = phi_prev,
                      Sigma_usage = Sigma_prev,
                      delta2 = delta_phi)
  
#  Sigma_curr = np.copy(Sigma)
  Sigma_curr = updt_Sigma(it = seed_v + j,
                          init_a_Sigma = a_Sigma,
                          init_b_Sigma = b_Sigma,
                          init_M0 = M0,
                          init_C0 = C0,
                          M_Y = Y,
                          M_X = X,
                          M_G = G,
                          Beta_usage = Beta_curr,
                          W_usage = W,
                          phi_usage = phi_curr,
                          V_usage = V_curr,
                          D_usage = D_prev)
  
#  D_curr = np.copy(D)
  D_curr, delta_D = updt_D(it = seed_v + j,
                           S_obs = S,
                           M_X = X,
                           M_Y = Y,
                           Sigma_usage = Sigma_curr,
                           Beta_usage = Beta_curr,
                           V_usage = V_curr,
                           phi_usage = phi_curr,
                           sigma2d_usage = sigma2d_usg,
                           R_usage = R_usg,
                           D_usage = D_prev,
                           delta_usage = delta_D,
                           warm_up = burn_in,
                           acc_now = acc_D)
  
  if abs(phi_curr - phi_prev) != 0:
    acc_phi += 1
  else:
    pass
  
  for r in range(0, 2):
    for n in range(0, N):
      if abs(D_curr[r,n] - D_prev[r,n]) != 0:
        acc_D[r,n] += 1
      else:
        pass

  if j > burn_in and j % thin == 0:
    e_V[ind] = V_curr
    e_Beta[ind] = Beta_curr
    e_phi[ind] = phi_curr
    e_Sigma[ind] = Sigma_curr
    e_D[ind] = D_curr
    ind += 1
  else:
    pass

  V_prev = V_curr
  phi_prev = phi_curr
  Beta_prev = Beta_curr
  Sigma_prev = Sigma_curr
  D_prev = D_curr

  print(j)
end = timer()
print(end - start)
