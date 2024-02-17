@njit
def updt_D(it: int64, S_obs: float64[:,:], M_X: float64[:,:,:],
           M_Y: float64[:,:,:], Beta_usage: float64[:,:,:], V_usage: float64,
           phi_usage: float64, sigma2d_usage: float64[:,:],
           D_usage: float64[:,:], R_usage: float64[:,:],
           Sigma_usage: float64[:,:]) -> float64[:,:]:
  N_obs = M_Y.shape[1]
  D_prev = D_usage.copy()
  w = 1
  m = +inf
  lower = -inf
  upper = +inf
  for i in range(0, 2):
    for j in range(2, N_obs):
      np.random.seed(it + i + j)
      D_prop = D_prev.copy()

      x0 = D_prop[i,j]
      x0_D = D_prop.copy()

      # Find the log density at the initial point, if not already known.
      gx0 = log_krnl_D(par = x0_D,
                       S_obs = S_obs,
                       M_X = M_X,
                       M_Y = M_Y,
                       Sigma_usage = Sigma_usage,
                       Beta_usage = Beta_usage,
                       V_usage = V_usage,
                       phi_usage = phi_usage,
                       sigma2d_usage = sigma2d_usage,
                       R_usage = R_usage)
      # Determine the slice level, in log terms.
      logy = gx0 - np.random.exponential(scale = 1.0,
                                         size = 1).item()
  
      # Find the initial interval to sample from.
      u = np.random.uniform(0, w)
      L = x0 - u
      R = x0 + (w - u) # should guarantee that x0 is in [L,R], even with roundoff
      
      L_D = D_prop.copy()
      L_D[i,j] = L
      R_D = D_prop.copy()
      R_D[i,j] = R

      # Expand the interval until its ends are outside the slice, or until
      # the limit on steps is reached.
      if isfinite(m) is False:
        while L > lower and logy < log_krnl_D(par = L_D,
                                              S_obs = S_obs,
                                              M_X = M_X,
                                              M_Y = M_Y,
                                              Sigma_usage = Sigma_usage,
                                              Beta_usage = Beta_usage,
                                              V_usage = V_usage,
                                              phi_usage = phi_usage,
                                              sigma2d_usage = sigma2d_usage,
                                              R_usage = R_usage):
          L = L - w
          L_D[i,j] = L
        while R < upper and logy < log_krnl_D(par = R_D,
                                              S_obs = S_obs,
                                              M_X = M_X,
                                              M_Y = M_Y,
                                              Sigma_usage = Sigma_usage,
                                              Beta_usage = Beta_usage,
                                              V_usage = V_usage,
                                              phi_usage = phi_usage,
                                              sigma2d_usage = sigma2d_usage,
                                              R_usage = R_usage):
          R = R + w
          R_D[i,j] = R
      else:
        pass
  
      # Shrink interval to lower and upper bounds.
      if L < lower:
        L = lower
        L_D[i,j] = L
      else:
        pass
      if R > upper:
        R = upper
        R_D[i,j] = R
      else:
        pass
  
      # Sample from the interval, shrinking it on each rejection.
      x1_D = D_prop.copy()
      x1 = np.random.uniform(L, R)
      x1_D[i,j] = x1
      gx1 = log_krnl_D(par = x1_D,
                       S_obs = S_obs,
                       M_X = M_X,
                       M_Y = M_Y,
                       Sigma_usage = Sigma_usage,
                       Beta_usage = Beta_usage,
                       V_usage = V_usage,
                       phi_usage = phi_usage,
                       sigma2d_usage = sigma2d_usage,
                       R_usage = R_usage)
  
      while gx1 < logy:
        x1 = np.random.uniform(L, R)
        x1_D[i,j] = x1
        gx1 = log_krnl_D(par = x1_D,
                         S_obs = S_obs,
                         M_X = M_X,
                         M_Y = M_Y,
                         Sigma_usage = Sigma_usage,
                         Beta_usage = Beta_usage,
                         V_usage = V_usage,
                         phi_usage = phi_usage,
                         sigma2d_usage = sigma2d_usage,
                         R_usage = R_usage)
        if x1 > x0:
          R = x1
        else:
          L = x1
  
      D_prev = x1_D.copy()
        
  return x1_D



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
V_prev = 1.0
Sigma_prev = 2.0*np.identity(q)
Beta_prev = np.zeros((T + 1, p, q))
D_prev = np.copy(S)
delta_phi = 100.0 # 2.0 if T = 10, 12.0 if T = 100, and 100.0 if T = 1000 



########################### Hybrid algorithm
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
  D_curr = updt_D(it = seed_v + j,
                  S_obs = S,
                  M_X = X,
                  M_Y = Y,
                  Sigma_usage = Sigma_curr,
                  Beta_usage = Beta_curr,
                  V_usage = V_curr,
                  phi_usage = phi_curr,
                  sigma2d_usage = sigma2d_usg,
                  R_usage = R_usg,
                  D_usage = D_prev)
  
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
