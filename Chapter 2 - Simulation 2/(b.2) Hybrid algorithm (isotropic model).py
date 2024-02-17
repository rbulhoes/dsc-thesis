########################### Creating objects to store the MCMC samples
burn_in = 5000
num_iter = 15000 + burn_in
thin = 15
smp_size = int((num_iter - burn_in)/thin)

e_phi = [0 for k in range(0, smp_size)]
e_V = [0 for k in range(0, smp_size)]
e_Beta = [np.zeros((T + 1, p, q)) for k in range(0, smp_size)]
e_Sigma = [np.zeros((q, q)) for k in range(0, smp_size)]
e_D = [np.zeros((2, N)) for k in range(0, smp_size)]



########################### Initial values for MCMC estimation
phi_prev = 2.0
V_prev = 2.0
Sigma_prev = 2.0*np.identity(q)
Beta_prev = np.zeros((T + 1, p, q))
D_prev = np.copy(S)
delta_phi = 130.0



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
  
  D_curr = np.copy(S)
  
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
