# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 23:41:23 2023

@author: rodri
"""

########################### Python libraries

import numpy as np
from scipy.stats import uniform, matrix_normal
from timeit import default_timer as timer
from math import pi, exp, log, sqrt, isfinite, inf, dist
from numba import njit, float64, boolean, uint8, int64, prange
from numba.core.errors import NumbaPerformanceWarning
from warnings import simplefilter
from scipy.spatial import distance, distance_matrix
simplefilter('ignore', category = NumbaPerformanceWarning)



########################### Auxiliary functions

# Obtaining the matrices B or R
@njit(parallel = True)
def nb_BR(matrix: float64[:,:], need_transp: boolean, scalar: float64,
          squared = False) -> float64[:,:]:
  if need_transp is True:
    X = matrix.T
  else:
    X = matrix
  n, m = X.shape
  M = np.ones((n, n))
  for i in prange(0, n - 1):
    for j in prange(i + 1, n):
      dist_ij = 0.0
      for k in prange(0, m):
        dist_ij += pow(X[i, k] - X[j, k], 2)
      if squared is False:
        M[i, j] = exp(-scalar * sqrt(dist_ij)) # This is for B
      else:
        M[i, j] = exp(-scalar * dist_ij) # This is for R
      M[j, i] = M[i, j]
      
  return M

# Kronecker product
@njit
def nb_kp(X: float64[:,:], Y: float64[:,:]) -> float64[:,:]:
  return np.kron(X, Y)

# Sampling from matrix-variate normal distribution (slow, but more robust)
def rMVN(n: np.uint8, avg: np.array, left: np.array,
         right: np.array, n_seed: np.int64) -> np.array:
  p, q = avg.shape
  vec_avg = np.matrix.flatten(avg, 'F')
  kron_cov = nb_kp(right, left)
  
  np.random.seed(seed = n_seed)
  vec_sample = np.random.multivariate_normal(mean = vec_avg,
                                             cov = kron_cov,
                                             size = n,
                                             check_valid = 'ignore')

  return np.reshape(vec_sample, (p, q), 'F')

# Sampling from matrix-variate normal distribution (fast, but less robust)
@njit
def chol_rMVN(n_seed: int64, avg: float64[:,:],
              left: float64[:,:], right: float64[:,:]) -> float64[:,:]:
  p, q = avg.shape
  vec_avg = (avg.T).flatten()
  kron_cov = np.kron(right, left)
  np.random.seed(n_seed)
  gen = np.linalg.cholesky(kron_cov) @ np.random.standard_normal(vec_avg.size)
  vec_sample = vec_avg + gen
  return vec_sample.reshape((q, p)).T

# Sampling from the inverse Wishart distribution
@njit
def chol_rIW(nu: float64, M: float64[:,:], n_seed: uint8) -> float64[:,:]:
  np.random.seed(n_seed)
  dim = M.shape[0]
  inv_M = np.linalg.inv(M)
  chol = np.linalg.cholesky(inv_M)
  df = nu + dim - 1
  foo = np.zeros((dim, dim))
  for i in range(dim):
    for j in range(i + 1):
      if i == j:
        foo[i, j] = sqrt(np.random.chisquare(df - (i + 1) + 1))
      else:
        foo[i, j] = np.random.standard_normal()
  return np.linalg.inv(chol @ (foo @ (foo.T @ chol.T)))

# Sampling from the Inverse-Normal Distribution
@njit
def rIN(n_seed: int64, mu: float64, sigma2: float64) -> float64:
  np.random.seed(n_seed)
  return np.random.wald(mu, sigma2)

# Sampling from the uniform distribution and applying the natural logarithm
@njit
def log_rUnif(n_seed: int64) -> float64:
  np.random.seed(n_seed)
  return log(np.random.rand())

# Logarithm of the PDF of the Inverse-Normal distribution
@njit
def logpdf_IN(x: float64, mean: float64, sigma2: float64) -> float64:
  term1 = 0.5*(log(sigma2) - log(2*pi*pow(x, 3)))
  term2 = -0.5*sigma2*pow(x - mean, 2) / (pow(mean, 2) * x)
  return term1 + term2



########################### Data generation

# Sites in [0, 1]²
S = np.vstack([[1/5, 1/5], [4/5, 4/5], [1/5, 4/5], [4/5, 1/5],
                           [1/5, 2/5], [1/5, 3/5], 
               [2/5, 1/5], [2/5, 2/5], [2/5, 3/5], [2/5, 4/5],
               [3/5, 1/5], [3/5, 2/5], [3/5, 3/5], [3/5, 4/5],
                           [4/5, 2/5], [4/5, 3/5]]).T
S_i = np.array([[0.7, 0.5],
                [0.5, 0.7],
                [0.3, 0.3]]).T
#S_i = np.array([[0.5200000, 0.7900000],
#                [0.2100000, 0.3200000],
#                [0.7900000, 0.5200000]]).T
#S_i = np.array([[0.5434049, 0.2783694],
#                [0.3500000, 0.6900000],
#                [0.7900000, 0.5200000]]).T
S_tot = np.hstack([S, S_i])

# The number of sites and the number of replications (or times)
N, N_i, N_tot = S.shape[1], S_i.shape[1], S_tot.shape[1]
T, T_p = 100, 10
T_tot = T + T_p
p, q = 2, 2 
seed_v = 500

V = 1.0
a_V, b_V = 0.001, 0.001

Sigma = np.array([[1.00, 0.80],
                  [0.80, 1.00]])
a_Sigma, b_Sigma = 0.001, 0.001*np.identity(q)

phi = 0.5
a_phi, b_phi = 0.001, 0.001

### B and B_i are specified by an anisotropic structure
xi_x = 1 # +0.85
xi_y = 3 # 2.75
Lambda1 = np.diag([xi_x, xi_y])
Sh_x = 0.05 # +0.04
Sh_y = 0 # +0.00
Lambda2 = np.array([[1,    Sh_x],
                   [Sh_y,    1]])
eta = -pi/6 # -pi/5
Lambda3 = np.array([[np.cos(eta), -np.sin(eta)],
                    [np.sin(eta),  np.cos(eta)]])
Lambda = Lambda1 @ Lambda2 @ Lambda3
A = Lambda.T @ Lambda
u = np.zeros(2)
D = np.zeros((2, N))
D_i = np.zeros((2, N_i))
for n in range(1, N_tot + 1):
  if n <= 2:
    D[:,n - 1] = S[:,n - 1]
  elif n >= 3 and n <= N:
    D[:,n - 1] = Lambda @ S[:,n - 1] + u
  else:
    D_i[:,N_tot - n - 1] = Lambda @ S_i[:,N_tot - n - 1] + u

B_aug = nb_BR(matrix = np.hstack([D, D_i]),
              need_transp = True, scalar = phi)

M0, C0 = np.zeros((p, q)), 1.0*np.identity(p)
W = 1.0*np.identity(p)
Beta_tot = np.zeros((T_tot + 1, p, q))
Beta, Beta_p = np.zeros((T + 1, p, q)), np.zeros((T_p, p, q))

G = np.zeros((T_tot + 1, p, p))
for t in range(0, T_tot + 1):
  if t == 0:
    G[t] = np.identity(p)
    Beta_tot[t] = rMVN(n = 1,
                        avg = M0,
                        left = V*C0,
                        right = Sigma,
                        n_seed = seed_v + t)
  else:
    G[t] = np.identity(p)
    Beta_tot[t] = rMVN(avg = np.matmul(G[t], Beta_tot[t-1]),
                       left = V*W,
                       right = Sigma,
                       n = 1,
                       n_seed = seed_v + t)

for t in range(0, T_tot + 1):
  if t <= T:
    Beta[t] = Beta_tot[t]
  else:
    Beta_p[t - T - 1] = Beta_tot[t] 

# X
X, X_i = np.zeros((T_tot + 1, N, p)), np.zeros((T_tot + 1, N_i, p))
for t in range(1, T_tot + 1):
  X[t] = np.vstack([np.ones(N),
                    uniform.rvs(loc = 0,
                                scale = 1,
                                size = N,
                                random_state = seed_v + t)]).T
  X_i[t] = np.vstack([np.ones(N_i),
                      uniform.rvs(loc = 0,
                                  scale = 1,
                                  size = N_i,
                                  random_state = seed_v + t)]).T
  
# Response matrices
Y_aug = np.zeros((T_tot + 1, N + N_i, q))
for t in range(1, T_tot + 1):
  Y_aug[t] = rMVN(n = 1,
                  avg = np.matmul(np.vstack([X[t], X_i[t]]), Beta_tot[t]),
                  left = V*B_aug,
                  right = Sigma,
                  n_seed = seed_v + t)

Y, Y_i = np.zeros((T+1,N,q)), np.zeros((T+1,N_i,q))
Y_p, Y_p_i = np.zeros((T_p,N,q)), np.zeros((T_p,N_i,q))
for t in range(1, T_tot + 1):
  if t <= T:
    Y[t] = Y_aug[t][0:N,:]
    Y_i[t] = Y_aug[t][N:N_tot,:]
  else:
    Y_p[t - T - 1] = Y_aug[t][0:N,:]
    Y_p_i[t - T - 1] = Y_aug[t][N:N_tot,:] 



########################### Full conditional densities

@njit(parallel = True)
def updt_V(it: int64, M_Y: float64[:,:,:], M_X: float64[:,:,:],
           M_G: float64[:,:,:], W_usage: float64[:,:], init_a_V: float64,
           init_b_V: float64, init_M0: float64[:,:], init_C0: float64[:,:], 
           D_usage: float64[:,:], phi_usage: float64,
           Beta_usage: float64[:,:,:], Sigma_usage: float64[:,:]) -> float64:
  
  T_obs = M_Y.shape[0] - 1
  N_obs = M_Y.shape[1]
  q_obs = M_Y.shape[2]
  p_obs = M_X.shape[2]
  
  B_usage = nb_BR(matrix = D_usage, need_transp = True, scalar = phi_usage)
  inv_B_usage = np.linalg.inv(B_usage)
  inv_init_C0 = np.linalg.inv(init_C0)
  inv_W_usage = np.linalg.inv(W_usage)
  inv_Sigma_usage = np.linalg.inv(Sigma_usage)

  df = init_a_V + 0.5*T_obs*N_obs*q_obs + 0.5*p_obs*q_obs*(T_obs + 1)
  par, aux = np.zeros((q_obs, q_obs)), np.zeros((q_obs, q_obs))
  for t in prange(1, T_obs + 1):
    dif_t = M_Y[t] - (M_X[t] @ Beta_usage[t])
    par += (dif_t.T @ inv_B_usage) @ (dif_t @ inv_Sigma_usage)
    sub_t = Beta_usage[t] - (M_G[t] @ Beta_usage[t - 1])
    aux += (sub_t.T @ inv_W_usage) @ (sub_t @ inv_Sigma_usage)
  sub_0 = Beta_usage[0] - init_M0
  aux_0 = (sub_0.T @ inv_init_C0) @ (sub_0 @ inv_Sigma_usage)

  quantity = init_b_V + 0.5*np.trace(aux_0 + par + aux)  

  np.random.seed(it)
  inv_sample = np.random.gamma(shape = df, scale = 1/quantity)

  return 1/inv_sample

start = timer()
updt_V(it = 1, init_a_V = a_V, init_b_V = b_V, M_Y = Y, M_X = X, M_G = G,
       D_usage = D, init_M0 = M0, init_C0 = C0, W_usage = W, phi_usage = phi,
       Beta_usage = Beta, Sigma_usage = Sigma)
end = timer()
print(end - start)

@njit(parallel = True)
def updt_Sigma(it: int64, Beta_usage: float64[:,:,:], M_Y: float64[:,:,:],
               M_G: float64[:,:,:], M_X: float64[:,:,:], W_usage: float64[:,:],
               init_M0: float64[:,:], init_C0: float64[:,:], V_usage: float64,
               init_a_Sigma: float64, init_b_Sigma: float64[:,:],
               D_usage: float64[:,:], phi_usage: float64) -> float64[:,:]:
  
  T_obs = M_Y.shape[0] - 1
  N_obs = M_Y.shape[1]
  q_obs = M_Y.shape[2]
  p_obs = M_X.shape[2]
  
  B_usage = nb_BR(matrix = D_usage, need_transp = True, scalar = phi_usage)
  
  inv_V_usage = 1 / V_usage
  inv_W_usage = np.linalg.inv(W_usage)
  inv_VW_usage = inv_V_usage * inv_W_usage
  inv_B_usage = np.linalg.inv(B_usage)
  inv_VB_usage = inv_V_usage * inv_B_usage
  inv_init_C0 = np.linalg.inv(init_C0)
  dif_0 = Beta_usage[0] - init_M0
  par_init = init_b_Sigma + ((dif_0.T @ (inv_V_usage * inv_init_C0)) @ dif_0)
  par_Beta = np.zeros((q_obs, q_obs))
  par_Y = np.zeros((q_obs, q_obs))
  for t in prange(1, T_obs + 1):
    dif_Beta_t = Beta_usage[t] - (M_G[t] @ Beta_usage[t - 1])
    dif_Y_t = M_Y[t] - (M_X[t] @ Beta_usage[t])
    par_Beta += ((dif_Beta_t.T @ inv_VW_usage) @ dif_Beta_t) 
    par_Y += ((dif_Y_t.T @ inv_VB_usage) @ dif_Y_t)
  df = init_a_Sigma + p_obs + T_obs*p_obs + N_obs*T_obs
  
  sample = chol_rIW(M = par_init + par_Beta + par_Y,
                    nu = df,
                    n_seed = it)
  
  return sample

start = timer()
updt_Sigma(it = 1, init_a_Sigma = a_Sigma, init_b_Sigma = b_Sigma, D_usage = D,
           init_M0 = M0, init_C0 = C0, M_Y = Y, M_X = X, M_G = G, W_usage = W,
           Beta_usage = Beta, phi_usage = phi, V_usage = V)
end = timer()
print(end - start)

@njit(parallel = True)
def log_krnl_phi(par: float64, init_a_phi: float64, init_b_phi: float64,
                 M_Y: float64[:,:,:], M_X: float64[:,:,:], 
                 D_usage: float64[:,:], V_usage: float64,
                 Beta_usage: float64[:,:,:], 
                 Sigma_usage: float64[:,:]) -> float64:
  T_obs = M_Y.shape[0] - 1
  q_obs = M_Y.shape[2]
  if par <= 0:
    return -inf
  else:
    log_prior = (init_a_phi - 1)*log(par) - init_b_phi*par
    B_usage = nb_BR(matrix = D_usage, need_transp = True, scalar = par)
    log_det_B_usage = np.linalg.slogdet(B_usage)[1]
    aux = np.zeros((q_obs, q_obs))
    left = V_usage * B_usage
    inv_left, inv_right = np.linalg.inv(left), np.linalg.inv(Sigma_usage)
    for t in prange(1, T_obs + 1):
      dif_t = M_Y[t] - (M_X[t] @ Beta_usage[t])
      aux += (dif_t.T @ inv_left) @ (dif_t @ inv_right)
    log_LH = - 0.5*np.trace(aux) - 0.5*T_obs*q_obs*log_det_B_usage
    return log_prior + log_LH

start = timer()
log_krnl_phi(par = phi, init_a_phi = a_phi, init_b_phi = b_phi, M_Y = Y, M_X = X,
             D_usage = D, V_usage = V, Beta_usage = Beta, Sigma_usage = Sigma)
end = timer()
print(end - start)

@njit
def updt_phi(it: int64, init_a_phi: float64, init_b_phi: float64,
             D_usage: float64[:,:], M_Y: float64[:,:,:], M_X: float64[:,:,:],
             Beta_usage: float64[:,:,:], V_usage: float64, phi_usage: float64,
             Sigma_usage: float64[:,:], delta2: float64) -> float64:
  
  # Proposing a new value from the inverse-normal distribution
  phi_prop = rIN(n_seed = it, mu = phi_usage, sigma2 = delta2)
  
  # Generating a random number between zero and one in log scale
  log_u = log_rUnif(it)
  
  # Computing the log of the Metropolis ratio
  log_num = +logpdf_IN(mean = phi_prop,
                       sigma2 = delta2,
                       x = phi_usage) + log_krnl_phi(par = phi_prop,
                                                     init_a_phi = init_a_phi,
                                                     init_b_phi = init_b_phi,
                                                     D_usage = D_usage,
                                                     M_Y = M_Y,
                                                     M_X = M_X,
                                                     Beta_usage = Beta_usage,
                                                     V_usage = V_usage,
                                                     Sigma_usage = Sigma_usage)
  log_den = +logpdf_IN(mean = phi_usage,
                       sigma2 = delta2,
                       x = phi_prop) + log_krnl_phi(par = phi_usage,
                                                    init_a_phi = init_a_phi,
                                                    init_b_phi = init_b_phi,
                                                    D_usage = D_usage,
                                                    M_Y = M_Y,
                                                    M_X = M_X,
                                                    Beta_usage = Beta_usage,
                                                    V_usage = V_usage,
                                                    Sigma_usage = Sigma_usage)
  log_rho = log_num - log_den
  log_eta = min([0, log_rho])
  
  if isfinite(log_eta) and log_u <= log_eta:
    return phi_prop
  else:
    return phi_usage

start = timer()
updt_phi(it = 3, init_a_phi = a_phi, init_b_phi = b_phi, D_usage = D, M_Y = Y,
         M_X = X, Beta_usage = Beta, V_usage = V, phi_usage = phi,
         Sigma_usage = Sigma, delta2 = 0.001)
end = timer()
print(end - start)

@njit
def FF(init_M0: float64[:,:], init_C0: float64[:,:], M_Y: float64[:,:,:],
       M_X: float64[:,:,:], M_G: float64[:,:,:], V_usage: float64,
       phi_usage: float64, D_usage: float64[:,:],
       W_usage: float64[:,:]) -> (float64[:,:,:], float64[:,:,:]):

  T_obs = M_Y.shape[0] - 1
  q_obs, p_obs = M_Y.shape[2], M_X.shape[2]  
  
  B_usage = nb_BR(matrix = D_usage, need_transp = True, scalar = phi_usage)  
  cov_usage = V_usage * B_usage
  M_M = np.zeros((T_obs + 1, p_obs, q_obs))
  M_C = np.zeros((T_obs + 1, p_obs, p_obs))
  for t in range(0, T_obs + 1):
    if t == 0:
      M_M[t] = init_M0
      M_C[t] = V_usage * init_C0
    else:
      E_t = V_usage * W_usage + (M_G[t] @ M_C[t - 1]) @ M_G[t].T
      Q_t = cov_usage + (M_X[t] @ E_t) @ M_X[t].T
      inv_Q_t = np.linalg.inv(Q_t)
      aux_prod_t = (E_t @ M_X[t].T) @ inv_Q_t
      aux_dif_t = M_Y[t] - (M_X[t] @ M_G[t]) @ M_M[t - 1]
      M_M[t] = (M_G[t] @ M_M[t - 1]) + (aux_prod_t @ aux_dif_t)
      M_C[t] = E_t - (aux_prod_t @ M_X[t]) @ E_t
      
  return M_M, M_C

@njit
def BS_t(C_now: float64[:,:], M_now: float64[:,:], G_now: float64[:,:],
         G_next: float64[:,:], Beta_next: float64[:,:], 
         V: float64, W: float64[:,:]) -> (float64[:,:], float64[:,:]):
  inv_VW = (1 / V) * np.linalg.inv(W)
  inv_C_now = np.linalg.inv(C_now)
  inv_H_now = inv_C_now + ((G_now.T @ inv_VW) @ G_now)
  H_now = np.linalg.inv(inv_H_now)
  h_now = (inv_C_now @ M_now) + ((G_next @ inv_VW) @ Beta_next)
  avg_now = H_now @ h_now
  return avg_now, H_now

@njit
def updt_Beta(it: int64, init_M0: float64[:,:], init_C0: float64[:,:],
              M_Y: float64[:,:,:], M_X: float64[:,:,:], M_G: float64[:,:,:],
              D_usage: float64[:,:], W_usage: float64[:,:], V_usage: float64,
              phi_usage: float64, Sigma_usage: float64[:,:]) -> float64[:,:,:]:

  T_obs = M_Y.shape[0] - 1
  q_obs, p_obs = M_Y.shape[2], M_X.shape[2]
  
  M_M, M_C = FF(init_M0 = init_M0,
                init_C0 = init_C0,
                M_Y = M_Y,
                M_X = M_X,
                M_G = M_G,
                V_usage = V_usage,
                W_usage = W_usage,
                D_usage = D_usage,
                phi_usage = phi_usage)
  
  sample = np.zeros((T_obs + 1, p_obs, q_obs))
  for x in range(0, T_obs + 1):
    t = T_obs - x
    if t == T_obs:
      sample[t] = chol_rMVN(avg = M_M[t],
                            left = M_C[t],
                            right = Sigma_usage,
                            n_seed = it)
    else:
      aux_BS_t = BS_t(C_now = M_C[t],
                      M_now = M_M[t],
                      G_now = M_G[t],
                      G_next = M_G[t + 1],
                      Beta_next = sample[t + 1],
                      W = W_usage,
                      V = V_usage)
      sample[t] = chol_rMVN(avg = aux_BS_t[0],
                            left = aux_BS_t[1],
                            right = Sigma_usage,
                            n_seed = it)
  
  return sample

start = timer()
updt_Beta(it = 5, init_M0 = M0, init_C0 = C0, M_Y = Y, M_X = X, M_G = G,
          D_usage = D, W_usage = W, phi_usage = phi, V_usage = V,
          Sigma_usage = Sigma)[20]
end = timer()
print(end - start)



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



########################### Acceptance ratio, given in %
round((acc_phi / num_iter)*100, 2)
np.round((acc_D / num_iter)*100, 2)



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



#### Plotting some figures

import matplotlib.pyplot as plt

## phi
# Vectorization
vec = e_phi[0]
for j in range(1, int(smp_size)):
  now = e_phi[j]
  vec = np.concatenate((vec, now), axis = None)

# Multiple line plot
plt.plot(np.arange(0, int(smp_size)),
         vec,
         linewidth = 1)
plt.xlabel('Iteration')
plt.ylabel('Sampled value',
           multialignment = 'center')
plt.title(r'$\phi$')

# Histogram
plt.hist(x = vec,
         density = True,
         rwidth = 0.9)
plt.ylabel('Probability')
plt.title(r'$\phi$')
plt.axvline(x = phi,
            color = 'r',
            linestyle = '-')

# Histogram
plt.hist(x = vec,
         density = True,
         rwidth = 0.9)
plt.ylabel('Probability')
plt.title(r'$\phi$')
plt.axvline(x = phi*sqrt((Lambda**2).sum()),
            color = 'r',
            linestyle = '-')

## V*Sigma
# Vectorization
r, s = 1, 1
vec = e_V[0]*e_Sigma[0][r,s]
for j in range(1, int(smp_size)):
  now = e_V[j]*e_Sigma[j][r,s]
  vec = np.concatenate((vec, now), axis = None)

# Multiple line plot
plt.plot(np.arange(0, int(smp_size)),
         vec,
         linewidth = 1)
plt.xlabel('Iteration')
plt.ylabel('Sampled value',
           multialignment = 'center')
plt.title(r'$V \cdot \Sigma$')

# Histogram
plt.hist(x = vec,
         density = True,
         rwidth = 0.9)
plt.ylabel('Probability')
plt.title(r'$V \cdot \Sigma$')
plt.axvline(x = V*Sigma[r,s],
            color = 'r',
            linestyle = '-')

## V*Sigma*phi
# Vectorization
r, s = 1, 1
vec = e_V[0]*e_Sigma[0][r,s]*e_phi[0]
for j in range(1, int(smp_size)):
  now = e_V[j]*e_Sigma[j][r,s]*e_phi[j]
  vec = np.concatenate((vec, now), axis = None)

# Multiple line plot
plt.plot(np.arange(0, int(smp_size)),
         vec,
         linewidth = 1)
plt.xlabel('Iteration')
plt.ylabel('Sampled value',
           multialignment = 'center')
plt.title(r'$V \cdot \Sigma$')

# Histogram
plt.hist(x = vec,
         density = True,
         rwidth = 0.9)
plt.ylabel('Probability')
plt.title(r'$V \cdot \Sigma \cdot \phi$')
plt.axvline(x = V*Sigma[r,s]*(phi)/sqrt((Lambda**2).sum()),
            color = 'r',
            linestyle = '-')

## Beta_t[k,i]
t = 50
k = 0
i = 2 - 1
# Vectorization
vec = e_Beta[0][t,k,i]
for j in range(1, int(smp_size)):
  now = e_Beta[j][t,k,i]
  vec = np.concatenate((vec, now),
                       axis = None)

# Multiple line plot
plt.plot(np.arange(0, int(smp_size)),
         vec,
         linewidth = 1)
plt.xlabel('Iteration')
plt.ylabel('Sampled value',
           multialignment = 'center')
plt.title(r'$\beta_t$')

# Histogram
plt.hist(x = vec,
         density = True,
         rwidth = 0.9)
plt.ylabel('Probability')
plt.title(r'$\beta_t$')
plt.axvline(x = Beta[t,k,i],
            color = 'r',
            linestyle = '-')
