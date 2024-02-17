# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 17:26:50 2021

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
S = np.vstack([[0.0, 0.0], [1.0, 1.0],
               [0.0, 1/3], [0.0, 2/3], [0.0, 1.0],
               [1/3, 0.0], [1/3, 1/3], [1/3, 2/3], [1/3, 1.0],
               [1/2, 1/2],
               [2/3, 0.0], [2/3, 1/3], [2/3, 2/3], [2/3, 1.0],
               [1.0, 0.0], [1.0, 1/3], [1.0, 2/3]]).T

# The number of sites and the number of replications (or times)
N = S.shape[1]
T = 100 # Run with 10, 100 and 1000
p, q = 2, 2
seed_v = 1000

D = np.zeros((2, N))
D[:,0:2] = S[:,0:2]

psi = -2*log(0.05)/np.max(distance.cdist(S.T, S.T, metric = 'sqeuclidean'))
R_d = nb_BR(matrix = S,
            need_transp = True,
            squared = True,
            scalar = psi)
R_12 = R_d[0:2,0:2]
R_3N = R_d[2:18,2:18]
R_star = R_d[2:18,0:2]
S_cond = S[:,2:(N + 1)] + np.matmul(D[:,0:2] - S[:,0:2],
                                    np.linalg.inv(R_12) @ R_star.T)
R_cond = R_3N - np.matmul(R_star, np.linalg.inv(R_12) @ R_star.T)

sigma2d = np.array([[0.500, 0.000],
                    [0.000, 0.500]])
 
D[:,2:(N + 1)] = rMVN(n = 1,
                      avg = S_cond,
                      left = sigma2d,
                      right = R_cond,
                      n_seed = seed_v)

V = 0.4
a_V, b_V = 0.001, 0.001

Sigma = np.array([[1.0, 0.7],
                  [0.7, 1.0]])
a_Sigma, b_Sigma = 0.001, 0.001*np.identity(q)

phi = 0.1
a_phi, b_phi = 0.001, 0.001

# Using the true deformations to specify B
B = nb_BR(matrix = D, need_transp = True, scalar = phi)

M0, C0 = np.zeros((p, q)), 1.0*np.identity(p)
W = 1.0*np.identity(p)
Beta = np.zeros((T + 1, p, q))

G = np.zeros((T + 1, p, p))
for t in range(0, T + 1):
  if t == 0:
    Beta[t] = rMVN(n = 1,
                   avg = M0,
                   left = V*C0,
                   right = Sigma,
                   n_seed = seed_v + t)
  else:
    G[t] = np.identity(p)
    Beta[t] = rMVN(avg = np.matmul(G[t], Beta[t-1]),
                   left = V*W,
                   right = Sigma,
                   n = 1,
                   n_seed = seed_v + t)

# Likelihood function
Y, X = np.zeros((T + 1, N, q)), np.zeros((T + 1, N, p))
for t in range(1, T + 1):
  X[t] = np.vstack([np.ones(N),
                    uniform.rvs(loc = 0,
                                scale = 1,
                                size = N,
                                random_state = seed_v + t)]).T
  Y[t] = rMVN(n = 1,
              avg = np.matmul(X[t], Beta[t]),
              left = V*B,
              right = Sigma,
              n_seed = seed_v + t)

# R correlation matrix
def comp_R(S_obs, psi_obs):
  if isfinite(psi_obs):
    R_N = nb_BR(matrix = S_obs,
                need_transp = True,
                squared = True,
                scalar = psi_obs)
  else:
    N_obs = S_obs.shape[1]
    R_N = np.identity(N_obs)
  return R_N

# Configuration of the distortion level
psi_usage = 5.0
R_usg = comp_R(S_obs = S, psi_obs = psi_usage)
sigma2d_usg = np.cov(S)



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
    log_det_B_usage = log(np.linalg.det(B_usage))
    aux = np.zeros((q_obs, q_obs))
    left = V_usage * B_usage
    inv_left, inv_right = np.linalg.inv(left), np.linalg.inv(Sigma_usage)
    for t in prange(1, T_obs + 1):
      dif_t = M_Y[t] - (M_X[t] @ Beta_usage[t])
      aux += (dif_t.T @ inv_left) @ (dif_t @ inv_right)
    log_LH = - 0.5*np.trace(aux) - 0.5*T_obs*q_obs*log_det_B_usage
    return log_prior + log_LH

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
  log_alpha = min([0, log_rho])
  
  if isfinite(log_alpha) and log_u <= log_alpha:
    return phi_prop
  else:
    return phi_usage

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

@njit(parallel = True)
def log_krnl_D(par: float64[:,:], S_obs: float64[:,:], M_X: float64[:,:,:],
               M_Y: float64[:,:,:], R_usage: float64[:,:],
               sigma2d_usage: float64[:,:], Sigma_usage: float64[:,:],
               V_usage: float64, Beta_usage: float64[:,:,:],
               phi_usage: float64) -> float64:
  T_obs = M_Y.shape[0] - 1
  q_obs = M_Y.shape[2]    
  dif = par - S_obs
  inv_R_usage = np.linalg.inv(R_usage)
  inv_sigma2d_usage = np.linalg.inv(sigma2d_usage)
  krnl_prior = (dif.T @ inv_sigma2d_usage) @ (dif @ inv_R_usage)
  log_prior = -0.5*np.trace(krnl_prior)

  B_usage = nb_BR(matrix = par, need_transp = True, scalar = phi_usage)
  log_det_B_usage = log(np.linalg.det(B_usage))

  aux = np.zeros((q_obs, q_obs))
  left = V_usage * B_usage
  inv_left, inv_right = np.linalg.inv(left), np.linalg.inv(Sigma_usage)
  for t in prange(1, T_obs + 1):
    dif_t = M_Y[t] - (M_X[t] @ Beta_usage[t])
    aux += (dif_t.T @ inv_left) @ (dif_t @ inv_right)
  log_LH = - 0.5*np.trace(aux) - 0.5*T_obs*q_obs*log_det_B_usage
  
  return log_prior + log_LH
