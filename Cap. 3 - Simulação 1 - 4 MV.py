# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 15:27:30 2023

@author: rodri
"""

########################### Python libraries

import numpy as np
from scipy.stats import uniform, multivariate_normal
from timeit import default_timer as timer
from math import pi, exp, log, sqrt, isfinite, inf, dist
from numba import njit, float64, boolean, bool_, uint8, int64, prange
from numba.core.errors import NumbaPerformanceWarning
from warnings import simplefilter
from scipy.spatial import distance, distance_matrix
from scipy.linalg import cholesky
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

# Obtaining Y_obs_t, Y_mis_t, N_obs_t, N_mis_t, P_t, L_obs_t and L_mis_t
@njit
def Permut(vector: float64[:]):
  N = len(vector)
  Positions = [0]*N
  Identify = [0]*N
  
  i = 0
  for n in range(0, N):
    if np.isnan(vector[n]):
      pass
    else:
      i += 1
      Positions[n] = i
      Identify[n] = 1
      
  N_o = i
  N_m = N - N_o
  
  j = i
  for n in range(0, N):
    if np.isnan(vector[n]):
      j += 1
      Positions[n] = j
    else:
      pass
  
  Seq = np.array(Positions, dtype = np.int64) - 1
  P_t = np.eye(N)[:,Seq]
  L_o_t = np.hstack((np.eye(N_o), np.zeros((N_o, N_m))))
  L_m_t = np.hstack((np.zeros((N_m, N_o)), np.eye(N_m)))
  Subseq_o = np.array(Identify, dtype = bool_)
  Subseq_m = np.logical_not(Subseq_o)
  Y_o_t = vector[Subseq_o].reshape(N_o, 1)
  Y_m_t = vector[Subseq_m].reshape(N_m, 1)
  
  return Y_o_t, Y_m_t, N_o, N_m, P_t, L_o_t, L_m_t, Subseq_o

# Kronecker product
@njit
def nb_kp(X: float64[:,:], Y: float64[:,:]) -> float64[:,:]:
  return np.kron(X, Y)

# Sampling from multivariate normal distribution (fast, but less robust)
@njit("float64[:,:](int64, float64[:,:], float64[:,:])")
def chol_rmvN(n_seed: int64, avg: float64[:,:],
              cov: float64[:,:]) -> float64[:,:]:
  p = cov.shape[0]
  np.random.seed(n_seed)
  gen = np.linalg.cholesky(cov) @ np.random.standard_normal(p).reshape(p, 1)
  sample = avg + gen
  return sample

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

# Sites in the unit square
S = np.vstack([[1/5, 1/5], [4/5, 4/5], [1/5, 4/5], [4/5, 1/5],
                           [1/5, 2/5], [1/5, 3/5], 
               [2/5, 1/5], [2/5, 2/5], [2/5, 3/5], [2/5, 4/5],
               [3/5, 1/5], [3/5, 2/5], [3/5, 3/5], [3/5, 4/5],
                           [4/5, 2/5], [4/5, 3/5]]).T

# The number of sites and the number of replications (or times)
N = S.shape[1]
T = 200
p, q = 2, 2 
seed_v = 500

V = 0.1
a_V, b_V = 0.001, 0.001

Sigma = np.array([[1.00, 0.75],
                  [0.75, 1.00]])
a_Sigma, b_Sigma = 0.001, 0.001*np.identity(q)

phi = 0.3
a_phi, b_phi = 0.001, 0.001

### B and B_i are specified by an anisotropic structure
D = np.zeros((2, N))
D[:,0:2] = S[:,0:2]

psi = 1/(2*np.max(distance.cdist(S.T, S.T, metric = 'sqeuclidean')))
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

sigma2d = np.array([[0.200, 0.000],
                    [0.000, 0.200]])
 
D[:,2:(N + 1)] = rMVN(n = 1,
                      avg = S_cond,
                      left = sigma2d,
                      right = R_cond,
                      n_seed = seed_v)

B = nb_BR(matrix = D, need_transp = True, scalar = phi)

M0, C0 = np.zeros((p, q)), 1.0*np.identity(p)
m0 = (M0.T).reshape((p*q, 1))
W = 1.0*np.identity(p)
Beta, beta = np.zeros((T + 1, p, q)), np.zeros((T + 1, p*q, 1))
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
    beta[t] = (Beta[t].T).reshape((p*q, 1))

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

# x e g
x, g = np.zeros((T+1,N*q,p*q)), np.zeros((T+1,p*q,p*q))
for t in range(0, T + 1):
  g[t] = nb_kp(X = np.eye(q), Y = G[t])
  x[t] = nb_kp(X = np.eye(q), Y = X[t])

# Inserting two some missing values in each column and vectorizing
copy_Y = np.copy(Y)
for t in range(1, T + 1):
  for i in range(0, q):
    np.random.seed(t + i)
    aux = np.random.choice(N, 4, replace = False)
    copy_Y[t][aux, i] = np.array([np.nan, np.nan, np.nan, np.nan])
                                            
y = np.zeros((T + 1, N*q, 1))
for t in range(1, T + 1):
  y[t] = (copy_Y[t].T).reshape((N*q, 1))

# Computing N_obs_t, P_t and Id_obs_t
N_o = np.zeros(T + 1, dtype = np.int32)
N_m = np.zeros(T + 1, dtype = np.int32)
P = np.zeros((T + 1, N*q, N*q))
inv_P = np.zeros((T + 1, N*q, N*q))
Id_o = np.zeros((T + 1, N*q), dtype = np.bool_)
for t in range(1, T + 1):
  aux_t = Permut(y[t])
  N_o[t] = aux_t[2]
  N_m[t] = aux_t[3]
  P[t] = aux_t[4]
  inv_P[t] = np.linalg.inv(P[t])
  Id_o[t] = aux_t[7]

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
psi_usage = 9.85 # 9.9
R_usg = comp_R(S_obs = S, psi_obs = psi_usage)
sigma2d_usg = np.cov(S)



########################### Full conditional densities

@njit(parallel = True)
def updt_y(it: int64, m_y: float64[:,:,:], M_x: float64[:,:,:],
           D_usage: float64[:,:], phi_usage: float64, V_usage: float64, 
           beta_usage: float64[:,:,:], Sigma_usage: float64[:,:],
           num_obs: int64[:], num_mis: int64[:], M_P: float64[:,:,:],
           M_inv_P: float64[:,:,:], Id_obs: bool_[:,:]) -> float64[:,:,:]:
  T_obs = beta_usage.shape[0] - 1
  q_obs = Sigma_usage.shape[0] 
  p_obs = int(beta_usage.shape[1] / q_obs)
  N_obs = D_usage.shape[1]
  B_usg = nb_BR(matrix = D_usage, need_transp = True, scalar = phi_usage)
  Sig_B_usg = np.kron(Sigma_usage, B_usg)
  vec_y = np.zeros((T_obs + 1, N_obs*q_obs, 1))
  M_y = np.zeros((T_obs + 1, N_obs, q_obs))
  for t in prange(1, T_obs + 1):
    if num_mis[t] == 0:
      vec_y[t] = m_y[t]
    elif num_mis[t] < N_obs*q_obs:
      vec_beta_t = beta_usage[t].reshape(p_obs*q_obs, 1)
      mu_t = M_P[t] @ (M_x[t] @ vec_beta_t)
      mu_o_t = mu_t[0:num_obs[t]]
      mu_m_t = mu_t[num_obs[t]:(N_obs*q_obs)]
      Delta_t = (M_P[t] @ Sig_B_usg) @ M_P[t].T
      Delta_oo_t = Delta_t[0:num_obs[t],0:num_obs[t]]
      Delta_om_t = Delta_t[0:num_obs[t],num_obs[t]:(N_obs*q_obs)]
      Delta_mo_t = Delta_om_t.T
      Delta_mm_t = Delta_t[num_obs[t]:(N_obs*q_obs),
                           num_obs[t]:(N_obs*q_obs)]
      inv_Delta_oo_t = np.linalg.inv(Delta_oo_t)
      y_obs_t = m_y[t][Id_obs[t]]
      dif_o_t = y_obs_t - mu_o_t
      aux_mult = Delta_mo_t @ inv_Delta_oo_t
      mu_cond_t = mu_m_t + (aux_mult @ dif_o_t)
      Delta_cond_t = Delta_mm_t - (aux_mult @ Delta_om_t)
      vec_m_j_t = chol_rmvN(n_seed = it + t,
                            avg = mu_cond_t,
                            cov = V_usage * Delta_cond_t)
      vec_y[t] = M_inv_P[t] @ np.vstack((y_obs_t, vec_m_j_t)) 
    else:
      vec_beta_t = beta_usage[t].reshape(p_obs*q_obs, 1)
      vec_y[t] = chol_rmvN(n_seed = it + t,
                           avg = M_x[t] @ vec_beta_t,
                           cov = V_usage * Sig_B_usg)
  return vec_y

@njit(parallel = True)
def updt_V(it: int64, M_x: float64[:,:,:], M_g: float64[:,:,:],
           W_usage: float64[:,:], init_a_V: float64, init_b_V: float64,
           init_m0: float64[:,:], init_C0: float64[:,:], D_usage: float64[:,:],
           phi_usage: float64, beta_usage: float64[:,:,:],
           Sigma_usage: float64[:,:], M_y: float64[:,:,:]) -> float64:
  T_obs = beta_usage.shape[0] - 1
  q_obs = Sigma_usage.shape[0]
  p_obs = W_usage.shape[0]
  N_obs = D_usage.shape[1]
  
  B_usage = nb_BR(matrix = D_usage, need_transp = True, scalar = phi_usage)
  Sig_B = np.kron(Sigma_usage, B_usage)
  Sig_C0 = np.kron(Sigma_usage, init_C0)
  Sig_W = np.kron(Sigma_usage, W_usage)
  inv_Sig_B = np.linalg.inv(Sig_B)
  inv_Sig_C0 = np.linalg.inv(Sig_C0)
  inv_Sig_W = np.linalg.inv(Sig_W)

  df = init_a_V + 0.5*(p_obs*q_obs + T_obs*p_obs*q_obs + T_obs*N_obs*q_obs)
  aux, par = 0.0, 0.0
  for t in prange(1, T_obs + 1):
    sub_t = beta_usage[t] - (M_g[t] @ beta_usage[t - 1])
    aux += ((sub_t.T @ inv_Sig_W) @ sub_t).item()
    mu_obs_t = M_x[t] @ beta_usage[t]
    y_obs_t = M_y[t]
    dif_obs_t = y_obs_t - mu_obs_t
    par += ((dif_obs_t.T @ inv_Sig_B) @ dif_obs_t).item()
  sub_0 = beta_usage[0] - init_m0
  aux_0 = ((sub_0.T @ inv_Sig_C0) @ sub_0).item()

  quantity = init_b_V + 0.5*(aux_0 + aux + par)

  np.random.seed(it)
  inv_sample = np.random.gamma(shape = df, scale = 1/quantity)

  return 1/inv_sample

@njit(parallel = True)
def updt_Sigma(it: int64, beta_usage: float64[:,:,:], M_y: float64[:,:,:],
               M_G: float64[:,:,:], M_X: float64[:,:,:], W_usage: float64[:,:],
               init_M0: float64[:,:], init_C0: float64[:,:], V_usage: float64,
               init_a_Sigma: float64, init_b_Sigma: float64[:,:],
               D_usage: float64[:,:], phi_usage: float64) -> float64[:,:]:
  
  T_obs = M_y.shape[0] - 1
  N_obs = D_usage.shape[1]
  q_obs = init_b_Sigma.shape[1]
  p_obs = M_X.shape[2]
  
  B_usage = nb_BR(matrix = D_usage, need_transp = True, scalar = phi_usage)
  
  M_Y = np.zeros((T_obs + 1, N_obs, q_obs))
  Beta_usage = np.zeros((T_obs + 1, p_obs, q_obs))
  M_Y[0] = (M_y[0].reshape(q_obs, N_obs)).T 
  Beta_usage[0] = (beta_usage[0].reshape(q_obs, p_obs)).T
  
  inv_V_usage = 1 / V_usage
  inv_W_usage = np.linalg.inv(W_usage)
  inv_VW_usage = inv_V_usage * inv_W_usage
  inv_B_usage = np.linalg.inv(B_usage)
  inv_VB_usage = inv_V_usage * inv_B_usage
  inv_init_C0 = np.linalg.inv(init_C0)
  dif_0 = Beta_usage[0] - init_M0
  par_init = init_b_Sigma + ((dif_0.T @ (inv_V_usage * inv_init_C0)) @ dif_0)
  par_beta = np.zeros((q_obs, q_obs))
  par_Y = np.zeros((q_obs, q_obs))
  for t in prange(1, T_obs + 1):
    M_Y[t] = (M_y[t].reshape(q_obs, N_obs)).T 
    Beta_usage[t] = (beta_usage[t].reshape(q_obs, p_obs)).T
    dif_beta_t = Beta_usage[t] - (M_G[t] @ Beta_usage[t - 1])
    dif_Y_t = M_Y[t] - (M_X[t] @ Beta_usage[t])
    par_beta += ((dif_beta_t.T @ inv_VW_usage) @ dif_beta_t) 
    par_Y += ((dif_Y_t.T @ inv_VB_usage) @ dif_Y_t)
  df = init_a_Sigma + p_obs + T_obs*p_obs + N_obs*T_obs
  
  sample = chol_rIW(M = par_init + par_beta + par_Y,
                    nu = df,
                    n_seed = it)
  
  return sample

@njit
def FF(init_m0: float64[:,:], init_C0: float64[:,:], M_y: float64[:,:,:],
       M_x: float64[:,:,:], M_g: float64[:,:,:], V_usage: float64,
       phi_usage: float64, D_usage: float64[:,:], Sigma_usage: float64[:,:],
       W_usage: float64[:,:]) -> (float64[:,:,:], float64[:,:,:]):

  T_obs = M_y.shape[0] - 1
  q_obs = Sigma_usage.shape[0]
  p_obs = W_usage.shape[0]
  N_obs = D_usage.shape[1]
  
  B_usage = nb_BR(matrix = D_usage, need_transp = True, scalar = phi_usage)
  VSig_W = V_usage * np.kron(Sigma_usage, W_usage)
  VSig_B = V_usage * np.kron(Sigma_usage, B_usage)
  
  cov_usage = V_usage * B_usage
  M_m = np.zeros((T_obs + 1, p_obs*q_obs, 1))
  M_C = np.zeros((T_obs + 1, p_obs*q_obs, p_obs*q_obs))
  M_m[0] = init_m0
  M_C[0] = V_usage * np.kron(Sigma_usage, init_C0)
  for t in range(1, T_obs + 1):
    E_t = VSig_W + (M_g[t] @ M_C[t-1]) @ M_g[t].T
    Q_t = VSig_B + (M_x[t] @ E_t) @ M_x[t].T
    inv_Q_t = np.linalg.inv(Q_t)
    a_t = M_g[t] @ M_m[t-1]
    aux_prod_t = (E_t @ M_x[t].T) @ inv_Q_t
    aux_dif_t = M_y[t] - (M_x[t] @ a_t)
    M_m[t] = a_t + (aux_prod_t @ aux_dif_t)
    M_C[t] = E_t - (aux_prod_t @ M_x[t]) @ E_t 
  return M_m, M_C

@njit
def BS_t(C_now: float64[:,:], m_now: float64[:,:], g_now: float64[:,:],
         g_next: float64[:,:], beta_next: float64[:,:], Sig_now: float64[:,:],
         V: float64, W: float64[:,:]) -> (float64[:,:], float64[:,:]):
  inv_VSig_W = (1 / V) * np.linalg.inv(np.kron(Sig_now, W))
  inv_C_now = np.linalg.inv(C_now)
  inv_H_now = inv_C_now + ((g_now.T @ inv_VSig_W) @ g_now)
  H_now = np.linalg.inv(inv_H_now)
  h_now = (inv_C_now @ m_now) + ((g_next @ inv_VSig_W) @ beta_next)
  avg_now = H_now @ h_now
  return avg_now, H_now

@njit
def updt_beta(it: int64, init_m0: float64[:,:], init_C0: float64[:,:],
              M_y: float64[:,:,:], M_x: float64[:,:,:], M_g: float64[:,:,:],
              D_usage: float64[:,:], W_usage: float64[:,:], V_usage: float64,
              phi_usage: float64, Sigma_usage: float64[:,:]) -> float64[:,:,:]:

  T_obs = M_y.shape[0] - 1
  q_obs = Sigma_usage.shape[0]
  p_obs = W_usage.shape[0]
  
  M_m, M_C = FF(init_m0 = init_m0,
                init_C0 = init_C0,
                M_y = M_y,
                M_x = M_x,
                M_g = M_g,
                V_usage = V_usage,
                W_usage = W_usage,
                D_usage = D_usage,
                phi_usage = phi_usage,
                Sigma_usage = Sigma_usage)
  
  sample = np.zeros((T_obs + 1, p_obs*q_obs, 1))
  for x in range(0, T_obs + 1):
    t = T_obs - x
    if t == T_obs:
      sample[t] = chol_rmvN(avg = M_m[t],
                            cov = M_C[t],
                            n_seed = it)
    else:
      aux_BS_t = BS_t(C_now = M_C[t],
                      m_now = M_m[t],
                      g_now = M_g[t],
                      g_next = M_g[t + 1],
                      Sig_now = Sigma_usage,
                      beta_next = sample[t + 1],
                      W = W_usage,
                      V = V_usage)
      sample[t] = chol_rmvN(avg = aux_BS_t[0],
                            cov = aux_BS_t[1],
                            n_seed = it)
  
  return sample

@njit(parallel = True)
def log_krnl_phi(par: float64, init_a_phi: float64, init_b_phi: float64,
                 M_y: float64[:,:,:], M_x: float64[:,:,:], V_usage: float64,
                 D_usage: float64[:,:], beta_usage: float64[:,:,:],
                 Sigma_usage: float64[:,:]) -> float64:
  T_obs = beta_usage.shape[0] - 1
  q_obs = Sigma_usage.shape[0]
  N_obs = D_usage.shape[1]
  log_lh = 0.0
  if par > 0:
    log_prior = (init_a_phi - 1)*log(par) - init_b_phi*par
    B_usage = nb_BR(matrix = D_usage, need_transp = True, scalar = par)
    ldet_B_usage = log(np.linalg.det(B_usage))
    Sig_B = np.kron(Sigma_usage, B_usage)
    inv_Sig_B = np.linalg.inv(Sig_B)
    for t in prange(1, T_obs + 1):
      mu_obs_t = M_x[t] @ beta_usage[t]
      y_obs_t = M_y[t]
      dif_obs_t = y_obs_t - mu_obs_t
      log_lh += ((dif_obs_t.T @ inv_Sig_B) @ dif_obs_t).item()
    return log_prior - 0.5*(1 / V_usage)*log_lh - 0.5*T_obs*q_obs*ldet_B_usage
  else:
    return -inf

@njit
def updt_phi(it: int64, init_a_phi: float64, init_b_phi: float64,
             D_usage: float64[:,:], M_y: float64[:,:,:], M_x: float64[:,:,:],
             beta_usage: float64[:,:,:], V_usage: float64, phi_usage: float64,
             Sigma_usage: float64[:,:], delta2: float64) -> float64:
  
  # Proposing a new value from the inverse-normal distribution
#  np.random.seed(it)
#  phi_prop = exp(np.random.normal(loc = log(phi_usage), scale = delta2))
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
                                                     M_y = M_y,
                                                     M_x = M_x,
                                                     beta_usage = beta_usage,
                                                     V_usage = V_usage,
                                                     Sigma_usage = Sigma_usage)
  log_den = +logpdf_IN(mean = phi_usage,
                       sigma2 = delta2,
                       x = phi_prop) + log_krnl_phi(par = phi_usage,
                                                    init_a_phi = init_a_phi,
                                                    init_b_phi = init_b_phi,
                                                    D_usage = D_usage,
                                                    M_y = M_y,
                                                    M_x = M_x,
                                                    beta_usage = beta_usage,
                                                    V_usage = V_usage,
                                                    Sigma_usage = Sigma_usage)
  log_rho = log_num - log_den
  log_alpha = min([0, log_rho])
  
  if isfinite(log_alpha) and log_u <= log_alpha:
    return phi_prop
  else:
    return phi_usage

@njit(parallel = True)
def log_krnl_D(par: float64[:,:], S_obs: float64[:,:], M_x: float64[:,:,:],
               M_y: float64[:,:,:], R_usage: float64[:,:],
               sigma2d_usage: float64[:,:], Sigma_usage: float64[:,:],
               V_usage: float64, beta_usage: float64[:,:,:],
               phi_usage: float64) -> float64:
  T_obs = beta_usage.shape[0] - 1
  q_obs = Sigma_usage.shape[0]
  N_obs = par.shape[1]
  log_lh = 0.0    
  dif = par - S_obs
  inv_R_usage = np.linalg.inv(R_usage)
  inv_sigma2d_usage = np.linalg.inv(sigma2d_usage)
  krnl_prior = (dif.T @ inv_sigma2d_usage) @ (dif @ inv_R_usage)
  log_prior = np.trace(krnl_prior)

  B_usage = nb_BR(matrix = par, need_transp = True, scalar = phi_usage)
  Sig_B = np.kron(Sigma_usage, B_usage)
  ldet_B = log(np.linalg.det(B_usage))
  inv_Sig_B = np.linalg.inv(Sig_B)
  for t in prange(1, T_obs + 1):
    mu_obs_t = M_x[t] @ beta_usage[t]
    y_obs_t = M_y[t]
    dif_obs_t = y_obs_t - mu_obs_t
    log_lh += ((dif_obs_t.T @ inv_Sig_B) @ dif_obs_t).item()
  
  return -0.5*(log_prior + T_obs*q_obs*ldet_B + (1 / V_usage)*log_lh)

@njit
def updt_D(it: int64, S_obs: float64[:,:], M_x: float64[:,:,:],
           M_y: float64[:,:,:], beta_usage: float64[:,:,:], V_usage: float64,
           phi_usage: float64, sigma2d_usage: float64[:,:],
           D_usage: float64[:,:], R_usage: float64[:,:],
           Sigma_usage: float64[:,:]) -> float64[:,:]:
  N_obs = D_usage.shape[1]
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
                       M_x = M_x,
                       M_y = M_y,
                       Sigma_usage = Sigma_usage,
                       beta_usage = beta_usage,
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
                                              M_x = M_x,
                                              M_y = M_y,
                                              Sigma_usage = Sigma_usage,
                                              beta_usage = beta_usage,
                                              V_usage = V_usage,
                                              phi_usage = phi_usage,
                                              sigma2d_usage = sigma2d_usage,
                                              R_usage = R_usage):
          L = L - w
          L_D[i,j] = L
        while R < upper and logy < log_krnl_D(par = R_D,
                                              S_obs = S_obs,
                                              M_x = M_x,
                                              M_y = M_y,
                                              Sigma_usage = Sigma_usage,
                                              beta_usage = beta_usage,
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
                       M_x = M_x,
                       M_y = M_y,
                       Sigma_usage = Sigma_usage,
                       beta_usage = beta_usage,
                       V_usage = V_usage,
                       phi_usage = phi_usage,
                       sigma2d_usage = sigma2d_usage,
                       R_usage = R_usage)
  
      while gx1 < logy:
        x1 = np.random.uniform(L, R)
        x1_D[i,j] = x1
        gx1 = log_krnl_D(par = x1_D,
                         S_obs = S_obs,
                         M_x = M_x,
                         M_y = M_y,
                         Sigma_usage = Sigma_usage,
                         beta_usage = beta_usage,
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
burn_in = 2000
num_iter = 10000 + burn_in
thin = 10
smp_size = int((num_iter - burn_in)/thin)

e_y = [np.zeros((T + 1, N, q)) for k in range(0, smp_size)]
e_phi = [0 for k in range(0, smp_size)]
e_V = [0 for k in range(0, smp_size)]
e_beta = [np.zeros((T + 1, p*q, 1)) for k in range(0, smp_size)]
e_Sigma = [np.zeros((q, q)) for k in range(0, smp_size)]
e_D = [np.zeros((2, N)) for k in range(0, smp_size)]



########################### Initial values for MCMC estimation
phi_prev = 2.0
V_prev = 2.0
Sigma_prev = 2.0*np.identity(q)
beta_prev = np.zeros((T + 1, p*q, 1))
D_prev = np.copy(S)
delta_phi = 21.0



########################### Metropolis-within-Gibbs algorithm
# ind will vary between 0 and smp_size - 1 (i.e. k = 1, ..., K)
ind = 0 
acc_phi, acc_D = 0, np.zeros((2, N))

# Algorithm and processing time in seconds
start = timer()
for j in range(1, num_iter + 1):

  y_curr = updt_y(it = seed_v + j,
                  m_y = y,
                  M_x = x,
                  D_usage = D_prev,
                  phi_usage = phi_prev,
                  V_usage = V_prev,
                  beta_usage = beta_prev,
                  Sigma_usage = Sigma_prev,
                  M_P = P,
                  M_inv_P = inv_P,
                  Id_obs = Id_o,
                  num_obs = N_o,
                  num_mis = N_m)
   
#  beta_curr = np.copy(beta)
  beta_curr = updt_beta(it = seed_v + j,
                        init_m0 = m0,
                        init_C0 = C0,
                        M_y = y_curr,
                        M_x = x,
                        M_g = g,
                        D_usage = D_prev,
                        W_usage = W,
                        phi_usage = phi_prev,
                        V_usage = V_prev,
                        Sigma_usage = Sigma_prev)

#  V_curr = np.copy(V).item()
  V_curr = updt_V(it = seed_v + j,
                  init_a_V = a_V,
                  init_b_V = b_V,
                  init_m0 = m0,
                  init_C0 = C0,
                  M_y = y_curr,
                  M_x = x,
                  M_g = g,
                  D_usage = D_prev,
                  beta_usage = beta_curr,
                  phi_usage = phi_prev,
                  Sigma_usage = Sigma_prev,
                  W_usage = W)

#  phi_curr = np.copy(phi).item()
  phi_curr = updt_phi(it = seed_v + j,
                      init_a_phi = a_phi,
                      init_b_phi = b_phi,
                      D_usage = D_prev,
                      M_y = y_curr,
                      M_x = x,
                      beta_usage = beta_curr,
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
                          M_y = y_curr,
                          M_X = X,
                          M_G = G,
                          beta_usage = beta_curr,
                          W_usage = W,
                          phi_usage = phi_curr,
                          V_usage = V_curr,
                          D_usage = D_prev)

#  D_curr = np.copy(D)  
  D_curr = updt_D(it = seed_v + j,
                  S_obs = S,
                  M_x = x,
                  M_y = y_curr,
                  Sigma_usage = Sigma_curr,
                  beta_usage = beta_curr,
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
    e_y[ind] = y_curr
    e_V[ind] = V_curr
    e_beta[ind] = beta_curr
    e_phi[ind] = phi_curr
    e_Sigma[ind] = Sigma_curr
    e_D[ind] = D_curr
    ind += 1
  else:
    pass

  y_prev = y_curr
  V_prev = V_curr
  phi_prev = phi_curr
  beta_prev = beta_curr
  Sigma_prev = Sigma_curr
  D_prev = D_curr

  print(j)
end = timer()
print(end - start)



########################### Acceptance ratio, given in %
round((acc_phi / num_iter)*100, 2)
np.round((acc_D / num_iter)*100, 2)



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

## V
# Vectorization
vec = e_V[0]
for j in range(1, int(smp_size)):
  now = e_V[j]
  vec = np.concatenate((vec, now), axis = None)

# Multiple line plot
plt.plot(np.arange(0, int(smp_size)),
         vec,
         linewidth = 1)
plt.xlabel('Iteration')
plt.ylabel('Sampled value',
           multialignment = 'center')
plt.title(r'$V$')

# Histogram
plt.hist(x = vec,
         density = True,
         rwidth = 0.9)
plt.ylabel('Probability')
plt.title(r'$V$')
plt.axvline(x = V,
            color = 'r',
            linestyle = '-')

## Sigma
# Vectorization
r, s = 0, 1
vec = e_Sigma[0][r,s]
for j in range(1, int(smp_size)):
  now = e_Sigma[j][r,s]
  vec = np.concatenate((vec, now), axis = None)

# Multiple line plot
plt.plot(np.arange(0, int(smp_size)),
         vec,
         linewidth = 1)
plt.xlabel('Iteration')
plt.ylabel('Sampled value',
           multialignment = 'center')
plt.title(r'$\Sigma$')

# Histogram
plt.hist(x = vec,
         density = True,
         rwidth = 0.9)
plt.ylabel('Probability')
plt.title(r'$\Sigma$')
plt.axvline(x = Sigma[r,s],
            color = 'r',
            linestyle = '-')

## V*Sigma
# Vectorization
r, s = 0, 0
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
r, s = 0, 0
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
plt.title(r'$V \cdot \phi \cdot \Sigma$')

# Histogram
plt.hist(x = vec,
         density = True,
         rwidth = 0.9)
plt.ylabel('Probability')
plt.title(r'$V \cdot \Sigma \cdot \phi$')
plt.axvline(x = V*Sigma[r,s]*phi,
            color = 'r',
            linestyle = '-')

## Some missing values
# Vectorization
t, n, i = 1, 3, 0 # 1, 5, 0
vec = e_Y[0][t][n,i]
for j in range(1, int(smp_size)):
  now = e_Y[j][t][n,i]
  vec = np.concatenate((vec, now), axis = None)

# Multiple line plot
plt.plot(np.arange(0, int(smp_size)),
         vec,
         linewidth = 1)
plt.xlabel('Iteration')
plt.ylabel('Sampled value',
           multialignment = 'center')
plt.title(r'$V \cdot \phi \cdot \Sigma$')

# Histogram
plt.hist(x = vec,
         density = True,
         rwidth = 0.9)
plt.ylabel('Probability')
plt.title(r'$V \cdot \Sigma \cdot \phi$')
plt.axvline(x = Y[t][n,i],
            color = 'r',
            linestyle = '-')
