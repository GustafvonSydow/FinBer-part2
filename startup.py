import numpy as np

def payoff_function(S, K):
    return [max(value - K, 0) for value in S]

def timestep(V, S, N, M, K, r, sigma, T, gamma, tau):
    s_max = 4*K
    delta_s = s_max/M
    delta_tau = T/N
    A = np.zeros((N - 1, N - 1))
    
    for j in range(2, N + 1):
        alpha = sigma**2 * S[j]**(2*gamma) * delta_tau/(2 * delta_s **2)
        d = 1 - r*delta_tau - 2*alpha
        A[j - 2, j - 2] = d
    
    for j in range(3, N + 1):
        alpha = sigma**2 * S[j]**(2*gamma) * delta_tau/(2 * delta_s **2)
        beta = r * S[j]*(delta_tau/(2* delta_s))
        l = alpha - beta
        A[j - 2, j - 3] = l
        
    for j in range(3, N + 1):
        alpha = sigma**2 * S[j]**(2*gamma) * delta_tau/(2 * delta_s **2)
        beta = r * S[j]*(delta_tau/(2* delta_s))
        u = alpha + beta
        A[j - 3, j - 2] = u
    
    alpha = sigma**2 * S[2]**(2*gamma) * delta_tau/(2 * delta_s **2)
    beta = r * S[2]*(delta_tau/(2* delta_s))
    l = alpha - beta
    
    alpha = sigma**2 * S[N]**(2*gamma) * delta_tau/(2 * delta_s **2)
    beta = r * S[N]*(delta_tau/(2* delta_s))
    u = alpha + beta
    
    boundary = np.zeros(M - 1)
    boundary[0] = l*0
    boundary[-1] = s_max - K*(np.e**(-r*(tau)))
    return np.matmul(A, V) + boundary
        
def callfunc(N = 250, M = 50, K = 15, r = 0.1, sigma = 0.25, T = 0.5, gamma = 1):
    s_max = 4*K
    delta_s = s_max/M
    S = [j*delta_s for j in range(0, M + 1)]
    V = payoff_function(S, K)
    V = V[1: -1]
    
    for n in range(0, N + 1):
        V = timestep(V, S, N, M, K, r, sigma, T, gamma, tau = n)
        
    final_V = np.zeros(M + 1)
    final_V[0] = 0
    final_V[-1] = s_max - K*(np.e**(-r*(T)))
    final_V[1:-1] = V
    
    
    print(final_V)
    print([])

def bsexact(sigma: float, R: float, K: float, T: float, s: float):
    from numpy import exp, log, sqrt
    from scipy.special import erf
    
    d1 = (log(s/K)+(R+0.5*sigma**2)*T)/(sigma*sqrt(T))
    d2 = d1-sigma*sqrt(T)
    F = 0.5*s*(1+erf(d1/sqrt(2)))-exp(-R*T)*K*0.5*(1+erf(d2/sqrt(2)))
    return F

    


def main():
    callfunc()

if __name__ == "__main__":
    main()