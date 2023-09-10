import numpy as np

'''
S: j = 1, 2, ... , N + 1
t: k = 1, 2, ... , M + 1
'''

'''
titta på indexering
'''

def payoff_function(S, K):
    return [max(value - K, 0) for value in S]

def callfunc(K, r = 0.1, sigma = 0.25, T = 0.5, gamma = 1, N = 10, M = 10):
    s_max = 4*K
    delta_t = T/N
    delta_s = s_max/M
    grid = np.zeros((M + 1, N + 1))
    S = np.linspace(0, s_max, N + 1)
    end_vals = (payoff_function(S, K))
    #initial conds
    grid[0, :] = end_vals
    #boundary conds
    grid[:, 0] = np.zeros(M + 1)
    grid[:, -1] = list(reversed([s_max - K*np.exp(-r*(T - delta_t*(n))) for n in range(M + 1)]))
    
    
    for k in range(1, M + 1):
        A = np.zeros((N-1, N-1))
        
        #diagonal
        for j in range(2, N + 1):
            alpha = ((sigma**2)*(S[j]**2)*delta_t)/(2*(delta_s**2))
            beta = r*S[j]*delta_t/(2*delta_s)
            d = 1 + r*delta_t + 2*alpha
            A[j-2, j-2] = d
        #lowerdiagonal
        for j in range(3, N + 1):
            alpha = ((sigma**2)*(S[j]**2)*delta_t)/(2*(delta_s**2))
            beta = r*S[j]*delta_t/(2*delta_s)
            l = - alpha + beta
            A[j-2, j-3] = l
        #upperdiagonal
        for j in range(2, N):
            alpha = ((sigma**2)*(S[j]**2)*delta_t)/(2*(delta_s**2))
            beta = r*S[j]*delta_t/(2*delta_s)
            u = - alpha - beta
            A[j-2, j-1] = u

        term = np.zeros(N - 1)
        #term[0] = 0
        alpha = ((sigma**2)*(S[N]**2)*delta_t)/(2*(delta_s**2))
        beta = r*S[N]*delta_t/(2*delta_s)
        u = - alpha - beta
        rhs = grid[k - 1, 1:-1] + term
        
        next_grid_part = np.linalg.solve(A, rhs)
        grid[k, 1:-1] = next_grid_part
        
    print(grid[-1, N//4])

def bsexact(sigma: float, R: float, K: float, T: float, s: float):
    from numpy import exp, log, sqrt
    from scipy.special import erf
    
    d1 = (log(s/K)+(R+0.5*sigma**2)*T)/(sigma*sqrt(T))
    d2 = d1-sigma*sqrt(T)
    F = 0.5*s*(1+erf(d1/sqrt(2)))-exp(-R*T)*K*0.5*(1+erf(d2/sqrt(2)))
    return F

def main():
    callfunc(K = 15)

if __name__ == "__main__":
    main()
