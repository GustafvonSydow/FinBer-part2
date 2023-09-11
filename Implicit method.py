import numpy as np
import matplotlib.pyplot as plt

'''
S: j = 0, 1, 2, ... , N
t: k = 0, 1, 2, ... , M
'''

def payoff_function(S, K):
    return [max(value - K, 0) for value in S]

def callfunc(K = 15, r = 0.1, sigma = 0.25, T = 0.5, gamma = 1, N = 1000, M = 1000):
    s_max = 4*K
    delta_t = T/M
    delta_s = s_max/N
    grid = np.zeros((M + 1, N + 1))
    S = np.linspace(0, s_max, N + 1)
    end_vals = (payoff_function(S, K))
    #initial conds
    grid[0, :] = end_vals
    #boundary conds
    grid[:, 0] = np.zeros(M + 1)
    grid[:, -1] = list(reversed([s_max - K*np.exp(-r*(T - delta_t*(n))) for n in range(M + 1)]))
    
    
    for k in range(0, M):
        A = np.zeros((N-1, N-1))
        
        #diagonal
        for j in range(1, N):
            alpha = ((sigma**2)*(S[j]**(2*gamma))*delta_t)/(2*(delta_s**2))
            beta = r*S[j]*delta_t/(2*delta_s)
            d = 1 + r*delta_t + 2*alpha
            A[j-1, j-1] = d
        #lowerdiagonal
        for j in range(2, N):
            alpha = ((sigma**2)*(S[j]**(2*gamma))*delta_t)/(2*(delta_s**2))
            beta = r*S[j]*delta_t/(2*delta_s)
            l = - alpha + beta
            A[j-1, j-2] = l
        #upperdiagonal
        for j in range(1, N - 1):
            alpha = ((sigma**2)*(S[j]**(2*gamma))*delta_t)/(2*(delta_s**2))
            beta = r*S[j]*delta_t/(2*delta_s)
            u = - alpha - beta  
            A[j-1, j] = u

        term = np.zeros(N - 1)
        #term[0] = 0
        alpha = ((sigma**2)*(S[-2]**2)*delta_t)/(2*(delta_s**2))
        beta = r*S[-2]*delta_t/(2*delta_s)
        u = - alpha - beta
        term[-1] = - u*grid[k + 1, -1]
        rhs = grid[k, 1:-1] + term
        
        next_grid_part = np.linalg.solve(A, rhs)
        grid[k + 1, 1:-1] = next_grid_part
    #hard coded    
    return grid[-1, N//4]

def bsexact(sigma: float, R: float, K: float, T: float, s: float):
    from numpy import exp, log, sqrt
    from scipy.special import erf
    
    d1 = (log(s/K)+(R+0.5*sigma**2)*T)/(sigma*sqrt(T))
    d2 = d1-sigma*sqrt(T)
    F = 0.5*s*(1+erf(d1/sqrt(2)))-exp(-R*T)*K*0.5*(1+erf(d2/sqrt(2)))
    return F

def delta_t_test():
    s = 15
    F = bsexact(sigma = 0.25, R = 0.1, K = 15, T = 0.5, s = s)
    M_array = [(2**exp) for exp in range(10)]
    delta_t_array = [0.5/M for M in M_array]
    E_array = np.zeros(len(M_array))
    
    for index, M in enumerate(M_array):
        print(np.abs(callfunc(M = M) - F))
        E_array[index] = np.abs(callfunc(M = M) - F)
        
    expected_errors_array = [delta_t for delta_t in delta_t_array]
    plt.loglog(delta_t_array, E_array, 'bo', label = 'Estimated errors')
    plt.loglog(delta_t_array, expected_errors_array, 'ro', label = 'Expected errors')
    plt.xlabel('delta t')
    plt.ylabel('error')
    plt.title('Error convergence for delta t')
    plt.savefig('Error convergence for delta t.jpg')
    plt.legend()
    plt.show()
    
def delta_s_test():
    s = 15
    F = bsexact(sigma = 0.25, R = 0.1, K = 15, T = 0.5, s = s)
    N_array = [8*(2**exp) for exp in range(7)]
    delta_t_array = [60/N for N in N_array]
    E_array = np.zeros(len(N_array))
    
    for index, N in enumerate(N_array):
        print(np.abs(callfunc(N = N) - F))
        E_array[index] = np.abs(callfunc(N = N) - F)
        
    expected_errors_array = [delta_t**2 for delta_t in delta_t_array]
    plt.loglog(delta_t_array, E_array, 'bo', label = 'Estimated errors')
    plt.loglog(delta_t_array, expected_errors_array, 'ro', label = 'Expected errors')
    plt.xlabel('delta s')
    plt.ylabel('error')
    plt.title('Error convergence for delta s')
    plt.savefig('Error convergence for delta s.jpg')
    plt.legend()
    plt.show()

def main():
    delta_s_test()

if __name__ == "__main__":
    main()
