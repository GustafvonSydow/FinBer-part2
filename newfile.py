import numpy as np

def payoff_function(S, K):
    return [max(value - K, 0) for value in S]

def callfunc(K = 15, r = 0.1, sigma = 0.25, T = 0.5, gamma = 0.8, M = 200, N = 1000):
    s_max = 4*K
    delta_s = s_max/M
    delta_t = T/N
    S = [j*delta_s for j in range(0, M + 1)]
    V = payoff_function(S, K)
    print(V)
    
    for n in reversed(range(1, N + 1)):
        t_n_minus_1 = (n - 1)*delta_t
        V[0] = 0
        V[-1] = s_max - K*np.exp(-r*(T - t_n_minus_1))
        
        for j in range(1, M):
            V[j] = V[j] + r * S[j]*(delta_t/(2*delta_s))*(V[j + 1] - V[j - 1]) + ((sigma**2)/2)*(S[j]**(2*gamma))*(delta_t/(delta_s**2))*(V[j + 1] - 2*V[j] + V[j - 1]) - delta_t*r*V[j]
            
    return V

def main():
    V = callfunc()
    print(V)
    
if __name__ == "__main__":
    main()