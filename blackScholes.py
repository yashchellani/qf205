import numpy as np
from scipy.linalg import solve
from scipy.stats import norm


# Initial setup
S0 = 20         # current price
K = 10          # exercise price
T = 0.25        # expiry time
r = 0.1         # no-risk interest rate
sigma = 0.4     # volatility of underlying asset

# potentially preset
N = 2000        # number of time steps 
M = 200         # number of space grids

# ensure S0 is well within the range
S_max = max(4 * K, S0 * 2)



def calculate_ftcs(S0, K, T, r, sigma):

    # potentially preset
    N = 2000        # number of time steps 
    M = 200         # number of space grids

    # ensure S0 is well within the range
    S_max = max(4 * K, S0 * 2)

    dt = T / N
    s = np.linspace(0, S_max, M+1)  # stock price range
    C = np.clip(s - K, 0, None)      # initial condition
    index = np.arange(1, M)
    
    for n in range(N):
        C[1:-1] = (0.5 * (sigma**2 * index**2 * dt - r * index * dt) * C[:-2] +
                   (1 - sigma**2 * index**2 * dt - r * dt) * C[1:-1] +
                   0.5 * (sigma**2 * index**2 * dt + r * index * dt) * C[2:])
    
    s_idx = np.searchsorted(s, S0)  # Find the index closest to S0
    return C[s_idx]  # Return the option price at S0

def calculate_crank_nicolson(S0, K, T, r, sigma):

    N = 2000
    M = 200        
    S_max = max(4 * K, S0 * 2)

    dt = T / N
    s = np.linspace(0, S_max, M+1)
    C = np.clip(s - K, 0, None)
    index = np.arange(1, M)

    alpha = dt / 4 * (r * index - sigma**2 * index**2)
    beta = dt / 2 * (r + sigma**2 * index**2)
    gamma = -dt / 4 * (r * index + sigma**2 * index**2)
    
    A = np.diag(1 + beta) + np.diag(gamma[:-1], k=1) + np.diag(alpha[1:], k=-1)
    
    for t in range(N):
        b = np.dot(np.diag(1 - beta) + np.diag(-gamma[:-1], k=1) + np.diag(-alpha[1:], k=-1), C[1:-1])
        b[-1] += -2 * gamma[-1] * (S_max - K)
        C[1:-1] = solve(A, b)
    
    s_idx = np.searchsorted(s, S0)  # Find the index closest to S0
    return C[s_idx]  # Return the option price at S0

def calculate_exact_black_scholes(S0, K, T, r, sigma):
    d1 = (np.log(S0 / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    C = S0 * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    return max(C, 0)

# Display results
if __name__ == '__main__':
    print("FTCS Option Price at S0:", calculate_ftcs(S0, K, T, r, sigma))
    print("Crank-Nicolson Option Price at S0:", calculate_crank_nicolson(S0, K, T, r, sigma))
    print("Exact Black-Scholes Option Price at S0:", calculate_exact_black_scholes(S0, K, T, r, sigma))