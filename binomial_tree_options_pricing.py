import math
from option_types import Option



# test parameters - can remove after the UI is implemented
S0 = 100      # initial stock price
K = 100       # strike price
T = 1         # time to maturity in years
r = 0.06      # annual risk-free rate
N = 3         # time steps
u = 1.1       # up-factor
d = 1/u       # ensure recombining tree
opttype = Option.PUT # Option Type 'C' or 'P'


def binomial_tree(K: int, T: int, S0: int, r: float, N: int, u: float, d: float, opttype) -> float | str: 
    
    dt = T / N
    q = (math.exp(r * dt) - d) / (u - d)
    disc = math.exp(-r * dt)

    # Initialise stock prices at maturity
    S = [0] * (N + 1)
    for j in range(0, N + 1):
        S[j] = S0 * u ** j * d ** (N - j)

    # Calculate option payoff
    C = [0] * (N + 1)
    for j in range(0, N + 1):
        if opttype == Option.PUT:
            C[j] = max(0, K - S[j])
        elif opttype == Option.CALL:
            C[j] = max(0, S[j] - K)
        else:
            return 'Invalid Option Type'

    # Backward traversal through the tree
    for i in range(N - 1, -1, -1):
        for j in range(0, i + 1):
            S_current = S0 * u ** j * d ** (i - j)
            C[j] = disc * (q * C[j + 1] + (1 - q) * C[j])
            if opttype == Option.PUT:
                C[j] = max(C[j], K - S_current)
            elif opttype == Option.CALL:
                C[j] = max(C[j], S_current - K)
            else:
                return 'Invalid Option Type'

    return C[0]

res = binomial_tree(K,T,S0,r,N,u,d, opttype)
print(f'Option Price: ${res:.2f}')
