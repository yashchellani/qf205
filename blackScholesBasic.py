from scipy.stats import norm
import numpy as np
from option_types import Option

# Function to calculate d1 and d2
def d1(S, K, T, r, sigma):
    return (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))

def d2(S, K, T, r, sigma):
    return d1(S, K, T, r, sigma) - sigma * np.sqrt(T)

# Greeks calculations
def delta_call(S, K, T, r, sigma):
    return norm.cdf(d1(S, K, T, r, sigma))

def delta_put(S, K, T, r, sigma):
    return norm.cdf(d1(S, K, T, r, sigma)) - 1

def gamma(S, K, T, r, sigma):
    return norm.pdf(d1(S, K, T, r, sigma)) / (S * sigma * np.sqrt(T))

def theta_call(S, K, T, r, sigma):
    d1_val = d1(S, K, T, r, sigma)
    d2_val = d2(S, K, T, r, sigma)
    term1 = - (S * norm.pdf(d1_val) * sigma) / (2 * np.sqrt(T))
    term2 = r * K * np.exp(-r * T) * norm.cdf(d2_val)
    return term1 - term2

def theta_put(S, K, T, r, sigma):
    d1_val = d1(S, K, T, r, sigma)
    d2_val = d2(S, K, T, r, sigma)
    term1 = - (S * norm.pdf(d1_val) * sigma) / (2 * np.sqrt(T))
    term2 = r * K * np.exp(-r * T) * norm.cdf(-d2_val)
    return term1 + term2

def vega(S, K, T, r, sigma):
    return S * norm.pdf(d1(S, K, T, r, sigma)) * np.sqrt(T)

def rho_call(S, K, T, r, sigma):
    return K * T * np.exp(-r * T) * norm.cdf(d2(S, K, T, r, sigma))

def rho_put(S, K, T, r, sigma):
    return -K * T * np.exp(-r * T) * norm.cdf(-d2(S, K, T, r, sigma))

def black_scholes_call(S, K, T, r, sigma):
    d1_val = d1(S, K, T, r, sigma)
    d2_val = d2(S, K, T, r, sigma)
    call_price = (S * norm.cdf(d1_val) - K * np.exp(-r * T) * norm.cdf(d2_val))
    return call_price

def black_scholes_put(S, K, T, r, sigma):
    d1_val = d1(S, K, T, r, sigma)
    d2_val = d2(S, K, T, r, sigma)
    put_price = (K * np.exp(-r * T) * norm.cdf(-d2_val) - S * norm.cdf(-d1_val))
    return put_price

# Example parameters for calculations

T = 0.25   # Time to maturity in years
r = 0.1  # Risk-free interest rate
sigma = 0.4  # Volatility
K = 10.  # Strike price
S = 4 * K  # Underlying asset price

# Calculate Greeks for a call option
delta_c = delta_call(S, K, T, r, sigma)
gamma_c = gamma(S, K, T, r, sigma)
theta_c = theta_call(S, K, T, r, sigma) / 365  # Annualize theta
vega_c = vega(S, K, T, r, sigma) / 100  # For a 1% change in volatility
rho_c = rho_call(S, K, T, r, sigma) / 100  # For a 1% change in interest rate

# Calculate Greeks for a put option
delta_p = delta_put(S, K, T, r, sigma)
gamma_p = gamma(S, K, T, r, sigma)
theta_p = theta_put(S, K, T, r, sigma) / 365  # Annualize theta
vega_p = vega(S, K, T, r, sigma) / 100  # For a 1% change in volatility
rho_p = rho_put(S, K, T, r, sigma) / 100  # For a 1% change in interest rate

call_price = black_scholes_call(S, K, T, r, sigma)
put_price = black_scholes_put(S, K, T, r, sigma)



print(delta_c, gamma_c, theta_c, vega_c, rho_c)
print(delta_p, gamma_p, theta_p, vega_p, rho_p)
print(call_price, put_price)