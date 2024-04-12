import numpy as np
from math import exp, sqrt

# Example usage
S0 = 100  # Initial stock price
K = 100   # Strike price
T = 1     # Time to maturity (years)
r = 0.06  # Risk-free interest rate
sigma = 0.2  # Volatility
num_simulations = 100  # Increase for better accuracy

def generate_stock_paths(S0, mu, sigma, T, dt, num_simulations):
    """
    Generates multiple stock price paths using Geometric Brownian Motion.

    Args:
        S0 (float): Initial stock price.
        mu (float): Expected drift (annualized return).
        sigma (float): Volatility (annualized standard deviation).
        T (float): Time to maturity of the option (in years).
        dt (float): Time increment for simulation.
        num_simulations (int): Number of simulated paths.

    Returns:
        numpy.ndarray: Array of simulated stock price paths.
    """

    num_steps = int(T / dt)
    paths = np.zeros((num_simulations, num_steps + 1))
    paths[:, 0] = S0

    rand = np.random.standard_normal((num_simulations, num_steps))  
    for t in range(1, num_steps + 1):
        paths[:, t] = paths[:, t - 1] * np.exp((mu - 0.5 * sigma ** 2) * dt + sigma * sqrt(dt) * rand[:, t - 1])

    return paths

def monte_carlo_option_price(S0, K, T, r, sigma, num_simulations):
    """
    Prices a European call option using Monte Carlo simulation.

    Args:
        S0 (float): Initial stock price.
        K (float): Strike price of the option.
        T (float): Time to maturity of the option (in years).
        r (float): Risk-free interest rate (annualized).
        sigma (float): Volatility (annualized standard deviation).
        num_simulations (int): Number of simulated paths.

    Returns:
        float: Estimated price of the call option.
    """

    dt = T / num_simulations  # Smaller time steps for better accuracy
    paths = generate_stock_paths(S0, r, sigma, T, dt, num_simulations)

    # Payoffs at maturity
    payoffs = np.maximum(paths[:, -1] - K, 0) 

    # Discounting back to present value and averaging
    option_price = exp(-r * T) * payoffs.mean()  

    return option_price



call_price = monte_carlo_option_price(S0, K, T, r, sigma, num_simulations)
print("Estimated price of the European call option:", call_price)