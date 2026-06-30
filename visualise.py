import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("out.csv", header=None, names=["N", "call_price"])

black_scholes_value = 0.678818  # replace this with your BSM output

plt.figure(figsize=(10, 6))
plt.plot(df["N"], df["call_price"], marker='o', markersize=3, label="Monte Carlo Estimate")
plt.axhline(y=black_scholes_value, color='r', linestyle='--', label="Black-Scholes (analytical)")

plt.xscale("log")
plt.xlabel("Number of Simulations (N)")
plt.ylabel("Call Option Price")
plt.title("Monte Carlo Convergence to Black-Scholes Price")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.3)

plt.savefig("convergence.png", dpi=150)
plt.show()