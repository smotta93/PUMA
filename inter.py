import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import numpy as np

# Load the data
file_path = './prob/a/mutau/prob_nue_file2.dat'


# Função para carregar dados de um arquivo, incluindo o valor do ângulo
data = pd.read_csv(file_path, delim_whitespace=True, header=None, names=["Energy", "Probability", "test_value"])

# Extract energy and probability
x = data['Energy']
y = data['Probability']

# Plot the original data
plt.figure(figsize=(10, 6))
plt.scatter(x, y, color='blue', label='Original Data', alpha=0.7, s=10)
plt.plot(x, y, color="black", label="Original Data", linewidth=2)
plt.xlabel('Energy (GeV)')
plt.ylabel('Probability')
plt.title('Oscillation Probability')
plt.grid(True)
plt.legend()
plt.show()

# Perform a non-linear interpolation
spl = UnivariateSpline(x, y, s=0.02)  # Adjust smoothness parameter 's' for a smoother curve
x_smooth = np.linspace(x.min(), x.max(), 500)
y_smooth = spl(x_smooth)

# Plot the original data and the smoothed curve
plt.figure(figsize=(10, 6))
plt.scatter(x, y, color='blue', label='Original Data', alpha=0.7, s=10)
plt.plot(x_smooth, y_smooth, color='red', label='Smoothed Curve', linewidth=2)
plt.xlabel('Energy (GeV)')
plt.ylabel('Probability')
plt.title('Smoothed Oscillation Probability Curve')
plt.grid(True)
plt.legend()
plt.show()

# Approximate equation (linear least-squares fit for representation)
from numpy.polynomial.polynomial import Polynomial

# Fit a simple polynomial (linear here for the line approximation)
poly_coeff = Polynomial.fit(x, y, 1).convert().coef
equation_text = f"y = {poly_coeff[1]:.4f}x + {poly_coeff[0]:.4f}"

# Plot with the linear approximation
plt.figure(figsize=(10, 6))
plt.scatter(x, y, color='blue', label='Original Data', alpha=0.7, s=10)
plt.plot(x_smooth, y_smooth, color='red', label='Smoothed Curve', linewidth=2)
plt.plot(x, poly_coeff[1] * x + poly_coeff[0], 'g--', label=f'Approximation: {equation_text}')
plt.xlabel('Energy (GeV)')
plt.ylabel('Probability')
plt.title('Smoothed Oscillation Probability with Linear Approximation')
plt.grid(True)
plt.legend()
plt.show()
