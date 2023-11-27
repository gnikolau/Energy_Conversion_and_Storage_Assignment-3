# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 15:11:52 2023

@author: Georgios Nikolaou, Viktor G. Johnsen, Prathamesh Panbude, Zhengxing Zhu
"""
from sympy import symbols, diff

# Define the symbol T_cover (cover temperature)
T_cover = symbols('T_cover')

#Constants of Newton Raphson
T_ptes = 80 + 273.15  #temperature of the water pit heat storage in ᵒC
T_a = 7.4 + 273.15  #temperature of ambient air in ᵒC
T_sky = T_a - 15  #sky temperature in ᵒC
K_insu = 0.045  # Thermal conductivity of PE insulation in W/m·K
t_top = 0.25  # Thickness of the top cover in meters
epsilon_cover = 0.3  # Emissivity of the cover
sigma = 5.6697e-8  # Stefan-Boltzmann constant in W/m²·K⁴
h_convection = 5.8  # Convective heat transfer coefficient in W/m²·K

# Define the balance function f(T_cover)
f_T_cover = ((K_insu / t_top) * T_ptes + h_convection * T_a - epsilon_cover * sigma * (T_cover + T_sky) *
(T_cover**2 + T_sky**2) * (T_cover - T_sky)) / (K_insu / t_top + h_convection) - T_cover

# Derivative of f(T_cover) with respect to T_cover
f_prime_T_cover = diff(f_T_cover, T_cover)

# Initial guess for T_cover in Celsius (could be the ambient temperature)
T_cover_old = T_a 

# Newton-Raphson iteration for T_cover
for i in range(1, 101):  # Allow up to 100 iterations for convergence
    f_T_cover_val = f_T_cover.subs(T_cover, T_cover_old).evalf()
    f_prime_T_cover_val = f_prime_T_cover.subs(T_cover, T_cover_old).evalf()
    
    # Update the estimate for T_cover
    T_cover_new = T_cover_old - f_T_cover_val / f_prime_T_cover_val
    
    # Check for convergence within a tolerance
    if abs(T_cover_new - T_cover_old) < 1e-5:
        break
    
    # The new estimated T_cover in Kelvin
    T_cover_old = T_cover_new   

print("the estimation of the root is:", T_cover_new - 273.15, "ᵒC")
print("f(root) is:", ((K_insu / t_top) * T_ptes + h_convection * T_a - epsilon_cover * sigma * (T_cover_new + T_sky) *
(T_cover_new**2 + T_sky**2) * (T_cover_new - T_sky)) / (K_insu / t_top + h_convection) - T_cover_new, "ᵒC")
print("Number of iterations:", i)