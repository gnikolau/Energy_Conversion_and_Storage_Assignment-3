# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 15:11:52 2023

@author: Georgios Nikolaou, Viktor G. Johnsen, Prathamesh Panbude, Zhengxing Zhu
"""
from sympy import symbols, diff, sqrt
import math

# Define the symbol L_bottom
L_bottom = symbols("L_bottom")

# Constants of Newton Raphson
H = 16  # The water depth of the PTES (m)
V_PTES = 64146.79  # Volume of PTES (m**3)
angle_degrees = 26.6  # The slope of side of the PTES (degrees)
L_bottom_old = 50 # Initial guess for L_bottom (m)

# Convert angle to radians
angle_radians = math.radians(angle_degrees)

# Calculate tan of the angle
tan_angle = math.tan(angle_radians)

# Define the balance function f(L_bottom)
f_L_bottom = V_PTES - H/3 * (L_bottom**2 + (L_bottom + 2*H/tan_angle)**2 + 
             sqrt((L_bottom + 2*H/tan_angle)**2 * L_bottom**2))

# Derivative of f(L_bottom) with respect to L_bottom
f_prime_L_bottom = diff(f_L_bottom, L_bottom)

# Newton-Raphson iteration for L_bottom
for i in range(1, 101):  # Allow up to 100 iterations for convergence
    f_L_bottom_val = f_L_bottom.subs(L_bottom, L_bottom_old).evalf()
    f_prime_L_bottom_val = f_prime_L_bottom.subs(L_bottom, L_bottom_old).evalf()
    
    # Update the estimate for L_bottom
    L_bottom_new = L_bottom_old - f_L_bottom_val / f_prime_L_bottom_val
    
    # Check for convergence within a tolerance
    if abs(L_bottom_new - L_bottom_old) < 1e-5:
        break
    
    # The new estimated L_bottom
    L_bottom_old = L_bottom_new   

print("Number of iterations:", i)
print("the estimation of the root is:", L_bottom_new, "m")
print("f(root) is:", V_PTES - H/3 * (L_bottom_new**2 + (L_bottom_new + 2*H/tan_angle)**2 + 
      sqrt((L_bottom_new + 2*H/tan_angle)**2 * L_bottom_new**2,).evalf()), "m**3")