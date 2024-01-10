'''
Practical Implementation of the Lifting Line Theory
Author: Lorenzo Cappello
Date: 10.01.2024

This tool computes lift and drag for a finite wing using the Lifting Line Theory and airfoil data generated with X-Foil for a NACA-2412 airfoil

'''

import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
from scipy import stats
import tkinter as tk
from tkinter import*
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy import integrate
import scipy
import scipy.special as sp
from scipy.integrate import simps

'------------------------------ Atmospheric Module -----------------------------------------------------'

# This module calculates atmospheric parameters according to ISA

rho0 = 1.225                    # ISA density at sea level [kg/m^3]
h = 1                           # altitude [km]
T_ref = 288.15                  # ISA temperature at sea level [K]
P_ref = 101325                  # ISA pressure at sea level [Pa]

if h <= 11:
    T = T_ref - 6.5*h           # temperature at reference altitude [K] 
    P = P_ref * (T/T_ref)**5.2561   # pressure at reference altitude [Pa]
elif 11 < h <= 20:
    T = 216.65                  # temperature at reference altitude [K] 
    P = 22630.6*10**(-((h-11)/14.596)) # pressure at reference altitude [Pa]
elif 20 < h < 32:
    T = 216.65 - (h-20)         # temperature at reference altitude [K] 
    P = P = 22630.6*10**-((h-11)/14.596) # pressure at reference altitude [Pa]
else:
    print("Please input an altitude between 0 and 32 kilometers")

mu_ref = 1.711*10**-5           # air dynamic viscosity at standard conditions [kg/ms]
T0 = 273.15 # 0 degrees Celsius [K]
mu = mu_ref*(T/T0)**(1.5)*((T0+110.4)/(T+110.4)) # air dynamic viscosity at reference conditions (Sutherland's relation) [kg/ms]
rho = rho0*(P/P_ref)*(T_ref/T)  # density at reference conditions [kg/m^3]

'------------------------------------------------- Input GUI ---------------------------------------------------'

# This module takes user inputs from a GUI

v_inf = 60                      # free stream velocity default value
Re = 1200000                    # Reynolds number defalut value
AoA = 6                         # default wing angle of attack
Data_selection = "XFoil"        # default data selection
c_root = 0.76                   # default chord length at wing root
c_tip = 0.29                    # default chord length at wing tip
b = 20                          # default wing span
n = 49                          # number of sections

def on_submit():
    global v_inf, Re, AoA, Data_selection, c_root, c_tip, b, n
    v_inf = int(v_inf_entry.get())
    Re = int(Re_entry.get())
    AoA = int(AoA_entry.get())
    Data_selection = data_selection_var.get()
    c_root = float(c_root_entry.get())
    c_tip = float(c_tip_entry.get())
    b = float(b_entry.get())
    n = int(n_entry.get())

window = tk.Tk()
window.title("Input Collection")

v_inf_label = tk.Label(window, text="Free Stream Velocity (m/s):")
v_inf_label.pack()
v_inf_entry = tk.Entry(window)
v_inf_entry.pack()

Re_label = tk.Label(window, text="Reynolds Number:")
Re_label.pack()
Re_entry = tk.Entry(window)
Re_entry.pack()

AoA_label = tk.Label(window, text="Angle of Attack:")
AoA_label.pack()
AoA_entry = tk.Entry(window)
AoA_entry.pack()

data_selection_label = tk.Label(window, text="Data Source (XFoil or Windtunnel):")
data_selection_label.pack()
data_selection_var = tk.StringVar()
data_selection_entry = tk.Entry(window, textvariable=data_selection_var)
data_selection_entry.pack()

c_root_label = tk.Label(window, text="Chord Length at Wing Root (m):")
c_root_label.pack()
c_root_entry = tk.Entry(window)
c_root_entry.pack()

c_tip_label = tk.Label(window, text="Chord Length at Wing Tip (m):")
c_tip_label.pack()
c_tip_entry = tk.Entry(window)
c_tip_entry.pack()

b_label = tk.Label(window, text="Wing Span (m):")
b_label.pack()
b_entry = tk.Entry(window)
b_entry.pack()

n_label = tk.Label(window, text="Number of spanwise sections: ")
n_label.pack()
n_entry = tk.Entry(window)
n_entry.pack()

submit_button = tk.Button(window, text="Submit", command=on_submit)
submit_button.pack()

window.mainloop()

b_half = b/2                        # span of only one wing [m]
TR = c_tip/c_root                   # taper ratio
Aw = b_half*(c_root+c_tip)          # wing area
AR = b**2/(2*Aw)                    # aspect ratio

path = r"C:\Users\Lorenzo\OneDrive - ZHAW\Studium\ZHAW Master of Science in Engineering\HS23\VT1\Model" # absolute file path

def dataextraction(path,file_name):
    data = {
    'alpha': [],
    'CL': [],
    'CD': [],
    'CDp': [],
    'CM': [], 
    'Top_Xtr': [],
    'Bot_Xtr': []
        }
    with open(path + file_name,'r') as file:
                lines = file.readlines()
                for line in lines[12:]:
                    if not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) == 7:
                        alpha, CL, CD, CDp, CM, Top_Xtr, Bot_Xtr = map(float, parts)
                        data['alpha'].append(alpha)
                        data['CL'].append(CL)
                        data['CD'].append(CD)
                        data['CDp'].append(CDp)
                        data['CM'].append(CM)
                        data['Top_Xtr'].append(Top_Xtr)
                        data['Bot_Xtr'].append(Bot_Xtr)
    return data

'---------------------------------------- Airfoil Polar Data Extraction ---------------------------------------------------'

# This module extracts airfoil polar data from text files

while True:
    if Data_selection == "XFoil":
        path = r"C:\Users\Lorenzo\OneDrive - ZHAW\Studium\ZHAW Master of Science in Engineering\HS23\VT1\Model\XFoil"
        file_name1 = f"\Re{Re}.txt"
        file_name2 = f"Re{Re}.txt"
        location = path + f"Re{Re}.txt"
        file_list = os.listdir(path)
        if file_name2 in file_list:
            data = dataextraction(path,file_name1)
        else:
            airfoil = "myairfoils/2412geometry" 
            AoA_i = -10 # define starting AoA
            AoA_f = 20 # define ending AoA
            AoA_step = 0.1 # define step size for AoA increase
            ITER = 100 # define number of iterations performed by X-Foil
            if os.path.exists("polar.txt"):
                os.remove("polar.txt")  # output file
                input = open("input_file.in", 'w') # creation of input file
                input.write("LOAD {0}.dat\n".format(airfoil)) # load airfoil
                input.write(airfoil + '\n') # put airfoil name into input file
                input.write("PANE\n") # Set current-airfoil panel nodes ( 140 ) based on curvature
                input.write("OPER\n") # Set direct operating point(s)
                input.write("Visc {0}\n".format(Re)) # Set the Reynolds number in the input file
                input.write("PACC\n") 
                input.write("polar.txt\n\n") # write the input data in the 
                input.write("ITER {0}\n".format(ITER))
                input.write("ASeq {0} {1} {2}\n".format(AoA_i, AoA_f,
                                            AoA_step)) # define flow conditions in input file
                input.write("\n\n")
                input.write("quit\n") # quit X-Foil after run
                input.close() # close the input file

                subprocess.call("xfoil.exe < input_file.in", shell=True) # execute X-Foil run

                polar_data = np.loadtxt("polar.txt", skiprows=12) # show file with polar data from run
                file_name = "polar.txt"
                data = {
                    'alpha': [],
                    'CL': [],
                    'CD': [],
                    'CDp': [],
                    'CM': [], 
                    'Top_Xtr': [],
                    'Bot_Xtr': []
                    }
                with open(file_name,'r') as file:
                    lines = file.readlines()
                    for line in lines[12:]:
                        if not line.strip():
                            continue
                        parts = line.split()
                        if len(parts) == 7:
                            alpha, CL, CD, CDp, CM, Top_Xtr, Bot_Xtr = map(float, parts)
                            data['alpha'].append(alpha)
                            data['CL'].append(CL)
                            data['CD'].append(CD)
                            data['CDp'].append(CDp)
                            data['CM'].append(CM)
                            data['Top_Xtr'].append(Top_Xtr)
                            data['Bot_Xtr'].append(Bot_Xtr)

        alpha_foil = np.array(data['alpha'])
        CL_foil = np.array(data['CL'])
        CD_foil = np.array(data['CD'])
        CDp_foil = np.array(data['CDp'])
        CM_foil = np.array(data['CM'])
        Top_Xtr = np.array(data['Top_Xtr'])
        Bot_Xtr = np.array(data['Bot_Xtr'])
        break

    elif Data_selection == "Windtunnel":
        path = r"C:\Users\Lorenzo\OneDrive - ZHAW\Studium\ZHAW Master of Science in Engineering\HS23\VT1\Model\Windtunnel"
        file_name = path + f"\Re{Re}.txt"
        file_list = os.listdir(path)
        if file_name in file_list:
            data = dataextraction(path,file_name)
        else:
            path = r"C:\Users\Lorenzo\OneDrive - ZHAW\Studium\ZHAW Master of Science in Engineering\HS23\VT1\Model\Windtunnel"
            file_list = os.listdir(path)
            valid_files = [file for file in file_list if file.startswith("Re") and file.endswith(".txt")]
            if not valid_files: 
                print("No valid files found in the folder.")
            else: 
                Re_numbers = [int(file[2:-4]) for file in valid_files]
                for i in range(len(Re_numbers)):
                    lower_difference = np.min(np.abs(Re-Re_numbers[i]))
                    upper_difference = np.min(np.abs(Re_numbers[i]-Re))
                    lower_Re = Re - lower_difference
                    upper_Re = Re + upper_difference       
                    closest_lower = f"Re{lower_Re}.txt"
                    closest_upper = f"Re{upper_Re}.txt"

            def read_data(filename):
                data = {
                    'alpha': [],
                    'CL': [],
                    'CD': [],
                    'CDp': [],
                    'CM': []
                }
                with open(os.path.join(path, filename),'r') as file:
                    lines = file.readlines()[12:]
                    for line in lines:
                        values = line.split()
                        data['alpha'].append(float(values[0]))
                        data['CL'].append(float(values[1]))
                        data['CD'].append(float(values[2]))
                        data['CDp'].append(float(values[3]))
                        data['CM'].append(float(values[4]))
                return data
            # extract parameters from upper and lower file
            lower_data = read_data(closest_lower)
            upper_data = read_data(closest_upper)
            alpha_lower = np.array(lower_data['alpha'])
            alpha_upper = np.array(upper_data['alpha'])
            alpha_foil = np.array(lower_data['alpha'])
            cl_lower = np.array(lower_data["CL"])
            cd_lower = np.array(lower_data["CD"])
            cm_lower = np.array(lower_data['CM'])      
            cl_upper = np.array(upper_data["CL"])
            cd_upper = np.array(upper_data["CD"])
            cm_upper = np.array(upper_data['CM']) 
            if len(cl_lower) == len(cl_upper):
                # check if there are same number of lift coefficients in the files
                CL_foil =  np.interp(alpha_foil,cl_lower,cl_upper)  
            else:
                common_alpha = np.intersect1d(alpha_lower,alpha_upper)
                mask1 = np.isin(alpha_lower, common_alpha)
                cl_lower = cl_lower[mask1]
                mask2 = np.isin(alpha_upper, common_alpha)
                cl_upper = cl_upper[mask2]
                CL_foil = np.interp(common_alpha,cl_lower,cl_upper)
            if len(cd_lower) == len(cd_upper):
                # check if there are same number of drag coefficients in the files
                CD_foil = np.interp(alpha_foil, cd_lower, cd_upper)    
            else: 
                common_alpha = np.intersect1d(alpha_lower,alpha_upper)
                mask1 = np.isin(alpha_lower, common_alpha)
                cd_lower = cd_lower[mask1]
                mask2 = np.isin(alpha_upper, common_alpha)
                cd_upper = cd_upper[mask2]
                CD_foil = np.interp(common_alpha,cd_lower,cd_upper)
            if len(cm_lower) == len(cm_upper):
                # check if there are same number of moment coefficients in the files
                CM_foil = np.interp(alpha_foil, cm_lower, cm_upper)   
            else: 
                common_alpha = np.intersect1d(alpha_lower,alpha_upper)
                mask1 = np.isin(alpha_lower, common_alpha)
                cm_lower = cm_lower[mask1]
                mask2 = np.isin(alpha_upper, common_alpha)
                cm_upper = cm_upper[mask2]
                CM_foil = np.interp(common_alpha,cm_lower,cm_upper)
            break
    else: 
        print("Invalid data source selection. Please enter either 'XFoil' or 'Windtunnel' as data source.")
        break

'----------------Aerodynamic Module: Lifting Line Theory--------------------------------------'

# This module contains the practical implementation of the classical LTT as described by Bertin and Cummings (2021)

if len(alpha_foil) == len(CL_foil):
    slope, intercept, r_value, p_value, std_err = stats.linregress(alpha_foil, CL_foil) # linear interpolation
else:
    slope, intercept, r_value, p_value, std_err = stats.linregress(common_alpha, CL_foil)
a0 = 2*np.pi # or slope*np.pi/180                               # airfoil cl slope
cl_index = np.argmin(np.abs(CL_foil))
alpha0 = -1.2*np.pi/180 # or alpha_foil[cl_index]*np.pi/180     # zero lift angle of attack
index = np.where(alpha_foil ==0)[0]
cl0 =  CL_foil[index[0]]                                        # airfoil lift coefficient at zero angle of attack
cd0 = CD_foil[index[0]]                                         # airfoil drag coefficient at zero angle of attack
iw = 0                                                          # wing incidence angle at the root 
alpha_tw = 0                                                    # twist angle
alpha_eff = AoA + iw                                            # effective angle of attack at wing root 
alpha = np.linspace(alpha_eff+alpha_tw,alpha_eff,n)*np.pi/180   # local angle of attack [rad]
alpha_deg = alpha*180/np.pi                                     # local angle of attack [deg]
index2 = np.where(alpha_foil == alpha_eff)[0]
cm = CM_foil[index2[0]]                                         # airfoil moment coefficient

if n==4:
    phi = np.array([22.5,45,67.5,90])
else:
    iv = 0.1e-8
    phi = np.linspace(iv,90,n)

phi_r = phi*np.pi/180                                           # section angle [rad]
y = (b/2)*np.cos(phi_r)                                         # spanwise stations - y-coordinate [m]
c = c_root*(1+(TR-1)*np.cos(phi_r))                             # local chord length [m]
geo = (c*a0)/(2*AR*c_root*(1+TR))                               # wing geometrical characteristics in one variable at different stations
Left = geo*(alpha-alpha0)*np.sin(phi_r)                         # left side of matrix
B = np.zeros((n,n))                                             # right side of matrix

for i in range(0,n):
    for j in range(0,n):
        B[i,j] = np.sin((2*(j+1)-1)*phi_r[i]) * (geo[i]*(2*(j+1)-1)+np.sin(phi_r[i])) 
A = np.linalg.solve(B,Left)                                     # matrix inversion

CL = np.pi * AR * A[0]                                          # wing lift coefficient

delta = 0 
for i in range(1,n):
    delta = delta + (2*(i+1)-1)*A[i]**2/A[0]**2

CD_ind = (CL**2/(np.pi*AR)) * (1+delta)                         # wing induced drag coefficient
CD = cd0 + CD_ind                                               # wing total drag coefficient (neglecting wave drag and viscous drag)
L = CL*0.5*rho*v_inf**2*Aw*2                                    # Lift [N]
x = np.linspace(-b/2, b/2, n)
D_ind = CD_ind*0.5*rho*v_inf**2*Aw                              # Induced Drag [N]
D = CD*0.5*rho*v_inf**2*Aw                                      # Total Drag [N] (neglecting wave drag and viscous drag)
gamma_temp = np.zeros((n))
for i in range(0,n):
    for j in range(0,n):
        gamma_temp[i] = gamma_temp[i] + A[j]*np.sin((2*(j+1)-1)*phi_r[i])

gamma = 2*b*v_inf*gamma_temp                                    # Circulation distribution [m^2/s]
tau = 0.05
aw =  a0/(1+(a0/(np.pi*AR))*(1+tau))                            # lift slope of finite wing [rad]
aw_deg = aw*180/np.pi                                           # lift slope of finite wing [deg]
L_dist = rho*v_inf*gamma                                        # lift distribution [N]
M_sum = np.zeros_like(L_dist)
y_new = y[::-1]
L_dist_new = L_dist[::-1]                                       # re-orderd lift distribution [N]
for i in range(len(y)):
    if i == 0:
        continue  # Skip the first element (tip of the wing)
    # Sum the contributions from the current point to the wing root
    M_sum[i] = np.sum(L_dist_new[:i] * np.diff(y_new[:i+1]))
M_dist = M_sum                                                  # bending moment distribution [Nm]

temp = np.zeros(n)
for i in range(0,n):
    for j in range(0,n):
        temp[i] = temp[i] + ((2*(j+1)-1)) * A[j] * np.sin((2*(j+1)-1) * phi_r[i]) / np.sin(phi_r[i])

dw = -v_inf*temp                                                # downwash [m/s]
alpha_i = dw/v_inf                                              # induced angle of attack [rad]

D_dist = L_dist*np.abs(alpha_i)                                 # induced drag distribution [N]

# alternative method to compute the circulation 

def compute_gamma(phi_r, A, b, v_inf):
    gamma_phi = 0
    for i in range(n):
        m = 2 * i + 1
        gamma_phi += A[i] * np.sin(m * phi_r)
    return 2 * b * v_inf * gamma_phi

gamma_phi = [compute_gamma(phi_r, A, b, v_inf) for phi_r in phi_r]

#print("Circulation distribution :", gamma)
#print("Circulation distribution alternative:", gamma_phi)
print("CL = ",CL)
#print("CD = ", CD)
print("CD_ind = ",CD_ind)
#print("Lift (N) = ", L) # discrepancy between L and L_alt
#print("Alternative Lift calculation: ", L_alt)
#print("Total Drag (N) = ",D)
#print("Induced drag (N) = ", D_ind)
#print("Lift distribution = ", L_dist)
#print("CL_alpha_wing = ", aw)
#print("Downwash = ",dw)
#print("Induced Angle of Attack = ", alpha_i*180/np.pi)

# Plots 

plt.plot(y,L_dist)
plt.xlabel('y (m)',fontsize=10)
plt.ylabel(r'L (N)',fontsize=10)
plt.title('Lift Distribution',fontsize=14)
plt.grid(False)
plt.show()

plt.plot(y,M_dist)
plt.xlabel('y (m)',fontsize=10)
plt.ylabel(r'M (Nm)',fontsize=10)
plt.title('Bending Moment Distribution',fontsize=14)
plt.grid(False)
plt.show()

plt.plot(y,dw)
plt.xlabel('y (m)',fontsize=10)
plt.ylabel(r'dw (m/s)',fontsize=10)
plt.title('Downwash Distribution',fontsize=14)
plt.grid(False)
plt.show()

plt.plot(y,gamma)
plt.xlabel('y (m)',fontsize=10)
plt.ylabel(r'Gamma (m^2/s)',fontsize=10)
plt.title('Circulation Distribution',fontsize=14)
plt.grid(False)
plt.show()

plt.plot(y,D_dist)
plt.xlabel('y (m)',fontsize=10)
plt.ylabel(r'Di (N)',fontsize=10)
plt.title('Induced Drag Distribution',fontsize=14)
plt.grid(False)
plt.show()


'----------------------------------------------------------Output file of classical LTT computation results --------------------------------------------------'

# This module writes the outputs of the computations performed with the classical LTT into an output text file

output = open("output.txt", "a")
out1 = CL
out2 = CD
out3 = CD_ind
out4 = L
out5 = D
out6 = D_dist
out7 = L_dist
out8 = M_dist
out9 = alpha_deg

def format_value(val, width=10, decimal_places=3):
    if isinstance(val, float):
        formatted_val = f"{val:.{decimal_places}f}"
    elif isinstance(val, (list, np.ndarray)) and len(val) > 0:
        formatted_val = format_value(val[0], width, decimal_places)
    else:
        formatted_val = str(val)

    return formatted_val.ljust(width)


def is_file_empty(file_path):
    return os.path.exists(file_path) and os.stat(file_path).st_size == 0

def format_array(array, width=10, decimal_places=3):
    formatted_array = []
    for val in array:
        formatted_array.append(format_value(val, width, decimal_places))
    return formatted_array

file_path = "output.txt"

with open(file_path, "a") as output:

    if is_file_empty(file_path):

        headers = ["AoA", "CL", "CD", "CD_ind", "Lift", "Drag"]
        formatted_headers = [h.ljust(10) for h in headers]
        output.write("".join(formatted_headers) + "\n")

    single_value_vars = [out1, out2, out3, out4, out5, out6]

    formatted_values = [format_value(val) for val in [out9] + single_value_vars]
    output.write("".join(formatted_values) + "\n")

    output.write("\nInduced Drag Distribution:\n")
    lift_dist_formatted = format_array(D_dist)
    for val in lift_dist_formatted:
        output.write(val + "\n")

    output.write("\nLift Distribution:\n")
    lift_dist_formatted = format_array(L_dist)
    for val in lift_dist_formatted:
        output.write(val + "\n")

    output.write("\nBending Moment Distribution:\n")
    moment_dist_formatted = format_array(M_dist)
    for val in moment_dist_formatted:
        output.write(val + "\n")

    
'----------------------------- Numerical Extension of the classical LTT -------------------------'

# This is a numerical extension of the classical LTT as described by Anderson (2017)
# It is able to account for the non-linearities in the lift curve occurring when approach and 
# exceeding the stall angle of attack of the wing

dy = np.zeros(len(y)-1)             

for i in range(len(dy)):
    dy[i] = y[i+1] - y[i]           # spacing between spanwise stations [m]

gamma_guess = np.array(gamma)       # input circulation guess [m^2/s]
c_n = c                             # local chord length [m]
S = Aw                              # wing surface area [m^2]
tolerance = 0.0001                  # tolerance for convergence check
D = 0.05                            # damping factor

def calculate_dGamma_dy(gamma_input, dy):  # circulation
    dGamma_dy = np.zeros_like(gamma_input)
    for i in range(len(gamma_input) - 1):
        dGamma_dy[i] = (gamma_input[i + 1] - gamma_input[i]) / dy[i]
    return dGamma_dy

def calculate_induced_angle_of_attack(y, dGamma_dy, v_inf):    # induced AoA
    k = len(y) 
    alpha_i_new = np.zeros(k)
    
    for i in range(len(dy)):
        delta_y = dy[i]

    for n in range(k):
        sum_integral = 0
        for j in range(1, k-1, 2):  
            terms = []
            for offset in [-1, 0, 1]: 
                idx = j + offset
                if idx < 0 or idx >= k:
                    continue  

                if y[n] == y[idx]:
                    
                    adjacent_terms = []
                    for adj_offset in [-1, 1]:
                        adj_idx = idx + adj_offset
                        if 0 <= adj_idx < k and adj_idx != n:
                            adjacent_terms.append(dGamma_dy[adj_idx] / (y[n] - y[adj_idx]))
                    if adjacent_terms:
                        term = np.mean(adjacent_terms)
                    else:
                        continue
                else:
                    term = dGamma_dy[idx] / (y[n] - y[idx])
                terms.append(term)

            weight = 4 if offset == 0 else 1
            sum_integral += weight * sum(terms)
        
        alpha_i_new[n] = 1 / (4 * np.pi * v_inf) * (delta_y / 3) * sum_integral

    return alpha_i_new

def effective_AoA(alpha, alpha_i_new):  # effective angle of attack

    alpha_eff_new = alpha - alpha_i_new

    return alpha_eff_new

def get_cl_n(alpha_eff_new, alpha_foil, CL_foil): # new airfoil lift coefficient
    
    cl_n = np.interp(alpha_eff_new, alpha_foil, CL_foil)

    return cl_n

def calculate_new_circulation(v_inf, c_n, cl_n): # new circulation distribution

    gamma_new = 0.5 * v_inf * c_n * cl_n

    return gamma_new

def update_circulation(gamma_old, gamma_new, D):  # updated circulation distribution

    gamma_input = gamma_old + D * (gamma_new - gamma_old)

    return gamma_input

def check_convergence(gamma_old, gamma_new, tolerance): # convergence check

    diff = np.abs(gamma_new - gamma_old)
    return np.all(diff <= tolerance)

max_iterations = 300                # maximum number of allowed iterations
converged = False                   
iterations = 0                      # number of iterations
 
gamma_old = gamma_guess.copy()      # initial circulation guess
gamma_input = gamma_guess.copy()    # input circulation

while iterations < max_iterations and not converged:  # convergence criterions
    
    dGamma_dy_new = calculate_dGamma_dy(gamma_input, dy) 
    alpha_i_new = calculate_induced_angle_of_attack(y, dGamma_dy_new, v_inf)
    alpha_eff_new = alpha - alpha_i_new
    cl_n = get_cl_n(alpha_eff, alpha_foil, CL_foil)
    gamma_new = calculate_new_circulation(v_inf, c_n, cl_n)
    gamma_input = update_circulation(gamma_old, gamma_new, D)

    if check_convergence(gamma_old, gamma_new, tolerance):
        converged = True
    else: 
        gamma_old = gamma_input.copy()

    iterations += 1

def trapezoidal_integration(y, values): # trapezoidal rule

    integral = 0.0
    for i in range(1, len(y)):
        delta_y = abs(y[i] - y[i-1])  
        integral += (values[i] + values[i-1]) * delta_y / 2
    return integral

def calculate_lift_coefficient_trapezoidal(gamma, v_inf, S, y): # new wing lift coefficient

    integral = trapezoidal_integration(y, gamma)
    C_L = 2 / (v_inf * S) * 2 * integral
    return C_L

def calculate_drag_coefficient_trapezoidal(gamma_new, alpha_i_new, v_inf, S, y): # new wing induced drag coefficient
    
    product = gamma_new * abs(alpha_i_new)
    integral = trapezoidal_integration(y, product)
    C_Di = 2 / (v_inf * S) * 2 * integral
    return C_Di

#print(f'Converged after {iterations} iterations')
#print("New change in circulation :", dGamma_dy_new)
#print("New Induced Angle of Attack :", alpha_i_new)
#print("New effective Angle of Attack :", alpha_eff_new)
#print("Input circulation: ", gamma_input)
#print("Old circulation = ", gamma_old)
#print("New circulation = ", gamma_new)

CL_new = calculate_lift_coefficient_trapezoidal(gamma_new, v_inf, S, y)                     # wing lift coefficient with numerical LTT
CD_i_new = calculate_drag_coefficient_trapezoidal(gamma_new, alpha_i_new, v_inf, S, y)      # wing induced drag coefficient with numerical LTT
print("New lift coefficient: ", CL_new)
print("New induced drag coefficient :", CD_i_new)
CD_new = cd0 + CD_i_new                                                                     # wing total drag coefficient (neglecting wave drag and viscous drag)
#print("New drag coefficient :", CD_new)

