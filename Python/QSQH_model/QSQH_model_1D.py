import os
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.io import loadmat, savemat
from progressbar import ProgressBar

# Size of the input friction field in wall unit
Lx = float(InputParameters[0])  # Length of the stream-wise domain
Lz = float(InputParameters[1])  # Length of the span-wise domain

# Choose the maximum spatial ranges in wall unit
# (xmax=0 & zmax=0 for 1D output field)
# xmax = float(InputParameters[2])  # xmax = {0, 204}
# zmax = float(InputParameters[3])  # zmax = {0, 120}
ymin = float(InputParameters[4])  # y = {0, 100}
ymax = float(InputParameters[5])

# Is this run a test? (true/false)
test = bool(int(InputParameters[6]))

# Starting loop (starting from 1 for a new run)
starting_loop = int(InputParameters[7])

# Ending loop
ending_loop = int(InputParameters[8])

# Is this a test?
if test:
    folderending = '_test'
    s = loadmat(os.path.join(workdir, 'QSQH_model', 'RNGseed.mat'))['s'][0, 0]
else:
    folderending = ''
    s = np.random.randint(2**32)

# The directory of the input friction field
dir_input = os.path.join(workdir, 'input_friction_field' + folderending)

# The directory of the reference wall-normal grid
dir_y_ref = os.path.join(workdir, 'QSQH_model')

# The directory of the 1D universal velocity field
dir_vel_tilde = os.path.join(workdir, 'universal_velocity_field_1D' + folderending)

# The folder that will store the 1D output of the synthetic model
dir_output = os.path.join(workdir, 'output_1D' + folderending)

# Data converter
lambda_x_plus = 2000 * np.pi  # Stream-wise cut-off
lambda_z_plus = 250 * np.pi   # Span-wise cut-off

# Corresponding cut-off wavenumbers
kx = int(Lx / lambda_x_plus)
kz = int(Lz / lambda_z_plus)

# Load the list of names of the input data
Files_input = [file for file in os.listdir(dir_input) if file.startswith('tau_')]

# Apply the large-scale filter to the input data
for filename in Files_input:
    LS_wall = loadmat(os.path.join(dir_input, filename))
    taux = LS_wall['taux']
    tauz = LS_wall['tauz']

    nx, nz = taux.shape
    taux_hat = np.fft.fft2(taux)  # 2D Fast Fourier Transform, stream-wise velocity
    tauz_hat = np.fft.fft2(tauz)  # 2D Fast Fourier Transform, span-wise velocity

    taux_hat[(kx + 2):(nx - kx), :] = 0  # Zero out the specified modes
    taux_hat[:, (kz + 2):(nz - kz)] = 0
    tauz_hat[(kx + 2):(nx - kx), :] = 0
    tauz_hat[:, (kz + 2):(nz - kz)] = 0

    tauxL = np.fft.ifft2(taux_hat).real  # 2D Inverse FFT
    tauzL = np.fft.ifft2(tauz_hat).real

    tauL = np.sqrt(tauxL ** 2 + tauzL ** 2)  # Resultant wall shear stress
    u_tauL_plus = np.sqrt(tauL)  # Large-scale friction velocity

    # Direction of the large-scale friction velocity
    theta = -np.arctan2(tauzL, tauxL)

    # Name of files of the filtered field
    newname = filename.replace('tau', 'u_tauL')

    # Store the output data
    savemat(os.path.join(dir_input, newname), {'u_tauL_plus': u_tauL_plus, 'theta': theta}, format='7.3')

# Load the list of names of the output data
Files_output = [file for file in os.listdir(dir_input) if file.startswith('u_tauL_')]

# Mean value of u_tauL_plus and theta
u_tauL_mean_loc = []
theta_mean_loc = []

for filename in Files_output:
    LS_wall = loadmat(os.path.join(dir_input, filename))
    u_tauL_plus = LS_wall['u_tauL_plus']
    theta = LS_wall['theta']

    u_tauL_mean_loc.append(np.mean(u_tauL_plus))
    theta_mean_loc.append(np.mean(theta))

u_tauL_plus_mean = np.mean(u_tauL_mean_loc)
theta_mean = np.mean(theta_mean_loc)

# Second moments of u_tauL_plus and theta
u_tauL_var_loc = []
theta_var_loc = []
theta_var_include_theta_mean_loc = []

for filename in Files_output:
    LS_wall = loadmat(os.path.join(dir_input, filename))
    u_tauL_plus = LS_wall['u_tauL_plus']
    theta = LS_wall['theta']

    u_tauL_var_loc.append(np.mean(((u_tauL_plus / u_tauL_plus_mean - 1) ** 2)))
    theta_var_loc.append(np.mean((theta ** 2)))
    theta_var_include_theta_mean_loc.append(np.mean(((theta - theta_mean) ** 2)))

# Variances of u_tauL and theta
u_tauL_var = np.mean(u_tauL_var_loc)
theta_var = np.mean(theta_var_loc)
theta_var_include_theta_mean = np.mean(theta_var_include_theta_mean_loc)

# Store the statistics of the output data
savemat(os.path.join(dir_input, 'utauL_stats.mat'), {'u_tauL_plus_mean': u_tauL_plus_mean,
                                                    'theta_mean': theta_mean,
                                                    'u_tauL_var': u_tauL_var,
                                                    'theta_var': theta_var,
                                                    'theta_var_include_theta_mean': theta_var_include_theta_mean},
        format='7.3')

# Load the input data
vel_tilde = loadmat(os.path.join(dir_vel_tilde, 'vel_tilde_1D.mat'))['vel_tilde']
y_ref = loadmat(os.path.join(dir_y_ref, 'y_ref.mat'))['y_ref']
stats = loadmat(os.path.join(dir_input, 'utauL_stats.mat'))
u_tauL_plus_mean = stats['u_tauL_plus_mean'][0, 0]
theta_mean = stats['theta_mean'][0, 0]

# Set up the rid
ymin_id = np.argmin(np.abs(ymin - y_ref))
ymax_id = np.argmin(np.abs(ymax - y_ref))
y_plus_ref = y_ref[ymin_id:ymax_id + 1]

# QSQH model
Files = [file for file in os.listdir(dir_input) if file.startswith('u_tauL_')]
pw_x = int(np.ceil(np.log2(u_tauL_plus.shape[0])))
pw_z = int(np.ceil(np.log2
