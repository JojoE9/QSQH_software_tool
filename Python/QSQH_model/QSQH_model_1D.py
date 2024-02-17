import os
import numpy as np
import h5py
from scipy.interpolate import CubicSpline
from tqdm import tqdm

# Size of the input friction field in wall unit
Lx = float(InputParameters[0])  # Length of the stream-wise domain
Lz = float(InputParameters[1])  # Length of the span-wise domain
ymin = float(InputParameters[4])  # y = {0, 100}
ymax = float(InputParameters[5])

# Is this run a test? (true/false)
test = bool(InputParameters[6])

# Starting loop (starting from 1 for a new run)
starting_loop = int(InputParameters[7])

# Ending loop
ending_loop = int(InputParameters[8])

# Is this a test?
if test == True:
    folderending = '_test'
    np.random.seed(42)
else:
    folderending = ''

# The directory of the input friction field
dir_input = os.path.join(current_directory, 'input_friction_field'+folderending)

# The directory of the reference wall-normal grid
dir_y_ref = os.path.join(current_directory, 'QSQH_model')

# The directory of the 1D universal velocity field
dir_vel_tilde = os.path.join(current_directory, 'universal_velocity_field_1D'+folderending)

# The folder that will store the 1D output of the synthetic model
dir_output = os.path.join(current_directory, 'output_1D'+folderending)

# Data converter
lambda_x_plus = 2000 * np.pi  # Stream-wise cut-off
lambda_z_plus = 250 * np.pi   # Span-wise cut-off

# Corresponding cut-off wavenumbers
Nx_cutoff = int(np.round(Lx / lambda_x_plus))
Nz_cutoff = int(np.round(Lz / lambda_z_plus))

# Load the list of names of the input data
Files_input = [file for file in os.listdir(dir_input) if file.startswith('tau_')]

# Apply the large-scale filter to the input data
for filename in Files_input:

    # Open the HDF5 file in read mode
    with h5py.File(os.path.join(dir_input,filename), 'r') as file:
        # List all the groups in the HDF5 file

        # Access a specific dataset within a group
        dataset_Lx = file['/Lx']
        dataset_Lz = file['/Lz']
        dataset_taux = file['/taux']
        dataset_tauz = file['/tauz']
    
        # Read the data from the dataset
        Lx = dataset_Lx[()] 
        Lz = dataset_Lz[()] 
        taux = dataset_taux[()] 
        tauz = dataset_tauz[()]  # Read all data into memory
        # Alternatively, you can specify a slice or indices to read a subset of the data
    
    taux_hat = np.fft.fft2(taux)  # 2D Fast Fourier Transform, stream-wise velocity
    tauz_hat = np.fft.fft2(tauz)  # 2D Fast Fourier Transform, span-wise velocity

    nz, nx = taux.shape

    # Zero out the specified modes
    taux_hat[(Nz_cutoff+2):(nz-Nz_cutoff), :] = 0
    taux_hat[:, (Nx_cutoff+2):(nx-Nx_cutoff)] = 0
    tauz_hat[(Nz_cutoff+2):(nz-Nz_cutoff), :] = 0
    tauz_hat[:, (Nx_cutoff+2):(nx-Nx_cutoff)] = 0

    tauxL = np.fft.ifft2(taux_hat).real  # 2D Inverse FFT
    tauzL = np.fft.ifft2(tauz_hat).real

    tauL = np.sqrt(tauxL ** 2 + tauzL ** 2)  # Resultant wall shear stress
    u_tauL_plus = np.sqrt(tauL)  # Large-scale friction velocity

    # Direction of the large-scale friction velocity
    theta = -np.arctan2(tauzL, tauxL)

    # Name of files of the filtered field
    newname = filename.replace('tau', 'u_tauL')

    # Open the HDF5 file in write mode
    with h5py.File(os.path.join(dir_input,newname), "w") as file:
        # Create a dataset within the HDF5 file and save the array
        file.create_dataset('u_tauL_plus', data=u_tauL_plus)
        file.create_dataset('theta', data=theta)

# Load the list of names of the output data
Files_output = [file for file in os.listdir(dir_input) if file.startswith('u_tauL_')]

# Mean value of u_tauL_plus and theta
u_tauL_mean_loc = []
theta_mean_loc = []

for filename in Files_output:

    # Open the HDF5 file in read mode
    with h5py.File(os.path.join(dir_input,filename), 'r') as file:
        # List all the groups in the HDF5 file

        # Access a specific dataset within a group
        dataset_u_tauL_plus = file['/u_tauL_plus']
        dataset_theta = file['/theta']
    
        # Read the data from the dataset
        u_tauL_plus = dataset_u_tauL_plus[()] 
        theta = dataset_theta[()]  # Read all data into memory
        # Alternatively, you can specify a slice or indices to read a subset of the data

    u_tauL_mean_loc.append(np.mean(u_tauL_plus))
    theta_mean_loc.append(np.mean(theta))

u_tauL_plus_mean = np.mean(u_tauL_mean_loc)
theta_mean = np.mean(theta_mean_loc)

# Second moments of u_tauL_plus and theta
u_tauL_var_loc = []
theta_var_loc = []
theta_var_include_theta_mean_loc = []

for filename in Files_output:
    
    # Open the HDF5 file in read mode
    with h5py.File(os.path.join(dir_input,filename), 'r') as file:
        # List all the groups in the HDF5 file

        # Access a specific dataset within a group
        dataset_u_tauL_plus = file['/u_tauL_plus']
        dataset_theta = file['/theta']
    
        # Read the data from the dataset
        u_tauL_plus = dataset_u_tauL_plus[()] 
        theta = dataset_theta[()]  # Read all data into memory
        # Alternatively, you can specify a slice or indices to read a subset of the data

    u_tauL_var_loc.append(np.mean(((u_tauL_plus / u_tauL_plus_mean - 1) ** 2)))
    theta_var_loc.append(np.mean((theta ** 2)))
    theta_var_include_theta_mean_loc.append(np.mean(((theta - theta_mean) ** 2)))

# Variances of u_tauL and theta
u_tauL_var = np.mean(u_tauL_var_loc)
theta_var = np.mean(theta_var_loc)
theta_var_include_theta_mean = np.mean(theta_var_include_theta_mean_loc)

# Store the statistics of the output data
with h5py.File(os.path.join(dir_input,'utauL_stats.h5'),'w') as file:
    file.create_dataset('u_tauL_plus_mean', data=u_tauL_plus_mean)
    file.create_dataset('theta_mean', data=theta_mean)
    file.create_dataset('u_tauL_var', data=u_tauL_var)
    file.create_dataset('theta_var', data=theta_var)
    file.create_dataset('theta_var_include_theta_mean', data=theta_var_include_theta_mean)

# Read the universal velocity field
with h5py.File(os.path.join(dir_vel_tilde,'vel_tilde_1D.h5'), 'r') as file:
    # Access a specific dataset within a group
    dataset_u_tilde = file['/u_tilde']
    dataset_v_tilde = file['/v_tilde']
    dataset_w_tilde = file['/w_tilde']
    dataset_y_tilde = file['/y_tilde']
    
    # Read the data from the dataset
    u_tilde = dataset_u_tilde[()] 
    v_tilde = dataset_v_tilde[()]  
    w_tilde = dataset_w_tilde[()]
    y_tilde = dataset_y_tilde[()]   # Read all data into memory
    # Alternatively, you can specify a slice or indices to read a subset of the data

# Read the reference y grid
with h5py.File(os.path.join(dir_y_ref,'y_ref.h5'), 'r') as file:
    # Access a specific dataset within a group
    dataset_y = file['/y']
    
    # Read the data from the dataset
    y = dataset_y[()]   # Read all data into memory
    # Alternatively, you can specify a slice or indices to read a subset of the data

# Read the utauL statistics
with h5py.File(os.path.join(dir_input,'utauL_stats.h5'), 'r') as file:
    # Access a specific dataset within a group
    dataset_u_tauL_plus_mean = file['/u_tauL_plus_mean']
    dataset_theta_mean = file['/theta_mean']
    
    # Read the data from the dataset
    u_tauL_mean = dataset_u_tauL_plus_mean[()]
    theta_mean = dataset_theta_mean[()]   # Read all data into memory
    # Alternatively, you can specify a slice or indices to read a subset of the data


# Set up the rid
idmin = np.argmin(np.abs(ymin - y))
idmax = np.argmin(np.abs(ymax - y))
y_plus_ref = y[:,idmin:idmax]
y_plus = y_plus_ref

# Store the y of the output data
with h5py.File(os.path.join(dir_output,'y_plus.h5'),'w') as file:
    file.create_dataset('y_plus', data=y_plus)

# QSQH model
Files = [file for file in os.listdir(dir_input) if file.startswith('u_tauL_')]

# Read the firs utauL
with h5py.File(os.path.join(dir_input,Files[0]), 'r') as file:
    # Access a specific dataset within a group
    dataset_u_tauL_plus = file['/u_tauL_plus']

    # Read the data from the dataset
    u_tauL_plus = dataset_u_tauL_plus[()] # Read all data into memory
    # Alternatively, you can specify a slice or indices to read a subset of the data

X = u_tauL_plus.shape[0]
Z = u_tauL_plus.shape[1]

pw_x = np.ceil(np.log2(X))
pw_z = np.ceil(np.log2(Z))
pw_t = np.ceil(np.log2(len(Files)))

Nx0 = 2
Nz0 = 2
T0 = 2
pw = starting_loop - 1

while pw < pw_x or pw < pw_z or pw < pw_t:
    
    pw = pw+1

    if pw > ending_loop:
        break

    if pw >= pw_x:
        Nx = X
    else:
        Nx = Nx0 ** pw

    if pw >= pw_z:
        Nz = Z
    else:
        Nz = Nz0 ** pw

    if pw >= pw_t:
        T = len(Files)
    else:
        T = T0 ** pw
    
    nx = np.round(np.linspace(0, X - 1, Nx)).astype(int)
    nz = np.round(np.linspace(0, Z - 1, Nz)).astype(int)
    tstep = np.round(np.linspace(0, len(Files)-1, T)).astype(int)

    u_plus_qsqh = np.zeros((y_plus_ref.size, Nx * Nz, T))
    v_plus_qsqh = np.zeros((y_plus_ref.size, Nx * Nz, T))
    w_plus_qsqh = np.zeros((y_plus_ref.size, Nx * Nz, T))

    for n in tqdm(tstep, desc=f"Loop {pw:03d}", leave=False):

        # Read utauL
        with h5py.File(os.path.join(dir_input,Files[n]), 'r') as file:
            dataset_u_tauL_plus = file['/u_tauL_plus']
            dataset_theta = file['/theta']

            u_tauL_plus = dataset_u_tauL_plus[()] 
            theta = dataset_theta[()]

        u_tauL = u_tauL_plus[nx[:, None], nz]
        theta = theta[nx[:, None], nz]

        u_tauL = u_tauL.reshape(-1,1)
        theta = theta.reshape(-1,1)

        if test == True:
            np.random.seed(42)
        else:
            np.random.seed()

        idx = np.random.randint(u_tilde.shape[1]+1, size = (len(u_tauL), 1))

        u_tilde_sample = u_tilde[:,idx.flatten()]
        v_tilde_sample = v_tilde[:,idx.flatten()]
        w_tilde_sample = w_tilde[:,idx.flatten()]
        y_tilde_sample = y_tilde[:,idx.flatten()]

        for num in range(u_tauL.size):
            u_plus = (u_tilde_sample[:,num] * np.cos(theta[num]) +  w_tilde_sample[:,num] * np.sin(theta[num])) * u_tauL[num] * u_tauL_plus_mean

            v_plus = v_tilde_sample[:,num] * u_tauL[num] * u_tauL_plus_mean

            w_plus = (-u_tilde_sample[:,num] * np.sin(theta[num]) + w_tilde_sample[:,num] * np.cos(theta[num])) * u_tauL[num] * u_tauL_plus_mean

            y_plus = (y_tilde_sample[:,num] / u_tauL[num]) / u_tauL_plus_mean

            # Interpolate using CubicSpline
            cs_u = CubicSpline(y_plus, u_plus)
            cs_v = CubicSpline(y_plus, v_plus)
            cs_w = CubicSpline(y_plus, w_plus)

            u_plus_interp = cs_u(y_plus_ref)
            v_plus_interp = cs_v(y_plus_ref)
            w_plus_interp = cs_w(y_plus_ref)

            u_plus_qsqh[:, n, num] = u_plus_interp
            v_plus_qsqh[:, n, num] = v_plus_interp
            w_plus_qsqh[:, n, num] = w_plus_interp

    # Save data
    
    with h5py.File(os.path.join(dir_output, f'vel_plus_qsqh_loop_{pw:03d}{folderending}'),'w') as file:
        file.create_dataset('u_plus_qsqh', data=u_plus_qsqh)
        file.create_dataset('v_plus_qsqh', data=v_plus_qsqh)
        file.create_dataset('w_plus_qsqh', data=w_plus_qsqh)
