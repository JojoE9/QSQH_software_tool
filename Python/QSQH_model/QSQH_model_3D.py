import os
import numpy as np
import h5py
from scipy.interpolate import griddata
import random

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
dir_vel_tilde = os.path.join(current_directory, 'universal_velocity_field_3D'+folderending)

# The folder that will store the 1D output of the synthetic model
dir_output = os.path.join(current_directory, 'output_3D'+folderending)

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

# Read the reference y grid
with h5py.File(os.path.join(dir_y_ref,'y_ref.h5'), 'r') as file:
    # Access a specific dataset within a group
    dataset_y = file['/y']
    
    # Read the data from the dataset
    y = dataset_y[()]   # Read all data into memory
    # Alternatively, you can specify a slice or indices to read a subset of the data


# Set up the uniform grid
x_plus_ref = np.arange(-xmax/2, xmax/2 + 12, 12)
x_plus_ref = x_plus_ref.reshape(-1,1)
z_plus_ref = np.arange(-zmax/2, zmax/2 + 6, 6)
z_plus_ref = z_plus_ref.reshape(-1,1)

idmin = np.argmin(np.abs(ymin - y))
idmax = np.argmin(np.abs(ymax - y))
y_plus_ref = y[:,idmin:idmax]
y_plus_ref = y_plus_ref.T
y_plus = y_plus_ref

Yref, Zref, Xref = np.meshgrid(y_plus_ref, z_plus_ref, x_plus_ref, indexing='ij')

x_plus = x_plus_ref
y_plus = y_plus_ref
z_plus = z_plus_ref

with h5py.File(os.path.join(dir_input,'boxsize.h5'),'w') as file:
    file.create_dataset('x_plus', data=x_plus)
    file.create_dataset('z_plus', data=z_plus)
    file.create_dataset('y_plus', data=y_plus)

# Read the firs utauL
Files = [file for file in os.listdir(dir_input) if file.startswith('u_tauL_')]

with h5py.File(os.path.join(dir_input,Files[0]), 'r') as file:
    # Access a specific dataset within a group
    dataset_u_tauL_plus = file['/u_tauL_plus']

    # Read the data from the dataset
    u_tauL_plus = dataset_u_tauL_plus[()] # Read all data into memory
    # Alternatively, you can specify a slice or indices to read a subset of the data

# Read the tilde box
with h5py.File(os.path.join(dir_vel_tilde,'xz_tilde.h5'), 'r') as file:
    # Access a specific dataset within a group
    dataset_x_tilde = file['/x_tilde']
    dataset_z_tilde = file['/z_tilde']

    # Read the data from the dataset
    x_tilde = dataset_x_tilde[()]
    z_tilde = dataset_z_tilde[()] # Read all data into memory
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

    u_plus_qsqh = np.zeros((T, Nx * Nz, y_plus_ref.size, z_plus_ref.size, x_plus_ref.size))
    v_plus_qsqh = np.zeros((T, Nx * Nz, y_plus_ref.size, z_plus_ref.size, x_plus_ref.size))
    w_plus_qsqh = np.zeros((T, Nx * Nz, y_plus_ref.size, z_plus_ref.size, x_plus_ref.size))

    all_uarrays = []
    all_varrays = []
    all_warrays = []

    for n in range(T):

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

        TildeFileList = [file for file in os.listdir(dir_vel_tilde) if file.startswith('vel_tilde_')]

        if test:
            random.seed(42)
        else:
            random.seed()

        index = random.choices(range(len(TildeFileList)), k=len(u_tauL))

        for k in range(len(u_tauL)):

            # Read the tilde field
            with h5py.File(os.path.join(dir_vel_tilde,TildeFileList[index[1]]), 'r') as file:
                # Access a specific dataset within a group
                dataset_u_tilde = file['/u_tilde']
                dataset_v_tilde = file['/v_tilde']
                dataset_w_tilde = file['/w_tilde']
                dataset_y_tilde = file['/y_tilde']

                # Read the data from the dataset
                u_tilde = dataset_u_tilde[()]
                v_tilde = dataset_v_tilde[()] 
                w_tilde = dataset_w_tilde[()]
                y_tilde = dataset_y_tilde[()] # Read all data into memory
                # Alternatively, you can specify a slice or indices to read a subset of the data

            u_plus = (u_tilde * np.cos(theta[k]) +  w_tilde * np.sin(theta[k])) * u_tauL[k] * u_tauL_plus_mean
            v_plus = v_tilde * u_tauL[k] * u_tauL_plus_mean
            w_plus = (-u_tilde * np.sin(theta[k]) + w_tilde * np.cos(theta[k])) * u_tauL[k] * u_tauL_plus_mean
            y_plus = (y_tilde / u_tauL[k]) / u_tauL_plus_mean

            x_plus_r = (x_tilde / u_tauL[k]) / u_tauL_plus_mean
            x_plus_r -= np.max(x_plus_r) / 2
            z_plus_r = (z_tilde / u_tauL[k]) / u_tauL_plus_mean
            z_plus_r -= np.max(z_plus_r) / 2

            Yq, Zq, Xq = np.meshgrid(y_plus, z_plus_r, x_plus_r, indexing='ij')

            Xrot = Xq * np.cos(theta[k]) + Zq * np.sin(theta[k])
            Zrot = -Xq * np.sin(theta[k]) + Zq * np.cos(theta[k])

            u_plus_qsqh_loc = griddata((Yq.flatten(), Zrot.flatten(), Xrot.flatten()), u_plus.flatten(), (Yref, Zref, Xref), method='linear')
            v_plus_qsqh_loc = griddata((Yq.flatten(), Zrot.flatten(), Xrot.flatten()), v_plus.flatten(), (Yref, Zref, Xref), method='linear')
            w_plus_qsqh_loc = griddata((Yq.flatten(), Zrot.flatten(), Xrot.flatten()), w_plus.flatten(), (Yref, Zref, Xref), method='linear')

            all_uarrays.append(u_plus_qsqh_loc)
            all_varrays.append(v_plus_qsqh_loc)
            all_warrays.append(w_plus_qsqh_loc)

    u_plus_qsqh = np.stack(all_uarrays, axis=0)
    v_plus_qsqh = np.stack(all_varrays, axis=0)
    w_plus_qsqh = np.stack(all_warrays, axis=0)

    with h5py.File(os.path.join(dir_output, f'vel_plus_qsqh_loop_{pw:03d}{folderending}'),'w') as file:
        file.create_dataset('u_plus_qsqh', data=u_plus_qsqh)
        file.create_dataset('v_plus_qsqh', data=v_plus_qsqh)
        file.create_dataset('w_plus_qsqh', data=w_plus_qsqh)


