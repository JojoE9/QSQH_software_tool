import os

# Extract values after comma ',' from the 'lines' list
InputParameters = lines.split(',')[1:]

# Convert strings to floating-point numbers
xmax = float(InputParameters[2])
zmax = float(InputParameters[3])

# Get the current working directory and file separator
workdir = os.getcwd()
file_separator = InputParameters[9]

if xmax == 0 and zmax == 0:
    script_path = os.path.join(workdir, 'QSQH_model', 'QSQH_model_1D')
else:
    script_path = os.path.join(workdir, 'QSQH_model', 'QSQH_model_3D')

# Execute the corresponding script
os.system(script_path)
