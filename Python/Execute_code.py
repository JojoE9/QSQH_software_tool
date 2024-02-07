# Extract values after comma ',' from the 'lines' list
InputParameters = [line.split(',')[1] for line in lines]

# Convert strings to floating-point numbers
xmax = float(InputParameters[2])
zmax = float(InputParameters[3])

if xmax == 0 and zmax == 0:
    script_path = os.path.join(current_directory, 'QSQH_model', 'QSQH_model_1D')
else:
    script_path = os.path.join(current_directory, 'QSQH_model','QSQH_model_3D')

# Execute the corresponding script
try:
    # Run the specified Python script
    subprocess.run(["python", script_path])
except Exception as e:
    print(f"An error occurred: {e}")
