# Interface of QSQH tool

import os
import subprocess

# Get the current working directory
current_directory = os.getcwd() 

# Specify the the input filename
InputFile_name = 'InputParameters.txt' 

# Construct the absolute path to the InputParameters.txt
InputFile_path = os.path.join(current_directory, InputFile_name)

# Check if the file exists before attempting to open it
if os.path.exists(InputFile_path):
    try:
        # Open and read the file
        with open(InputFile_path, 'r') as file:
            lines = file.readlines()

        # Do something with the file content
        print(lines)

    except FileNotFoundError:
        print(f"The file '{InputFile_path}' was not found.")
    except IOError as e:
        print(f"An error occurred while reading the file: {e}")
else:
    print(f"The file '{InputFile_path}' does not exist.")

# Specify the the execution filename
Execution_name = 'Execution_code.py' 

# Construct the absolute path to the Execution_code.py
Execution_path = os.path.join(current_directory,Execution_name)

try:
    # Run the specified Python script
    subprocess.run(["python", Execution_path])
except Exception as e:
    print(f"An error occurred: {e}")