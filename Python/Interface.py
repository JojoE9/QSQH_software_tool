# Interface of QSQH tool

import os

# Get the current working directory
current_directory = os.getcwd() 

# Specify the subfolder name and the filename
InputFile_name = 'InputParameters.txt'     # Replace with your text file name

# Construct the absolute path to the text file
InputFile_path = os.path.join(current_directory, file_name)

# Check if the file exists before attempting to open it
if os.path.exists(file_path):
    try:
        # Open and read the file
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Do something with the file content
        print(file_content)

    except FileNotFoundError:
        print(f"The file '{file_path}' was not found.")
    except IOError as e:
        print(f"An error occurred while reading the file: {e}")
else:
    print(f"The file '{file_path}' does not exist.")

# Execute the code
exec(open('Execute_code.py').read())