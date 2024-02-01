import os

# Interface of QSQH tool

# Clear the screen
os.system('clear')  # On Windows, you can use 'cls' instead of 'clear'

# User input parameters
with open('InputParameters.txt', 'r') as file:
    lines = file.readlines()

# Execute the code
exec(open('Excute_code.py').read())
