# QSQH synthetic turbulence model

## 📁 Project Structure
```bash

.
├── docs/                      # User manual
├── model/                     # QSQH synthetic turbulence model
│   └── core/                  # Core code
│   └── Interface.m            # User interface
│   └── InputParameters.txt    # User input of the model
│   └── loop_calculator.m      # Calculator of the size of output
│   └── Stats_calculator.m     # Calculator of output statistics
└── validation/                # Benchmark datasets
```

## Description of Important Files/Folders


|        Files/Folders        |                     Utility                    |
|-----------------------------|------------------------------------------------|
|           'core/'           |             Core codes of QSQH STM             |
|         'Interface'         |           User Interface of the tool           |
|      'loop_calculator'      |       Determine the size of output field       |
|      'Stats_calculator'     |    Calculate the statistics of output field    |