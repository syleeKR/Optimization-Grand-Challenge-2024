# Optimization-Grand-Challenge-2024

**Winner, Team DMS**  
With Teammate [@prisem123](https://github.com/prisem123)

For more information, visit the [Official Site](https://optichallenge.com)


## Getting started
On your linux platform,
```bash
$ git clone https://github.com/syleeKR/Optimization-Grand-Challenge-2024
$ conda env create -f ogc2024_env.yml    
$ conda activate ogc2024    
$ pip install pybind11  
```
*Gurobi License is required.*


## Files

- **src/**: Contains source code files
  - **myalgorithm.py**: Main interface file
  - **setup.py** : for pybind
  - **build.sh** : compile 
  - **GurobiSP.py** : Gurobi for solving set-partion problem
  - **c_backend/**: c++ files implementing the actual algorithm
    - **runner.hpp**: overall ALNS 
    - **main_iteration.hpp**: main iteration of ALNS
    - **repair.hpp**: Insertion of ALNS
    - **repair_tools/**: insertion method implementation 
    - **destroy_solution.hpp**: Destruction of ALNS
    - **destroy_tools/**: removal methods implementation 
    - **annealer.hpp**: simulated annealing
    - **hyperparam.hpp**: hyperparameters

