# Author: Matthew Hancock
# Date: 10/18/21
# Description: This bash script times a single evaluation of the Xtal restraint under differential conditions.

# Test debug vs release.
source activate testenv
conda activate testenv

# Activate debug environment.
source /home/matthew/imp_development/debug/setup_environment.sh
python /home/matthew/xtal_benchmark/evaluation_time/src/single_evaluation.py

# Activate release environment (no debug compiler flags).
source /home/matthew/imp_development/release/setup_environment.sh
python /home/matthew/xtal_benchmark/evaluation_time/src/single_evaluation.py
