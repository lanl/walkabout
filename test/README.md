# Walkabout Verification and Acceptance Test Suite

General test directory for all versions of walkabout. Compares output to Scott Painter's original walkabout results (See Verification and Acceptance Tests in the user manual.) Test suite does not test amanzi capabilities, runtime changes, nor other critical failure conditions. OpenMP is currently set to 1 thread in test mode.

Instructions:

Run the python script proceeding with the path to the desire walkabout executable.

python walkabout_test.py /home/user/bin/walkabout
