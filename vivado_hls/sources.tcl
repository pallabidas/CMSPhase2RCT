## Set the top level module
set_top algo_unpacked
##
## Add source code
add_files src/algo_unpacked.cpp
add_files src/adder_tree.cpp
add_files src/et_3by3.cpp
add_files src/PU_LUT.cpp
add_files src/boostedjet.cpp
add_files src/bitonicFun64.cpp
add_files src/bitonicSort64.cpp
add_files src/bitonic32Inc.cpp
add_files src/bitonic32Dec.cpp
##
## Add testbench files
add_files -tb src/algo_unpacked_tb.cpp

add_files -tb data/test_random112_inp.txt
add_files -tb data/test_random112_out_ref.txt
