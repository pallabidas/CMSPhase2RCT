## Set the top level module
#set_top UCTSummaryCard
set_top algo_unpacked
##
## Add source code
add_files src/algo_unpacked.cpp
#add_files src/UCTSummaryCard.cc
add_files src/adder_tree.cpp
#add_files src/am_sort_256x8.cpp
add_files src/egamma.cpp
add_files src/tau.cpp
add_files src/jet.cpp
add_files src/et_3by3.cpp
add_files src/PU_LUT.cpp
add_files src/boostedjet.cpp
add_files src/superindices.cpp 
##
## Add testbench files
#add_files -tb src/data.txt
#add_files -tb src/UCTSummaryCard_tb.cc
add_files -tb src/algo_unpacked_tb.cpp

add_files -tb data/test_random_set2_inp.txt
add_files -tb data/test_random_set2_out_ref.txt

add_files -tb data/tv_RCT1_3CS_set1_inp.txt
add_files -tb data/tv_RCT1_3CS_set1_out_ref.txt

add_files -tb data/tv_RCT2_3CS_set1_inp.txt
add_files -tb data/tv_RCT2_3CS_set1_out_ref.txt
