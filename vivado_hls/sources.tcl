## Set the top level module
set_top UCTSummaryCard
##
## Add source code
add_files src/UCTSummaryCard.cc
add_files src/adder_tree.cpp
add_files src/am_sort_256x8.cpp
add_files src/egamma.cpp
add_files src/tau.cpp
add_files src/jet.cpp
add_files src/et_3by3.cpp
add_files src/PU_LUT.cpp
##
## Add testbench files
add_files -tb src/data.txt
add_files -tb src/UCTSummaryCard_tb.cc

