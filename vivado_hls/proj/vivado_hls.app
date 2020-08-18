<AutoPilot:project xmlns:AutoPilot="com.autoesl.autopilot.project" projectType="C/C++" name="proj" top="algo_unpacked">
    <Simulation argv="tv_RCT1_3CS_set1">
        <SimFlow name="csim" setup="false" optimizeCompile="false" clean="false" ldflags="" mflags=""/>
    </Simulation>
    <includePaths/>
    <libraryFlag/>
    <files>
        <file name="../../data/tv_RCT2_3CS_set1_out_ref.txt" sc="0" tb="1" cflags="  -Wno-unknown-pragmas"/>
        <file name="../../data/tv_RCT2_3CS_set1_inp.txt" sc="0" tb="1" cflags="  -Wno-unknown-pragmas"/>
        <file name="../../data/tv_RCT1_3CS_set1_out_ref.txt" sc="0" tb="1" cflags="  -Wno-unknown-pragmas"/>
        <file name="../../data/tv_RCT1_3CS_set1_inp.txt" sc="0" tb="1" cflags="  -Wno-unknown-pragmas"/>
        <file name="../../data/test_random_set2_out_ref.txt" sc="0" tb="1" cflags="  -Wno-unknown-pragmas"/>
        <file name="../../data/test_random_set2_inp.txt" sc="0" tb="1" cflags="  -Wno-unknown-pragmas"/>
        <file name="../../src/algo_unpacked_tb.cpp" sc="0" tb="1" cflags="  -Wno-unknown-pragmas"/>
        <file name="src/boostedjet.cpp" sc="0" tb="false" cflags=""/>
        <file name="src/PU_LUT.cpp" sc="0" tb="false" cflags=""/>
        <file name="src/et_3by3.cpp" sc="0" tb="false" cflags=""/>
        <file name="src/jet.cpp" sc="0" tb="false" cflags=""/>
        <file name="src/tau.cpp" sc="0" tb="false" cflags=""/>
        <file name="src/egamma.cpp" sc="0" tb="false" cflags=""/>
        <file name="src/am_sort_256x8.cpp" sc="0" tb="false" cflags=""/>
        <file name="src/adder_tree.cpp" sc="0" tb="false" cflags=""/>
        <file name="src/UCTSummaryCard.cc" sc="0" tb="false" cflags=""/>
        <file name="src/algo_unpacked.cpp" sc="0" tb="false" cflags=""/>
    </files>
    <solutions>
        <solution name="solution1" status=""/>
    </solutions>
</AutoPilot:project>

