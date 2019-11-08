# useful macros/fcts
# imported from waves project
# Authors Romain Boman


MACRO(MACRO_AddTest srcDir)
    #message(STATUS "Adding test directory ${srcDir}")
    if(PYTHON_WRAPPER)
        set(inputFile ".*_fluid.py")
        set(runCommand ${PYTHON_EXECUTABLE})
        set(runPy "run.py")
    else()
        set(inputFile "infile.*arp")
        set(runCommand "./bin/VLM")
        set(runPy "")
    endif()
    
    file(GLOB tfiles RELATIVE ${srcDir} ${srcDir}/*)
    #message(STATUS "tfiles=${tfiles}")
    foreach(tfile ${tfiles})
        set(spath ${srcDir}/${tfile})
        if(NOT IS_DIRECTORY ${spath} AND ${spath} MATCHES ${inputFile})
            string(REPLACE "${PROJECT_SOURCE_DIR}/" "" strip ${spath}) 
			file(TO_NATIVE_PATH "${strip}" strip)
            message(STATUS "Adding test ${strip}")
            add_test(NAME ${strip} 
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} 
                     COMMAND ${runCommand} ${runPy} ${strip})
        else()
            MACRO_AddTest(${srcDir}/${tfile})
        endif()
    endforeach()
ENDMACRO()
