cmake_minimum_required (VERSION 2.4)
project (OpenMalaria-doc)


# our markdown -> html convertor
find_program (MD_CONVERTOR NAMES pandoc)
find_program (SED NAMES sed)
set (SED_FILE "${CMAKE_SOURCE_DIR}/resource/finish.sed")

# dummy target
add_custom_target (MARKDOWN2HTML ALL VERBATIM)

# easy way to find files...
file (GLOB_RECURSE MD_FILES "${CMAKE_SOURCE_DIR}/*.markdown")

foreach (MD_FILE ${MD_FILES})
    # set HTML_FILE
    get_filename_component (FILE_NAME_WE ${MD_FILE} NAME_WE)
    get_filename_component (FILE_PATH ${MD_FILE} PATH)
    #message (STATUS "input:  ${MD_FILE}")
    #message (STATUS "output: ${FILE_PATH}/${FILE_NAME_WE}.html")
    set (HTML_FILE "${FILE_PATH}/${FILE_NAME_WE}.html")
    
    add_custom_command (
        TARGET MARKDOWN2HTML
        POST_BUILD
        COMMAND ${MD_CONVERTOR} ${MD_FILE} | ${SED} -f ${SED_FILE} > ${HTML_FILE}
        COMMENT "Converting ${MD_FILE}"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        VERBATIM
    )
endforeach (MD_FILE)
