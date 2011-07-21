SchemaTranslator
================

This is a tool for updating scenario (XML) files to use a later schema version.


Building
--------

1.  Build OpenMalaria with CMake, then look in build/util/SchemaTranslator. Run
    from that folder or copy jar here.
2.  This folder has eclipse project files. In theory, you should be able to
    create a workspace in util, then add this directory as an existing project.


Usage
-----

Copy any XML files to update into the scenarios/ folder, then run the tool:

1.  From the jar: make sure the working directory is the one containing the
    "scenarios" folder of input files, then run
    
        java -jar SchemaTranslator.jar
2.  From the class files: as above, but use the command:
    
        java SchemaTranslator
3.  From eclipse: run as java application (you may need to make sure the
    working directory is correct).

