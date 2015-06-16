SchemaTranslator
================

This is a tool for updating scenario (XML) files to use a later schema version.


Building
--------

As of July 2013, the schema translator has been partially converted to Kotlin.

To use, run the jar or see build instructions:

java -jar bin/schema-translator.jar

To build, first install Kotlin (via http://kotlin.jetbrains.org/ or an IDE
plugin or from the source). Then use the compile script:

./compile /PATH/TO/DIRECTORY/kotlinc

In case you need to do the compilation on a different platform and don't
understand my shell script, it's only four lines at the end which are important
(invoke kotlinc-jvm, javac, kotlinc-jvm again, then call jar).

For development the options are probably to use IntelliJ or a simple text
editor (there may be some Eclipse support).


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

