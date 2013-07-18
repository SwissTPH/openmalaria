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
