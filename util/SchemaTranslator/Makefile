# Script to build the SchemaTranslator under linux.
# (If somebody wants to write configuration for a proper build tool or
# something compatible with other platforms, feel free, but those are beyond
# my requirements.)

# Default location to look for Kotlin
KOTLIN_DIR=/home/install/kotlinc
KOTLINC=$(KOTLIN_DIR)/bin/kotlinc-jvm
KOTLIN_RUNTIME=$(KOTLIN_DIR)/lib/kotlin-runtime.jar

# Assume javac and jar are available:
JAVAC=/usr/bin/javac
JAR=/usr/bin/jar

# Location for built classes
CP=classes
# Location of sources
SRC=src

# Location for output jar
ST_JAR=bin/schema-translator.jar
KT_JAR=bin/kotlin-runtime.jar

TRK_SRC=$(SRC)/SchemaTranslator/TranslatorKotlin.kt
TRJ_SRC=$(SRC)/SchemaTranslator/TranslatorJava.java
ST_SRC=$(SRC)/SchemaTranslator/SchemaTranslator.kt

# Note that these are a subset of generated classes; they are used as
# representatives of whether or not a rebuild is needed
TRK_CLASS=$(CP)/SchemaTranslator/TranslatorKotlin.class
TRJ_CLASS=$(CP)/SchemaTranslator/TranslatorJava.class
ST_CLASS=$(CP)/SchemaTranslator/SchemaTranslatorPackage.class


all :	jar
jar : $(ST_JAR)
	@echo "Build complete. You should now be able to run the SchemaTranslator with:"
	@echo "java -jar" $(ST_JAR)

clean :
	rm -rf $(CP) $(ST_JAR) $(KT_JAR)

prepare :
	mkdir -p $(CP)

$(TRK_CLASS) :	$(TRK_SRC)
	$(KOTLINC)  -cp $(CP) -d $(CP) $(TRK_SRC)

$(TRJ_CLASS) :	$(TRK_CLASS) $(TRJ_SRC)
	$(JAVAC) -cp $(CP):$(KOTLIN_RUNTIME) -d $(CP) $(TRJ_SRC)

$(ST_CLASS) :	$(TRJ_CLASS) $(TRK_CLASS) $(ST_SRC)
	$(KOTLINC)  -cp $(CP) -d $(CP) $(ST_SRC)

$(ST_JAR) :	$(TRK_CLASS) $(TRJ_CLASS) $(ST_CLASS)
	mkdir -p bin
	cp $(KOTLIN_RUNTIME) $(KT_JAR)
	$(JAR) cfm $(ST_JAR) META-INF/MANIFEST.MF -C $(CP) SchemaTranslator
# Note: the Java-based jar and fastjar appear to handle options differently.
# E.g. fastjar accepts c -f $ST_JAR -m META-INF/MANIFEST.MF -C ...
