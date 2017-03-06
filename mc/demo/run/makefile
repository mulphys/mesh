#!/bin/sh

SRCPATH   = ../src
JAVAC     = $(JAVA_HOME)/bin/javac 
JAVAC_OPT = -sourcepath $(SRCPATH) -d . 
SRC       = $(SRCPATH)/Mesh.java

#$(JAVAC) -sourcepath ../src -d . ../src/Mesh.java

all : $(SRCPATH)/*.java
	$(JAVAC) $(JAVAC_OPT) $(SRC)

clean :
	find . -type f \( -name "*.class" -o -name "*.[ou]" -o -name "*~" \) \
	  -exec rm -f {} \;
