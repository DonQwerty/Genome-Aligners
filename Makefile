
    #######################################
    #     	    MAKEFILE				  #
    # 									  #
    # AUTHOR: Daniel Gimenez Llorente	  #
    # 	                              	  #
    #######################################

# DIRECTORIES
IDIR		= include
SDIR		= src
ODIR		= obj
BDIR		= bin

# FILES
TARGET		= aligners
SOURCES		:= $(wildcard $(SDIR)/*.cpp)
INCLUDES 	:= $(wildcard $(IDIR)/*.h)
OBJECTS  	:= $(SOURCES:$(SDIR)/%.cpp=$(ODIR)/%.o)
CHK_SOURCES	:= $(SOURCES) $(INCLUDES)

# VARIABLES
# Compiler
CC			= g++ -c
CFLAGS		= -Wall -I$(IDIR) -ansi -pedantic -std=c++11
# Linker
LINKER		= g++ -o
LFLAGS		= -Wall -I$(IDIR) -ansi -pedantic -std=c++11
# Others
rm			= rm -f

# Compiles all programs in project
compile: $(BDIR)/$(TARGET)
	@echo "INFO: All targets compiled"

$(BDIR)/$(TARGET): $(OBJECTS)
	@$(LINKER) $@ $(LFLAGS) $(OBJECTS)
	@echo "INFO: Linking completed"

$(OBJECTS):$(ODIR)/%.o: $(SDIR)/%.cpp
	@$(CC) $(CFLAGS) $< -o $@
	@echo "INFO: Compiled "$<" successfully"

# Displays the help for this makefile
help:
	@echo "TODO: help"

.PHONEY: clean
clean:
	@$(rm) $(OBJECTS)
	@echo "INFO: Cleaning completed"
