#MAKEFILE

##################################################################

DEBUG = 1

##################################################################

# compiler
# (CXX should be considered the c++ compiler automatically; like g++ on GNU/Linux)
CCC = $(CXX)
# CCC = clang++-10
CCFLAGS = -O3 -m64 -fPIC -fexceptions -DIL_STD -std=c++2a
# CCFLAGS = -Wall -std=c++2a -g
TARGETDIR=src

# adds directory the compiler looks for "gurobi_c++.h"
GRB_INCS = -I$(GUROBI_HOME)/include
# adds directory  the linker looks for libraries "libgurobi90.so"
# (or "libgurobi_c++.a" which might need to be recompiled)
GRB_LIBS = -L$(GUROBI_HOME)/lib -lgurobi_g++8.5 -lgurobi110

CCFLAGS += $(GRB_INCS)

ifeq ($(DEBUG),1)
	CCFLAGS += -O0 -g3 -DDEBUG
else
	CCFLAGS += -O3 -g
endif

OBJS2 = \
        $(TARGETDIR)/Data.o \
        $(TARGETDIR)/Environment.o \
        $(TARGETDIR)/WQAssign.o \
        $(TARGETDIR)/main.o

##################################################################

# dummy rule so "make" uses "all" as default target
# (needs to be first rule in the makefile)
.PHONY: default
default: all

##################################################################

# Pattern rules
# "$<" : name of first prerequisite
# "$@" : filename of target
# "$^" : all prerequisites minus copies
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp
	$(CCC) $(CCFLAGS) -c $< -o $@
	
##################################################################

.PHONY: all
all : onlineAssignment

$(TARGETDIR)/WQAssign.o: $(TARGETDIR)/WQAssign.cpp $(TARGETDIR)/Data.o
	$(CCC) $(CCFLAGS) -c $< -o $@

onlineAssignment: $(OBJS2)
	$(CCC) $(CCFLAGS) $^ -o $@ $(GRB_LIBS)

# test: onlineAssignment
# 	./onlineAssignment ../../../

##################################################################

RM = rm -rf
.PHONY: clean
clean:
	$(RM) $(TARGETDIR)/*.o $(TARGETDIR)/*.a onlineAssignment

#END
