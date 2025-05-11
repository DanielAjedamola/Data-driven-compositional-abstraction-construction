#
# optimiser
GRBROOT = /Library/gurobi1000
PLATFORM = macos_universal2
TWOUP    = $(GRBROOT)/$(PLATFORM)
INC      = $(TWOUP)/include/
CPP      = g++
CARGS    = -m64 -g
CPPLIB   = -L$(TWOUP)/lib -lgurobi_c++ -lgurobi100

#
GRBAPP   = dotnet
DOTNETFRAMEWORK = --framework=netcoreapp6.0

# compiler
#
CCL        = g++
CC       	  = clang++
CXXLFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG -DSCOTS_BDD

#
# scots 
#
SCOTSROOT		= ../../..
SCOTSINC		= -I$(SCOTSROOT)/bdd -I$(SCOTSROOT)/srcc -I$(SCOTSROOT)/utils

#
# cudd 
#
CUDDPATH		=  $(SCOTSROOT)/cudd-3.0.0
CUDDINC 		= -I$(CUDDPATH)
CUDDLIBS		= -lcudd 
CUDDLPATH   	= -L$(CUDDPATH)/lib

#
# cudd 
#
CUDDLINC 		= -I$(CUDDPATH)
CUDD		    =  $(CUDDLINC) -L$(CUDDPATH)/lib -lcudd


.PHONY: consensus

TARGET = consensus

all: $(TARGET).cc
	$(CC) $(CARGS) $(CXXFLAGS) -o $@ $< -I$(INC) $(CUDDINC) $(SCOTSINC) $(CPPLIB) -lm $(CUDDLPATH) $(CUDDLIBS)

$(TARGET): $(TARGET).o
	$(CC) $(CXXFLAGS) $(SCOTSINC) -o $(TARGET) $(TARGET).o $(CUDDLPATH) $(CUDDLIBS)
# consensus:
# 	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) consensus.cc -o consensus.o

# .PHONY: M_bar

# TARGET1 = M_bar

# all1: $(TARGET1).cc
# 	$(CC) $(CARGS) $(CXXFLAGS) -o $@ $< -I$(INC) $(CUDDINC) $(SCOTSINC) $(CPPLIB) -lm $(CUDDLPATH) $(CUDDLIBS)

# $(TARGET1): $(TARGET1).o
# 	$(CC) $(CXXFLAGS) $(SCOTSINC) -o $(TARGET1) $(TARGET1).o $(CUDDLPATH) $(CUDDLIBS)



clean:
	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm *.dSYM; \
	if [ -d $(GRBAPP) ]; then \
		cd $(GRBAPP); \
		find . ! \( -name "gurobi*.csproj" -o -name . \) -exec rm -rf {} +; \
	fi

clean:
	rm consensus

# clean:
# 	rm M_bar