# Make file for ldhot
# Author: Adam Auton

# Flag used to turn on multi-threading
ifndef MULTI
    MULTI = 0
endif

CPP = g++
#CPP = g++-mp-4.8
CPPFLAGS = -Wall -Wextra -O2 -m64 -mtune=native -std=c++11 
ifeq ($(MULTI), 1)
    CPPFLAGS += -fopenmp
endif
#CPPFLAGS = -Wall -Wextra -m64 -O2 
#CPPFLAGS = -g

all: ldhot ldhot_summary

ldhot: ldhot.cpp sim.cpp tools.cpp parameters.cpp output_log.cpp gpd_fit.cpp
	$(CPP) $(CPPFLAGS) ldhot.cpp parameters.cpp output_log.cpp tools.cpp sim.cpp gpd_fit.cpp -o ldhot

ldhot_summary: summary.cpp tools.cpp output_log.cpp 
	$(CPP) $(CPPFLAGS) summary.cpp tools.cpp output_log.cpp -o my_ldhot_summary
	
clean:
	rm -rf *.o
	rm -rf *~
	rm ldhot ldhot_summary
