# Make file for my_ldhot
# Author: Adam Auton

CPP = g++
CPPFLAGS = -Wall -Wextra -Ofast -m64 -mtune=native 
#CPPFLAGS = -Wall -Wextra -m64 -O3 
#CPPFLAGS = -g

all: my_ldhot my_ldhot_summary

my_ldhot: my_ldhot.cpp sim.cpp tools.cpp parameters.cpp output_log.cpp gpd_fit.cpp
	$(CPP) $(CPPFLAGS) my_ldhot.cpp parameters.cpp output_log.cpp tools.cpp sim.cpp gpd_fit.cpp -o my_ldhot

my_ldhot_summary: summary.cpp tools.cpp output_log.cpp gpd_fit.cpp
	$(CPP) $(CPPFLAGS) summary.cpp tools.cpp output_log.cpp gpd_fit.cpp -o my_ldhot_summary
	
clean:
	rm -rf *.o
	rm -rf *~
	rm my_ldhot my_ldhot_summary
