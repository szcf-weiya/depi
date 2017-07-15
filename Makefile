CPP = g++
OFLAG = -o
GSLFLAG = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm
OPENMPFLAG = -fopenmp

source = main.cpp

all:\
	mpout \
	out \



mpout: $(source)
	$(CPP) $(source) $(OFLAG) mpout -O3 $(GSLFLAG) $(OPENMPFLAG)

out: $(source)
	$(CPP) $(source) $(OFLAG) out $(GSLFLAG) 

