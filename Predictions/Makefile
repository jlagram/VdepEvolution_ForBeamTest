GCC = g++
LD = g++
CXXFLAGS = $(shell root-config --cflags) 
LDFLAGS=-g $(shell root-config --ldflags)
LIBS = $(shell root-config --libs)
SRCS = ComputeVdepEvolution.C
OBJS = ComputeVdepEvolution.o
PROGRAM = Exec_ComputeVdepEvolution


all:	compile link

compile:	$(OBJS)

link:	$(OBJS)
	@echo "Linking object files ..."
	$(LD) $(LDFLAGS) $(LIBS) $(OBJS) -o $(PROGRAM)


clean:; @rm -f $(OBJS) $(PROGRAM) *~ */*~

