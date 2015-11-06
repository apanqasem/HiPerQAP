CC=nvcc
CXX=g++

CXXFLAGS=-O2

INCPATH=-I.

ARCH=-arch sm_35
FLAGS=$(ARCH) -O3 -w -ccbin g++ -DPROFILE #--ptxas-options=-v

LFLAGS=-lcudadevrt 
DFLAGS= -DDEBUG

EXECS=2opt 2opt_shared 2opt_dyn 2opt_dyn_diverse tabu tabu_dyn anneal

all: $(EXECS)


#
# core algorithms 
#
2opt: 2opt.cu qapio.o 2optlib.o
	$(CC) -o 2opt $(INCPATH) $(FLAGS) 2opt.cu qapio.o 2optlib.o

anneal: anneal.cu qapio.o 2optlib.o
	$(CC) -o anneal $(INCPATH) $(FLAGS) anneal.cu qapio.o 2optlib.o

tabu: tabu.cu qapio.o 2optlib.o
	$(CC) -o tabu $(INCPATH) -rdc=true $(FLAGS) tabu.cu  $(LFLAGS) qapio.o 2optlib.o


#
# variants with dynamic threads and shared memory 
#
tabu_dyn: tabu_dyn.cu qapio.o  
	$(CC) -o tabu_dyn $(INCPATH) -rdc=true $(FLAGS) tabu_dyn.cu qapio.o $(LFLAGS)

2opt_shared: 2opt_shared.cu qapio.o 2optlib.o
	$(CC) -o 2opt_shared $(INCPATH) $(FLAGS) 2opt_shared.cu qapio.o 2optlib.o

2opt_dyn: 2opt_dyn.cu  qapio.o 2optlib.o
	$(CC) -o 2opt_dyn $(INCPATH) -rdc=true $(FLAGS) 2opt_dyn.cu qapio.o 2optlib.o $(LFLAGS)


2opt_dyn_diverse: 2opt_dyn_diverse.cu  qapio.o 2optlib.cu
	$(CC) -o 2opt_dyn_diverse $(INCPATH) -rdc=true $(FLAGS) 2opt_dyn_diverse.cu qapio.o 2optlib.o $(LFLAGS)

#
# builds with debug symbols 
#
2opt_D: 2opt.cu qapio.o 2optlib.o
	$(CC) -o 2opt $(INCPATH) $(FLAGS) $(DFLAGS) 2opt.cu qapio.o 2optlib.o

2opt_shared_D: 2opt_shared.cu qapio.o 2optlib.o
	$(CC) -o 2opt_shared $(INCPATH) $(FLAGS) $(DFLAGS) 2opt_shared.cu qapio.o 2optlib.o

%.o: %.cpp
	$(CXX) $(INCPATH) -c $(CXXFLAGS) $< -o $@
clean:
	rm -rf *~  *.o $(EXECS)


