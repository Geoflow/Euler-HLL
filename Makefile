# Compilateur utilisé
CC=g++

# Options en mode optimisé - La variable DEBUG est définie comme fausse
OPTIM_FLAG = -O3 -DNDEBUG 
# Options en mode debug - La variable est DEBUG est définie comme vraie
DEBUG_FLAG = -g -DDEBUG -Wall --leak-check=full 

CXX_FLAGS= $(OPTIM_FLAG)

# Librairies

LIB=-L/usr/lib -lm -fopenmp 
INC=-I/usr/include -lm -fopenmp 


# Le nom de l'exécutable
PROG = exe

# Les fichiers source à compiler
SRC = main.cpp euler.cpp mesh.cpp cell.cpp


# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

# La commande complète : compile seulement si un fichier a été modifié
#$(PROG) : $(SRC)$(CC) $(SRC) $(CXX_FLAGS) -o $(PROG)

exe : main.o mesh.o euler.o
	g++  main.cpp mesh.cpp euler.cpp cell.cpp  $(CXXFLAGS) -lm -fopenmp  -o   exe
main.o : main.cpp mesh.h euler.h 
	g++ -o main.o $(LIB) $(INC) -c main.cpp $(CXXFLAGS) -lm -fopenmp 
cell.o: cell.cpp cell.h
	g++ -c cell.cpp  $^ $(CXXFLAGS)  -o cell.o -lm -fopenmp 
mesh.o: mesh.cpp mesh.h cell.h
	g++  -c mesh.cpp $(CXXFLAGS)  -o mesh.o	-lm -fopenmp 
euler.o: euler.cpp euler.h  
	g++ -c euler.cpp $(CXXFLAGS)  -o euler.o -lm -fopenmp 

clean :
	rm *.o exe *.txt
