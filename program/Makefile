CC = g++
OBJ = fem.o

default : main

$(OBJ) : %.o : %.cpp %.h
	$(CC) -c $< -o $@ -I. -lm -ldl -Wall

main.o : main.cpp
	$(CC) -c -o main.o main.cpp -I. -lm -ldl -Wall

main : $(OBJ) main.o
	$(CC) -o main $(OBJ) main.o -I. -lm -ldl -Wall

clean :
	rm -f $(OBJ) main.o
	rm -f main

.PHONY : clean
