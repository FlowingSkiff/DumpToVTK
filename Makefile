OBJS	= main.o
SOURCE	= main.cpp
HEADER	= OutputVTK.cpp OutputVTK.h
OUT	= LammpsToVtk
CC	 = g++
FLAGS	 = -g -c -Wall
LFLAGS	 = 

all: $(OBJS)
		$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

main.o: main.cpp
		$(CC) $(FLAGS) main.cpp 

OutputVTK.o: OutputVTK.cpp
		$(CC) $(FLAGS) OutputVTK.cpp 


clean:
		rm -f $(OBJS) $(OUT)
