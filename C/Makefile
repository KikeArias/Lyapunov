CC= gcc -g -DDEBUG
LIBS = -lm 

lyap: main.c
	$(CC) -c main.c -o main.o
	$(CC) -c memoryfun.c -o memoryfun.o
	$(CC) -c runge.c -o runge.o
	$(CC) -c duffing.c -o duffing.o
	$(CC) -c clambda.c -o clambda.o
	$(CC) -o lyap main.o memoryfun.o runge.o duffing.o clambda.o $(LIBS)

clean: 
	rm *.o
	rm lyap
