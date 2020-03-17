CC= gcc -DDEBUG
LIBS = -lm 

lyap: main.c
	$(CC) -g -c main.c -o main.o
	$(CC) -g -c memoryfun.c -o memoryfun.o
	$(CC) -g -c runge.c -o runge.o
	$(CC) -g -c duffing.c -o duffing.o
	$(CC) -g -c clambda.c -o clambda.o
	$(CC) -o lyap main.o memoryfun.o runge.o duffing.o clambda.o $(LIBS)

clean: 
	rm *.o
	rm lyap
