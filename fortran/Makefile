CC= ifort
LIBS = 

lyap: lyapJuanjo2015font3Marzo.for
	$(CC) -g -c lyapJuanjo2015font3Marzo.for -o main.o
	$(CC) -o lyap main.o $(LIBS)

clean: 
	rm *.o
	rm lyap
