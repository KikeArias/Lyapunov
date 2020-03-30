CC= gfortran
LIBS = 

lyap: lyapJuanjo2015font3Marzo.for
	$(CC) -g -c clambda.for -o clambda.o
	$(CC) -g -c lyapJuanjo2015font3Marzo.for -o main.o
	$(CC) -o lyap main.o clambda.o $(LIBS)

clean: 
	rm *.o
	rm lyap
