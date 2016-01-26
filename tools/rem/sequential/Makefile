
rem: replica.o main.o
	g++ -g -o rem replica.o main.o

main.o: main.cc
	g++ -g -c -Wall main.cc

replica.o: replica.cc replica.h
	g++ -g -c -Wall replica.cc

clean:
	@- rm -fv *.o

realclean: clean
	@- rm -fv rem

.DELETE_ON_ERROR:

