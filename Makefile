CC = gcc
LIBS = -lgsl -lgslcblas -lm

default:
	$(CC) -ISWSH -I/usr/local/include -c *.c
	$(CC) -L/usr/local/lib -o SWSH *.o $(LIBS)
	rm *.o