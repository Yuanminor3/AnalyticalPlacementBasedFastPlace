#CFLAGS=-O3
CFLAGS=-g

all: outdir suraj_parser.o main.o
	g++ $(CFLAGS) -o suraj_parser suraj_parser.o main.o

suraj_parser.o: suraj_parser.cpp
	g++ $(CFLAGS) -c suraj_parser.cpp

main.o: main.cpp suraj_parser.h
	g++ $(CFLAGS) -c main.cpp

clean: 
	rm -f *.o suraj_parser
	rm -rf output

.PHONY: outdir

outdir: 
	mkdir -p output

testtoy01:
	./suraj_parser toy01

testtoy02:
	./suraj_parser toy02

testibm01:
	./suraj_parser ibm01

testibm03:
	./suraj_parser ibm03

testibm10:
	./suraj_parser ibm010

testibm16:
	./suraj_parser ibm16

	


