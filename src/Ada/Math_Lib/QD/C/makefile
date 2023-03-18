all: dd_test dd_funtest td_test qd_test pd_test od_test da_test hd_test

BIN=../../../../bin

qd_test:
	gcc -o $(BIN)/qd_test qd_test.c quad_double.c double_double.c -lm

dd_test:
	gcc -o $(BIN)/dd_test dd_test.c double_double.c -lm

dd_funtest:
	gcc -c double_double.c -o double_double.o
	gcc -o $(BIN)/dd_funtest dd_funtest.c double_double.o -lm

td_test:
	gcc -c double_double.c -o double_double.o
	gcc -c triple_double.c -o triple_double.o
	gcc -o $(BIN)/td_test td_test.c \
               double_double.o triple_double.o -lm

pd_test:
	gcc -c double_double.c -o double_double.o
	gcc -c penta_double.c -o penta_double.o
	gcc -o $(BIN)/pd_test pd_test.c \
               double_double.o penta_double.o -lm

od_test:
	gcc -c double_double.c -o double_double.o
	gcc -c octo_double.c -o octo_double.o
	gcc -o $(BIN)/od_test od_test.c \
               double_double.o octo_double.o -lm

da_test:
	gcc -c double_double.c -o double_double.o
	gcc -c deca_double.c -o deca_double.o
	gcc -o $(BIN)/da_test da_test.c \
               double_double.o deca_double.o -lm

hd_test:
	gcc -c double_double.c -o double_double.o
	gcc -c hexa_double.c -o hexa_double.o
	gcc -o $(BIN)/hd_test hd_test.c \
               double_double.o hexa_double.o -lm

clean:
	del *.o *~
