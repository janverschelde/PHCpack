qd_test:
	gcc -o /tmp/qd_test qd_test.c quad_double.c double_double.c -lm

dd_test:
	gcc -o /tmp/dd_test dd_test.c double_double.c -lm

dd_funtest:
	g++ -c double_double.c -o /tmp/double_double.o
	g++ -o /tmp/dd_funtest dd_funtest.cpp /tmp/double_double.o -lm

td_test:
	gcc -c double_double.c -o /tmp/double_double.o
	gcc -c triple_double.c -o /tmp/triple_double.o
	gcc -o /tmp/td_test td_test.c \
               /tmp/double_double.o /tmp/triple_double.o -lm

pd_test:
	gcc -c double_double.c -o /tmp/double_double.o
	gcc -c penta_double.c -o /tmp/penta_double.o
	gcc -o /tmp/pd_test pd_test.c \
               /tmp/double_double.o /tmp/penta_double.o -lm

od_test:
	gcc -c double_double.c -o /tmp/double_double.o
	gcc -c octo_double.c -o /tmp/octo_double.o
	gcc -o /tmp/od_test od_test.c \
               /tmp/double_double.o /tmp/octo_double.o -lm

da_test:
	gcc -c double_double.c -o /tmp/double_double.o
	gcc -c deca_double.c -o /tmp/deca_double.o
	gcc -o /tmp/da_test da_test.c \
               /tmp/double_double.o /tmp/deca_double.o -lm

hd_test:
	gcc -c double_double.c -o /tmp/double_double.o
	gcc -c hexa_double.c -o /tmp/hexa_double.o
	gcc -o /tmp/hd_test hd_test.c \
               /tmp/double_double.o /tmp/hexa_double.o -lm

clean:
	/bin/rm -f -r *.o
