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

clean:
	/bin/rm -f -r *.o
