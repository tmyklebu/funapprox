CFLAGS = -O3 -g -ansi -Wall -Wshadow -W -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wpointer-arith -Wnested-externs -Wundef -Wcast-qual -Wcast-align -Wwrite-strings 

demo_qs: demo_qs.o
	gcc $(CFLAGS) -o demo_qs demo_qs.o ../lib/QSopt_ex.a -lm -lgmp

# DO NOT DELETE THIS LINE -- make depend depends on it.

demo_qs.o:  demo_qs.c  ../include/QSopt_ex.h

