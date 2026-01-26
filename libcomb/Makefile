all: libft.so

ft_iteration.o: ft_iteration.h ft_iteration.c
	gcc -O3 -c -Wall -Werror -fPIC ft_iteration.c

cnot_plus_t_iteration.o: cnot_plus_t_iteration.c cnot_plus_t_iteration.h
	gcc -O3 -c -Wall -Werror -fPIC cnot_plus_t_iteration.c

libft.so: ft_iteration.o cnot_plus_t_iteration.o
	gcc -O3 -shared -o libft.so ft_iteration.o cnot_plus_t_iteration.o

clean:
	rm -f *.o libft.so test

test: clean all
	gcc test.c -L. -lft -o test
	LD_LIBRARY_PATH=. ./test
