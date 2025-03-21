all: libft.so

libft.so: ft_iteration.h ft_iteration.cpp comb.hpp comb.cpp
	g++ -O3 -c -Wall -Werror -fPIC ft_iteration.cpp comb.cpp
	g++ -O3 -shared -o libft.so ft_iteration.o comb.o

test-libft: libft.so ft_iteration.h ft_iteration.cpp test-libft.c
	g++ -O3 -Wall -Werror -Lft -o test-libft ft_iteration.cpp comb.cpp test-libft.c
	./test-libft

test-comb: comb.hpp comb.cpp test-comb.cpp
	g++ -O3 -Wall -Werror -o test-comb comb.cpp test-comb.cpp
	./test-comb

format:
	clang-format-20 -i ft_iteration.h ft_iteration.cpp comb.hpp comb.cpp test-libft.c test-comb.cpp

clean:
	rm -f *.o *.so test-libft test-comb
