all:
	@cc \
		-std=c99 -pedantic \
		-Wall -Wextra \
		-O2 -funroll-loops \
		rng.c \
		u512.s fp.s \
		mont.c \
		fp.c\
		mont_own.c\
		csidh.c \
		main.c \
		-o main
bench:
	@cc \
		-std=c99 -pedantic \
		-Wall -Wextra \
		-O2 -funroll-loops \
		rng.c \
		u512.s fp.s \
		fp.c\
		mont.c \
		mont_own.c\
		csidh.c \
		bench.c \
		-o bench
		
test:
	@cc \
		-std=c99 -pedantic \
		-Wall -Wextra \
		-O2 -funroll-loops \
		rng.c \
		u512.s fp.s \
		fp.c\
		mont.c \
		mont_own.c\
		csidh.c \
		main_test.c \
		-o test

me_a_nice_signature_please:
	@cc \
		-std=c99 -pedantic \
		-Wall -Wextra \
		-O2 -funroll-loops \
		rng.c \
		u512.s fp.s \
		mont.c \
		fp.c\
		mont_own.c\
		csidh.c \
		tools_sign.c \
		sign.c \
		main_sign.c \
		-o main_sign

debug:
	cc \
		-std=c99 -pedantic \
		-Wall -Wextra \
		-g \
		rng.c \
		u512.s fp.s \
		mont.c \
		csidh.c \
		main.c \
		-o main

clean:
	rm -f main bench test

