utest : unitTest.o question1.o
		gcc unitTest.o question1.o -o utest -lm

unitTest.o : unitTest.c question1.c
		gcc -c unitTest.c -o unitTest.o -lm

question1.o : question1.c
		gcc -c question1.c -o question1.o -lm


vtest : verifTest.o question1.o
		gcc verifTest.o question1.o -o vtest -lm

verifTest.o : verifTest.c question1.c
		gcc -c verifTest.c -o verifTest.o -lm

question1.o : question1.c
		gcc -c question1.c -o question1.o -lm
