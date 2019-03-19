pca: lab2_io.c lab2_omp.c main_omp.c
	g++ -fopenmp -lm lab2_io.c lab2_omp.c main_omp.c -o pca

pcag: lab2_io.c lab2_omp.c main_omp.c
	g++ -g -fopenmp -lm lab2_io.c lab2_omp.c main_omp.c -o pcag

.PHONY: clean
	rm -f pca pcag