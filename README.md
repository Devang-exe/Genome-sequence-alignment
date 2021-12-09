# Genome-sequence-alignment
Genome sequence alignment using Smith waterman Algorithm. Parallel and distributed computing project (PDC)
Used Openmp


## To run the program
->serial <br />
gcc serial_sw_file.c -o serial_sw_file <br />
./serial_sw_file <br />
->parallel <br />
gcc parallel_sw_file.c -o parallel_sw_file -fopenmp -DDEBUG <br />
./parallel_sw_file 3 <br />
