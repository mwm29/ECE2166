#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
//#include <mpi.h>

void print_usage(char *filename){
    printf("\nIncorrect number of args:\n");
    printf("usage:\n");
    printf("\t%s <x_max> <delta_x> <boundary_condition> <num_iter>\n\n", filename);
}

void validate_inputs(int argc, char **argv){
    if(argc < 5){
        print_usage(argv[0]);
        exit(1);
    }
}

int main ( int argc, char **argv ){
    validate_inputs(argc, argv);
    clock_t start, end;
    double cpu_time_used;

    start = clock();

    double x_max = atof(argv[1]);
    double delta_x = atof(argv[2]);
    double b_cond = atof(argv[3]);
    int num_iter = atoi(argv[4]);

    int num_points = x_max/delta_x+1;

    double *u = calloc(num_points,sizeof(double));
    double *u_guess = calloc(num_points,sizeof(double));
    double *temp;

    u_guess[0] = b_cond;
    u[0] = b_cond;
    for(int iter = 0; iter < num_iter; iter++){
        for(int i = 1; i < num_points; i++){
            u[i] = (u_guess[i-1] + delta_x*u_guess[i]*u_guess[i])/(1+2*delta_x*u_guess[i]);
        }
        temp = u;
        u = u_guess;
        u_guess = temp;
    }

    // FILE *outfile;
    // outfile = fopen("./results.txt", "w");
    //
    // if (!outfile) {
    //     fprintf(stderr, "Something went wrong when opening outfile\n");
    //     exit(3);
    // }

    fprintf(stdout, "%s\t%s\n", "x", "result");
    for(int i = 0; i < num_points; i++){
        fprintf(stdout, "%g\t%g\n", i*delta_x,u_guess[i]);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    free(u);
    free(u_guess);
    //fclose(outfile);
    fprintf(stdout, "Took %g seconds.\n",cpu_time_used);
}
