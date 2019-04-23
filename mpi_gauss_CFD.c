#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <mpi.h>

int sparse_find(const int *S_row, const int *S_col, int S_size, int row, int col){
    for(int i = 0; i < S_size; i++){
        if(S_row[i] == row && S_col[i] == col){
            return i;
        }
    }
    return -1;
}
double sparse_row_sum(const int *S_row, const double *S, int S_size, int row){
    double sum = 0;
    for(int i = 0; i < S_size; i++){
        if(S_row[i] == row){
            sum += S[i];
        }
    }
    return sum;
}
void gaussian(double *A, double *b, int nrows, int world_rank, int world_size){
    //forward elimination
    for(int i = 0; i < nrows; i++){
        MPI_Bcast(&(A[i*nrows]), nrows, MPI_DOUBLE, i%world_size, MPI_COMM_WORLD);
        MPI_Bcast(&(b[i]), 1, MPI_DOUBLE, i%world_size, MPI_COMM_WORLD);
        //fprintf(stdout,"proc%d start: %d\n",world_rank,i-i%world_size+world_rank);
        for(int j = i-i%world_size+world_rank; j < nrows; j+= world_size){
            if(i<j){
                double scale = A[j*nrows+i]/A[i*nrows+i];
                //fprintf(stdout,"proc%d scale: %g\n",world_rank,scale);
                if(fabs(scale) > 10e-7){//only do subtraction if scale is not zero
                    for(int k = 0; k < nrows; k++){
                        A[j*nrows+k] -= A[i*nrows+k]*scale;
                    }
                    b[j] -= b[i]*scale;
                }
            }
        }
    }
    //back substitution
    for(int i = nrows-1; i >= 0; i--){
        b[i] = b[i] / A[i*nrows+i];
        A[i*nrows+i] = 1;
        //MPI_Bcast(&(A[i*nrows+i]), 1, MPI_DOUBLE, i%world_size, MPI_COMM_WORLD);
        MPI_Bcast(&(b[i]), 1, MPI_DOUBLE, i%world_size, MPI_COMM_WORLD);
        for(int j = i-i%world_size+world_rank; j >= 0; j-= world_size){
            if(i>j){
                double scale = A[j*nrows+i]/A[i*nrows+i];
                //fprintf(stdout,"proc%d %d scale: %g\n",world_rank,j,scale);
                A[j*nrows+i] = 0;
                //fprintf(stdout,"proc%d set A[%d] to 0 %g\n",world_rank,j*nrows+i,A[j*nrows+i]);
                b[j] -= b[i]*scale;
            }
        }
    }
}

int main ( int argc, char **argv ){
    //clock_t start, end;
    //double cpu_time_used;
    //start = clock();
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    double t1, t2;
    t1 = MPI_Wtime();

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //space variables (ENTER)
    int nx = 100;           //number of columns
    int ny = 100;           //number of rows
    //nx = 10;ny = 10;
    int nxy=nx*ny;
    //discretization variables (ENTER)
    double dx = 1;       //x-grid size
    double dy = 1;       //y-grid size
    double dxi = 1/dx;
    double dyi = 1/dy;
    //fluid density & viscosity (ENTER)
    double density = 1000;
    double viscosity = 1;
    double ratio = density/viscosity;       //ratio of density and dynamic viscosity

    //size and dimension of pressure and velocity (AUTO)
    double *p = calloc(nxy,sizeof(double));
    double *u = calloc((nx+1)*(ny+1),sizeof(double));
    for(int i = 0; i < (nx+1)*(ny+1); i++){
        u[i] = 0.1;
    }
    double *v = calloc((nx+2)*(ny+1),sizeof(double));
    double *residual = calloc(nxy,sizeof(double));
    double *dp = calloc(nxy,sizeof(double));

    // p = zeros(Ny,Nx);      %pressure
    // u = zeros(Ny+1,Nx+1)+0.1;  %x-velocity
    // v = zeros(Ny+1,Nx+2);  %y-velocity
    // residual = zeros(Ny*Nx,1); %residuals from continuity
    // dp = zeros(Ny*Nx,1);  %changes in pressures

    //initial conditions (ENTER)
    //u = zeros(Ny+1,Nx+1)+0.1;   //constant value of 0.1 (initializes velocity field)

    //temporary variables (AUTO)
    double *u1 = calloc((nx+1)*(ny+1),sizeof(double));
    for(int i = 0; i < (nx+1)*(ny+1); i++){
        u1[i] = 0.1;
    }
    double *v1 = calloc((nx+2)*(ny+1),sizeof(double));
    //dp1=zeros(Nx,Ny);
    double *dp1 = calloc(nxy,sizeof(double));
    //residual1=zeros(Nx,Ny);
    double *residual1 = calloc(nxy,sizeof(double));

    //timestep value, relaxation factor, number of iterations (ENTER)
    double dt = 1;
    double relaxation_factor = 0.5;
    int total_iterations = 1000;
    //total_iterations = 1;
    //residual_max = zeros(total_iterations,1);
    double *residual_max = calloc(total_iterations,sizeof(double));

    //check CFL criteria (CHECK!)
    //CFL_x = max(max(u))*dt/dx;
    //CFL_y = max(max(u))*dt/dy;

    t2 = MPI_Wtime();
    fprintf(stdout, "Setup finished at %g seconds.\n",t2-t1);

    //calculate sparse matrix (AUTO)
    double J_a = 2*(1/(dx*dx)+1/(dy*dy));
    double J_b = -1/(dy*dy);
    double J_c = -1/(dx*dx);

    double *J = calloc(nxy*nxy,sizeof(double));
    double *J_cpy = calloc(nxy*nxy,sizeof(double));
    for(int i = 0; i < nxy-1; i++){
        J[i*nxy+i+1] = J_b/J_a;
    }
    for(int i = 1; i < nxy; i++){
        J[i*nxy+i-1] = J_b/J_a;
    }
    for(int i = 0; i < nxy-nx; i++){
        J[i*nxy+i+nx] = J_c/J_a;
    }
    for(int i = 0; i < nxy-nx; i++){
        J[(i+nx)*nxy+i] = J_c/J_a;
    }
    for(int i = 1; i < ny; i++){
        J[(i*nx)*nxy+i*nx-1] = 0.0;
        J[(i*nx-1)*nxy+i*nx] = 0.0;
    }
    double sum;
    for(int i = 0; i < nxy; i++){
        sum = 0.0;
        for(int j = 0; j < nxy; j++){
            sum += J[i*nxy+j];
        }
        J[i*nxy+i] = -sum;
    }



    t2 = MPI_Wtime();
    fprintf(stdout, "J matrix built at %g seconds.\n",t2-t1);

    int *send_sizesu = calloc(world_size,sizeof(int));
    int *send_sizesv = calloc(world_size,sizeof(int));
    int *send_offsetsu = calloc(world_size,sizeof(int));
    int *send_offsetsv = calloc(world_size,sizeof(int));
    int leftovers = (ny-1) - ((ny-1) / world_size)*world_size;
    for(int ii = 0; ii < world_size; ii++){
        send_sizesu[ii] = (ny-1) / world_size;
        send_sizesv[ii] = (ny-1) / world_size;
        if(ii < leftovers){
            send_sizesu[ii] += 1;
            send_sizesv[ii] += 1;
        }
        // if(ii == 0){
        //     send_sizesu[ii] -= 1;
        //     send_sizesv[ii] -= 1;
        // }
        send_sizesu[ii] *= (nx+1);
        send_sizesv[ii] *= (nx+2);
        for(int jj = ii+1; jj < world_size; jj++){
            send_offsetsu[jj] += send_sizesu[ii];
            send_offsetsv[jj] += send_sizesv[ii];
        }
    }
    int start_row = send_offsetsu[world_rank] / (nx+1)+1;
    int rows_to_do = send_sizesu[world_rank] / (nx+1);
    int end_row = start_row + rows_to_do;
    // if(world_rank == 0 || world_size == 1){
    //     start_row++;
    //     rows_to_do--;
    // }
    //fprintf(stderr,"proc%d start %d rtd %d end %d\n",world_rank,start_row,rows_to_do,end_row);

    //calculate velocity field using Navier-Stokes equations
    for(int j = start_row; j < end_row; j++){
        for(int i = 1; i < nx; i++){
            double a1 = -(u[j*(nx+1)+i+1]*u[j*(nx+1)+i+1] - u[j*(nx+1)+i-1]*u[j*(nx+1)+i-1])/(2*dx) - (u[(j+1)*(nx+1)+i]*(v[(j+1)*(nx+2)+i+1]+v[(j+1)*(nx+2)+i]) - u[j*(nx+1)+i]*(v[j*(nx+2)+i+1]+v[j*(nx+2)+i]))/(4*dy);
            double a3 = (u[j*(nx+1)+i+1]-2*u[j*(nx+1)+i]+u[j*(nx+1)+i-1])/(dx*dx);
            double a4 = (u[(j+1)*(nx+1)+i]-2*u[j*(nx+1)+i]+u[(j-1)*(nx+1)+i])/(dy*dy);

            double A = a1+(a3+a4)/ratio;

            u1[j*(nx+1)+1] = u[j*(nx+1)+1] + dt*(A-(p[j*(nx)+i]-p[j*(nx)+i-1])/(dx));
        }
    }

    for(int j = start_row; j < end_row; j++){
        for(int i = 1; i <= nx; i++){
            double b1 = -(v[(j+1)*(nx+2)+i]*v[(j+1)*(nx+2)+i] - v[(j-1)*(nx+2)+i]*v[(j-1)*(nx+2)+i])/(2*dx) - ((v[j*(nx+2)+i+1]*(u[(j-1)*(nx+1)+i]+u[j*(nx+1)+i]))-(v[j*(nx+2)+i-1]*(u[(j-1)*(nx+1)+i-1]+u[j*(nx+1)+i-1])))/(4*dx);
            double b3 = (v[(j+1)*(nx+2)+i]-2*v[j*(nx+2)+i]+v[(j-1)*(nx+2)+i])/(dy*dy);
            double b4 = (v[j*(nx+2)+i+1]-2*v[j*(nx+2)+i]+v[j*(nx+2)+i-1])/(dx*dx);

            double B = b1+(b3+b4)/ratio;

            v1[j*(nx+2)+i] = v1[j*(nx+2)+i] + dt*(B-(p[j*(nx)+i-1]-p[(j-1)*(nx)+i-1])/dy);
        }
    }

    MPI_Allgatherv(&(u1[start_row*(nx+1)]), send_sizesu[world_rank], MPI_DOUBLE, &(u1[nx+1]), send_sizesu, send_offsetsu, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(&(v1[start_row*(nx+2)]), send_sizesv[world_rank], MPI_DOUBLE, &(v1[nx+2]), send_sizesv, send_offsetsv, MPI_DOUBLE, MPI_COMM_WORLD);


    //apply boundary conditions
    // flow around square
    for(int j = nx/2-nx/5+1; j < nx/2+nx/5; j++){
        for(int i = ny/2-ny/5+1; i < ny/2+ny/5; i++){
            u1[j*(nx+1)+i] = 0.0;
            v1[j*(nx+2)+i] = 0.0;
        }
    }

    int uv_start_row = start_row;
    int uv_end_row = end_row;
    int *send_sizes = calloc(world_size,sizeof(int));
    int *send_offsets = calloc(world_size,sizeof(int));
    leftovers = (ny) - ((ny) / world_size)*world_size;
    for(int ii = 0; ii < world_size; ii++){
        send_sizes[ii] = (ny) / world_size;
        if(ii < leftovers){
            send_sizes[ii] += 1;
        }
        // if(ii == 0){
        //     send_sizesu[ii] -= 1;
        //     send_sizesv[ii] -= 1;
        // }
        send_sizes[ii] *= (nx);
        for(int jj = ii+1; jj < world_size; jj++){
            send_offsets[jj] += send_sizes[ii];
        }
    }
    start_row = send_offsets[world_rank] / (nx);
    rows_to_do = send_sizes[world_rank] / (nx);
    end_row = start_row + rows_to_do;
    //fprintf(stdout, "proc%d start %d end %d\n",world_rank,start_row,end_row);
    int iteration;
    //gsl_vector_view b,x;
    //iterate for pressure and velocity corrections
    for(iteration = 0; iteration < total_iterations; iteration++){

        for(int j = start_row; j < end_row; j++){
            for(int i = 0; i < nx; i++){
                residual1[j*nx+i] = (u1[j*(nx+1)+i+1] - u1[j*(nx+1)+i] + v1[(j+1)*(nx+2)+i] - v1[j*(nx+2)+i])/(-J_a*dt);
                //fprintf(stdout,"%g\n",residual1[j*nx+i]);
                if(residual1[j*nx+i] > residual_max[iteration]){
                    residual_max[iteration] = residual1[j*nx+i];
                }
            }
        }
        MPI_Allgatherv(&(residual1[start_row*(nx)]), send_sizes[world_rank], MPI_DOUBLE, residual1, send_sizes, send_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &(residual_max[iteration]), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // x = gsl_vector_view_array(dp, nxy);
        // b = gsl_vector_view_array(residual1, nxy);
        // gsl_linalg_LU_solve(mat, perm, &b.vector, &x.vector);
        memcpy(J_cpy, J, nxy*nxy*sizeof(double));
        memcpy(dp, residual1, nxy*sizeof(double));
        gaussian(J_cpy, dp, nxy, world_rank, world_size);



        for(int j = uv_start_row; j < uv_end_row; j++){
            for(int i = 1; i < nx; i++){
                u1[j*(nx+1)+i] = u1[j*(nx+1)+i] + relaxation_factor*(dp[j*nx+i-1]-dp[j*nx+i])*dt/dx;
            }
        }


        for(int j = uv_start_row; j < uv_end_row; j++){
            for(int i = 1; i <= nx; i++){
                v1[j*(nx+2)+i] = v1[j*(nx+2)+i] + relaxation_factor*(dp[(j-1)*nx+i-1]-dp[j*nx+i-1])*dt/dy;
            }
        }


        for(int j = start_row; j < end_row; j++){
            for(int i = 0; i < nx; i++){
                p[j*nx+i] = p[j*nx+i] + relaxation_factor * dp[j*nx+i];
            }
        }

        MPI_Allgatherv(&(u1[uv_start_row*(nx+1)]), send_sizesu[world_rank], MPI_DOUBLE, &(u1[nx+1]), send_sizesu, send_offsetsu, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(&(v1[uv_start_row*(nx+2)]), send_sizesv[world_rank], MPI_DOUBLE, &(v1[nx+2]), send_sizesv, send_offsetsv, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(&(p[start_row*(nx)]), send_sizes[world_rank], MPI_DOUBLE, p, send_sizes, send_offsets, MPI_DOUBLE, MPI_COMM_WORLD);

        if(residual_max[iteration] < 0.00001*nx){
            break;
        }
    }//iteration ends
    fprintf(stdout,"Exited after %d iterations\n",iteration);




    t2 = MPI_Wtime();
    MPI_Finalize();

    if(world_rank == 0){
        FILE *out_fp = fopen("p_CFD_Res.csv", "w");
        for(int j = 0; j < ny; j++){
            for(int i = 0; i < nx; i++){
                fprintf(out_fp,"%g",p[j*nx+i]);
                if(i == nx-1){
                    fprintf(out_fp,"\n");
                }else{
                    fprintf(out_fp,",");
                }
            }
        }
        fclose(out_fp);

        out_fp = fopen("u_CFD_Res.csv", "w");
        for(int j = 0; j < ny+1; j++){
            for(int i = 0; i < nx+1; i++){
                fprintf(out_fp,"%g",u1[j*(nx+1)+i]);
                if(i == nx){
                    fprintf(out_fp,"\n");
                }else{
                    fprintf(out_fp,",");
                }
            }
        }
        fclose(out_fp);

        out_fp = fopen("v_CFD_Res.csv", "w");
        for(int j = 0; j < ny+1; j++){
            for(int i = 0; i < nx+2; i++){
                fprintf(out_fp,"%g",v1[j*(nx+2)+i]);
                if(i == nx+1){
                    fprintf(out_fp,"\n");
                }else{
                    fprintf(out_fp,",");
                }
            }
        }
        fclose(out_fp);

        out_fp = fopen("r_CFD_Res.csv", "w");
        for(int i = 0; i < iteration; i++){
            fprintf(out_fp,"%g",residual_max[i]);
            if(i == iteration-1){
                fprintf(out_fp,"\n");
            }else{
                fprintf(out_fp,",");
            }
        }
        fclose(out_fp);
    }

    fprintf(stdout, "Program ended at %g seconds.\n",t2-t1);
}
