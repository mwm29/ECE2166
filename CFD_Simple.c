#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
//#include <mpi.h>
#define USE_LIB
#ifdef USE_LIB
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#endif

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
#ifdef USE_LIB
void gsl_mat_inv(gsl_matrix *A, double **A_inv, int nxy){
    gsl_matrix * gA_t = gsl_matrix_alloc (nxy, nxy);
    gsl_matrix_transpose_memcpy (gA_t, A);


    gsl_matrix * U = gsl_matrix_alloc (nxy, nxy);
    gsl_matrix * V = gsl_matrix_alloc (nxy, nxy);
    gsl_vector * S = gsl_vector_alloc (nxy);

    gsl_vector * work = gsl_vector_alloc (nxy);
    gsl_linalg_SV_decomp (gA_t, V, S, work);
    gsl_vector_free(work);

    gsl_matrix_memcpy (U, gA_t);




    gsl_matrix * Sp = gsl_matrix_alloc (nxy, nxy);
    gsl_matrix_set_zero (Sp);
    for (int i = 0; i < nxy; i++)
    	gsl_matrix_set (Sp, i, i, gsl_vector_get(S, i));	// Vector 'S' to matrix 'Sp'

    gsl_permutation * p = gsl_permutation_alloc (nxy);
    int signum;
    gsl_linalg_LU_decomp (Sp, p, &signum);				// Computing the LU decomposition

    // Compute the inverse like in the MATLAB script

    gsl_matrix * SI = gsl_matrix_calloc (nxy, nxy);

    for (int i = 0; i < nxy; i++) {
      if (gsl_vector_get (S, i) > 0.0000000001)
        gsl_matrix_set (SI, i, i, 1.0 / gsl_vector_get (S, i));
    }

    gsl_matrix * VT = gsl_matrix_alloc (nxy, nxy);
    gsl_matrix_transpose_memcpy (VT, V);					// Tranpose of V


    //THE PSEUDOINVERSE//
    //----------------------------------------------------------
    //Computation of the pseudoinverse of trans(A) as pinv(A) = U·inv(S).trans(V)   with trans(A) = U.S.trans(V)
    //----------------------------------------------------------
    gsl_matrix * SIpVT = gsl_matrix_alloc (nxy, nxy);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,				// Calculating  inv(S).trans(V)
                	1.0, SI, VT,
                	0.0, SIpVT);


    gsl_matrix * pinv = gsl_matrix_alloc (nxy, nxy);	// Calculating  U·inv(S).trans(V)
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                	1.0, U, SIpVT,
                	0.0, pinv);

    //copy to c matrix
    for(int i = 0; i < nxy; i++){
        gsl_vector_view row = gsl_matrix_row(pinv, i);
        for(int j = 0; j < nxy; j++){
            (*A_inv)[i*nxy+j] = row.vector.data[j];
            //fprintf(stdout,"%4.2g\t",row.vector.data[j]);
        }
        //fprintf(stdout,"\n");
    }


    gsl_matrix_free(pinv);
    gsl_matrix_free(SIpVT);
    gsl_matrix_free(VT);
    gsl_matrix_free(SI);
    gsl_permutation_free(p);
    gsl_matrix_free(Sp);
    gsl_vector_free(S);
    gsl_matrix_free(V);
    gsl_matrix_free(U);
    gsl_matrix_free(gA_t);
}
#endif

int main ( int argc, char **argv ){
    clock_t start, end;
    double cpu_time_used;
    start = clock();
    //space variables (ENTER)
    int nx = 50;           //number of columns
    int ny = 50;           //number of rows
    //nx = 10;ny = 10;
    int nxy=nx*ny;
    //discretization variables (ENTER)
    double dx = 1;       //x-grid size
    double dy = 1;       //y-grid size
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

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    fprintf(stdout, "Setup finished at %g seconds.\n",cpu_time_used);

    //calculate sparse matrix (AUTO)
    double J_a = 2*(1/(dx*dx)+1/(dy*dy));
    double J_b = -1/(dy*dy);
    double J_c = -1/(dx*dx);
    //int *J_row = malloc(((nx-2)*(ny-2)*4+nx*ny)*sizeof(int));
    //int *J_col = malloc(((nx-2)*(ny-2)*4+nx*ny)*sizeof(int));
    //double *J = malloc(((nx-2)*(ny-2)*4+nx*ny)*sizeof(double));//this size was too small
    #ifdef USE_SPARSE
    int *J_row = malloc(nxy*5*sizeof(int));
    int *J_col = malloc(nxy*5*sizeof(int));
    double *J = malloc(nxy*5*sizeof(double));
    for(int i = 0; i < nxy*5; i++){
        J_row[i] = -1;
        J_col[i] = -1;
    }
    int J_num = 0;
    //J=spalloc(Nx*Ny,Nx*Ny,(Nx-2)*(Ny-2)*4+Nx*Ny);

    //for i=1:Nx*Ny-1
    //    J(i,i+1)=J_b/J_a;
    //end
    for(int i = 0; i < nxy-1; i++){
        // if(sparse_find(J_row, J_col, J_num, i, i+1) != -1){
        //     fprintf(stderr,"repeat1\n");
        // }
        J_row[J_num] = i;
        J_col[J_num] = i+1;
        J[J_num] = J_b/J_a;
        J_num++;
    }
    //for i=2:Nx*Ny
        //J(i,i-1)=J_b/J_a;
    //end
    for(int i = 1; i < nxy; i++){
        // if(sparse_find(J_row, J_col, J_num, i, i-1) != -1){
        //     fprintf(stderr,"repeat2\n");
        // }
        J_row[J_num] = i;
        J_col[J_num] = i-1;
        J[J_num] = J_b/J_a;
        J_num++;
    }
    //for i=1:Nx*Ny-Nx
        //J(i,i+Nx)=J_c/J_a;
    //end
    for(int i = 0; i < nxy-nx; i++){
        // if(sparse_find(J_row, J_col, J_num, i, i+nx) != -1){
        //     fprintf(stderr,"repeat3\n");
        // }
        //fprintf(stdout,"%d %d\n",i,i+nx);
        J_row[J_num] = i;
        J_col[J_num] = i+nx;
        J[J_num] = J_c/J_a;
        J_num++;
    }
    //for i=1:Nx*Ny-Nx
        //J(i+Nx,i)=J_c/J_a;
    //end
    for(int i = 0; i < nxy-nx; i++){
        // if(sparse_find(J_row, J_col, J_num, i+nx, i) != -1){
        //     fprintf(stderr,"repeat4\n");
        // }
        J_row[J_num] = i+nx;
        J_col[J_num] = i;
        J[J_num] = J_c/J_a;
        J_num++;
    }
    //for i=1:Ny-1
    //J(i*Nx+1,i*Nx)=0;
    //J(i*Nx,i*Nx+1)=0;
    //end
    for(int i = 1; i < ny; i++){
        int index = sparse_find(J_row, J_col, J_num, i*nx, i*nx-1);
        // if(index == -1){
        //     fprintf(stderr,"not found1\n");
        // }
        J_num--;
        J_row[index] = J_row[J_num];
        J_col[index] = J_col[J_num];
        J[index] = J[J_num];
        J_row[J_num] = -1;
        J_col[J_num] = -1;

        index = sparse_find(J_row, J_col, J_num, i*nx-1, i*nx);
        // if(index == -1){
        //     fprintf(stderr,"not found2\n");
        // }
        J_num--;
        J_row[index] = J_row[J_num];
        J_col[index] = J_col[J_num];
        J_row[J_num] = -1;
        J_col[J_num] = -1;
        J[index] = J[J_num];
    }

    //for i=1:Nx*Ny
        //J(i,i)=-sum(J(i,:));
    //end
    for(int i = 0; i < nxy; i++){
        // if(sparse_find(J_row, J_col, J_num, i, i) != -1){
        //     fprintf(stderr,"repeat5\n");
        // }
        J_row[J_num] = i;
        J_col[J_num] = i;
        J[J_num] = -1.0*sparse_row_sum(J_row, J, J_num, i);
        J_num++;
    }
    // fprintf(stdout,"%d %d\n", 4*nxy, J_num);
    // for(int i = 0; i < nxy; i++){
    //     for(int j = 0; j < nxy; j++){
    //         int index = sparse_find(J_row, J_col, J_num, i, j);
    //         //if((index != -1) && (J[index] > 0.00000001 || J[index] < -0.00000001)){
    //         if((index != -1) ){
    //             fprintf(stdout,"%4.2g\t",J[index]);
    //         }else{
    //             fprintf(stdout,"0.0\t");
    //         }
    //     }
    //     fprintf(stdout,"\n");
    // }
    // exit(0);
    //%spy(J)   %for checking sparse matrix
    #endif
    #ifndef USE_SPARSE
    #ifndef USE_LIB
    double *J = calloc(nxy*nxy,sizeof(double));
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
    #endif
    #ifdef USE_LIB
    double *J_inv = calloc(nxy*nxy,sizeof(double));
    gsl_matrix *mat = gsl_matrix_calloc(nxy, nxy);
    //gsl_matrix *inv = gsl_matrix_calloc(nxy, nxy);
    //gsl_permutation * perm = gsl_permutation_alloc(nxy);
    //int signum;
    for(int i = 0; i < nxy-1; i++){
        //J[i*nxy+i+1] = J_b/J_a;
        gsl_matrix_set(mat, i, i+1, J_b/J_a);
    }
    for(int i = 1; i < nxy; i++){
        //J[i*nxy+i-1] = J_b/J_a;
        gsl_matrix_set(mat, i, i-1, J_b/J_a);
    }
    for(int i = 0; i < nxy-nx; i++){
        //J[i*nxy+i+nx] = J_c/J_a;
        gsl_matrix_set(mat, i, i+nx, J_c/J_a);
    }
    for(int i = 0; i < nxy-nx; i++){
        //J[(i+nx)*nxy+i] = J_c/J_a;
        gsl_matrix_set(mat, i+nx, i, J_c/J_a);
    }
    for(int i = 1; i < ny; i++){
        //J[(i*nx)*nxy+i*nx-1] = 0.0;
        //J[(i*nx-1)*nxy+i*nx] = 0.0;
        gsl_matrix_set(mat, i*nx, i*nx-1, 0.0);
        gsl_matrix_set(mat, i*nx-1, i*nx, 0.0);
    }
    double sum;
    for(int i = 0; i < nxy; i++){
        //J[i*nxy+i] = sum;
        gsl_vector_view row = gsl_matrix_row(mat, i);
        sum = 0.0;
        for(int j = 0; j < nxy; j++){
            sum += row.vector.data[j];
        }
        gsl_matrix_set(mat, i, i, -sum);
    }
    gsl_mat_inv(mat, &J_inv, nxy);
    //gsl_linalg_LU_decomp(mat, perm, &signum);
    //gsl_linalg_LU_invert(mat, perm, inv);
    // for(int i = 0; i < nxy; i++){
    //     gsl_vector_view row = gsl_matrix_row(mat, i);
    //     for(int j = 0; j < nxy; j++){
    //         //J_inv[i*nxy+j] = row.vector.data[j];
    //         fprintf(stdout,"%4.2g\t",row.vector.data[j]);
    //     }
    //     fprintf(stdout,"\n");
    // }


    gsl_matrix_free(mat);
    //gsl_matrix_free(inv);
    //gsl_permutation_free(perm);
    #endif

    #endif

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    fprintf(stdout, "Mat inv finished at %g seconds.\n",cpu_time_used);

    //calculate velocity field using Navier-Stokes equations
    for(int j = 1; j < ny; j++){
        for(int i = 1; i < nx; i++){
            double a1 = -(u[j*(nx+1)+i+1]*u[j*(nx+1)+i+1] - u[j*(nx+1)+i-1]*u[j*(nx+1)+i-1])/(2*dx) - (u[(j+1)*(nx+1)+i]*(v[(j+1)*(nx+2)+i+1]+v[(j+1)*(nx+2)+i]) - u[j*(nx+1)+i]*(v[j*(nx+2)+i+1]+v[j*(nx+2)+i]))/(4*dy);
            double a3 = (u[j*(nx+1)+i+1]-2*u[j*(nx+1)+i]+u[j*(nx+1)+i-1])/(dx*dx);
            double a4 = (u[(j+1)*(nx+1)+i]-2*u[j*(nx+1)+i]+u[(j-1)*(nx+1)+i])/(dy*dy);

            double A = a1+(a3+a4)/ratio;

            u1[j*(nx+1)+1] = u[j*(nx+1)+1] + dt*(A-(p[j*(nx)+i]-p[j*(nx)+i-1])/(dx));
        }
    }

    for(int j = 1; j < ny; j++){
        for(int i = 1; i <= nx; i++){
            double b1 = -(v[(j+1)*(nx+2)+i]*v[(j+1)*(nx+2)+i] - v[(j-1)*(nx+2)+i]*v[(j-1)*(nx+2)+i])/(2*dx) - ((v[j*(nx+2)+i+1]*(u[(j-1)*(nx+1)+i]+u[j*(nx+1)+i]))-(v[j*(nx+2)+i-1]*(u[(j-1)*(nx+1)+i-1]+u[j*(nx+1)+i-1])))/(4*dx);
            double b3 = (v[(j+1)*(nx+2)+i]-2*v[j*(nx+2)+i]+v[(j-1)*(nx+2)+i])/(dy*dy);
            double b4 = (v[j*(nx+2)+i+1]-2*v[j*(nx+2)+i]+v[j*(nx+2)+i-1])/(dx*dx);

            double B = b1+(b3+b4)/ratio;

            v1[j*(nx+2)+i] = v1[j*(nx+2)+i] + dt*(B-(p[j*(nx)+i-1]-p[(j-1)*(nx)+i-1])/dy);
        }
    }

    // for j=2:Ny
    //     for i=2:Nx
    //
    //         a1=-((u(j,i+1))^2-(u(j,i-1))^2)/(2*dx)  -  (u(j+1,i)*(v(j+1,i+1)+v(j+1,i))-u(j-1,i)*(v(j,i+1)+v(j,i)))/(4*dy);
    //         a3=(u(j,i+1)-2*u(j,i)+u(j,i-1))/(dx^2);            //solving N-S Eq. for u velocity
    //         a4=(u(j+1,i)-2*u(j,i)+u(j-1,i))/(dy^2);
    //
    //         A=a1+(a3+a4)/ratio;
    //
    //         u1(j,i)=u(j,i)+dt*(A-(p(j,i)-p(j,i-1))/dx);
    //
    //     end
    // end
    // for j=2:Ny
    //     for i=2:Nx+1
    //
    //         b1=-((((v(j+1,i))^2-(v(j-1,i))^2)/(2*dy))  -  ((v(j,i+1)*(u(j-1,i)+u(j,i)))-(v(j,i-1)*(u(j-1,i-1)+u(j,i-1))))/(4*dx));
    //         b3=(v(j+1,i)-2*v(j,i)+v(j-1,i))/(dy^2);            //solving N-S for v Velocity
    //         b4=(v(j,i+1)-2*v(j,i)+v(j,i-1))/(dx^2);
    //
    //         B=b1+(b3+b4)/ratio;
    //
    //         v1(j,i)=v(j,i)+dt*(B-(p(j,i-1)-p(j-1,i-1))/dy);
    //
    //     end
    // end





    //apply boundary conditions
    // flow around square
    //u1(Nx/2-Nx/5+2:Nx/2+Nx/5,Ny/2-Ny/5+2:Ny/2+Ny/5) = 0.0;
    //v1(Nx/2-Nx/5+2:Nx/2+Nx/5,Ny/2-Ny/5+2:Ny/2+Ny/5) = 0.0;
    for(int j = nx/2-nx/5+1; j < nx/2+nx/5; j++){
        for(int i = ny/2-ny/5+1; i < ny/2+ny/5; i++){
            u1[j*(nx+1)+i] = 0.0;
            v1[j*(nx+2)+i] = 0.0;
        }
    }


    int iteration;
    //iterate for pressure and velocity corrections
    //for iteration=1:total_iterations           % Iteration Loop
    for(iteration = 0; iteration < total_iterations; iteration++){

        for(int j = 0; j < ny; j++){
            for(int i = 0; i < nx; i++){
                // index order here?
                residual1[j*nx+i] = (u1[j*(nx+1)+i+1] - u1[j*(nx+1)+i] + v1[(j+1)*(nx+2)+i] - v1[j*(nx+2)+i])/(-J_a*dt);
                //fprintf(stdout,"%g\n",residual1[j*nx+i]);
                if(residual1[j*nx+i] > residual_max[iteration]){
                    residual_max[iteration] = residual1[j*nx+i];
                }
            }
        }
        // for j=1:Ny
        //     for i=1:Nx
        //
        //         residual1(j,i)=(u1(j,i+1)-u1(j,i)+v1(j+1,i)-v1(j,i))/(-J_a*dt);    %calculate residuals from continuity
        //
        //     end
        // end

        //skipped this because residual1 is already a vector
        // for j=1:Ny
        //     for i=1:Nx
        //
        //         residual(Nx*(j-1)+i,1)=residual1(j,i);                          %converting residual from a matrix to a vector
        //
        //     end
        // end

        //this is non trivial
        //dp=J\residual;                                              %changes in pressure field
        for(int i = 0; i < nxy; i++){
            sum = 0.0;
            for(int j = 0; j < nxy; j++){
                sum += J_inv[i*nxy+j] * residual1[j];
                // if(residual1[j] > 0.1 || residual1[j] < -0.1){
                //     fprintf(stdout,"%d %d %g %g %g %g\n",i,j,J_inv[i*nxy+j],residual[j],J_inv[i*nxy+j] * residual[j],sum);
                // }
            }
            dp[i] = sum;
            //fprintf(stdout,"%d %g\n",i,sum);
        }

        //skipped for the same reason as above
        // for j=1:Ny
        //     for i=1:Nx
        //         dp1(j,i)=dp(Nx*(j-1)+i,1);                             %converitng changes in pressure field from a vector to a matrix
        //     end
        // end

        for(int j = 1; j < ny; j++){
            for(int i = 1; i < nx; i++){
                u1[j*(nx+1)+i] = u1[j*(nx+1)+i] + relaxation_factor*(dp[j*nx+i-1]-dp[j*nx+i])*dt/dx;
            }
        }
        // for j=2:Ny
        //     for i=2:Nx
        //         u1(j,i)=u1(j,i)+relaxation_factor*(dp1(j,i-1)-dp1(j,i))*dt/dx;                 %u velocity correction
        //     end
        // end

        for(int j = 1; j < ny; j++){
            for(int i = 1; i <= nx; i++){
                v1[j*(nx+2)+i] = v1[j*(nx+2)+i] + relaxation_factor*(dp[(j-1)*nx+i-1]-dp[j*nx+i-1])*dt/dy;
            }
        }
        // for j=2:Ny
        //     for i=2:Nx+1
        //         v1(j,i)=v1(j,i)+relaxation_factor*(dp1(j-1,i-1)-dp1(j,i-1))*dt/dy;             %v velocity correction
        //     end
        // end

        for(int j = 0; j < ny; j++){
            for(int i = 0; i < nx; i++){
                p[j*nx+i] = p[j*nx+i] + relaxation_factor * dp[j*nx+i];
            }
        }
        //p = p + relaxation_factor*dp1;                                      %pressure field correction

        //are these even used?
        //u = u1;
        //v = v1;

        //handled in residual loop
        //residual_max(iteration) = max(abs(residual));                  %output maximum value of residual

        if(residual_max[iteration] < 0.0001){
            break;
        }
    }//iteration ends
    fprintf(stdout,"Exited after %d iterations\n",iteration);
    //     if residual_max(iteration) < 1.0e-4                            %stop on convergance
    //         break
    //     end
    // end                  %iterations ends



    /*
    figure(1)
    contourf (p)  %plot pressure field
    hold on
    UU = u(2:Ny-1,3:Nx);  %select u velocity field (adjust for staggered grid)
    VV = v(2:Ny-1,3:Nx);  %select v velocity field (adjust for staggered grid)
    [X,Y]=meshgrid(2:1:Nx-1,2:1:Ny-1);   %vector plot
    q=quiver(X,Y,UU,VV,1);
    q.Color = 'black';
    axis equal;
    %draw a square
    v_ = [Nx/2-Nx/5+1 Ny/2-Ny/5+1.5;  Nx/2-Nx/5+1 Ny/2+Ny/5+0.5; Nx/2+Nx/5 Ny/2-Ny/5+1.5; Nx/2+Nx/5 Ny/2+Ny/5+.5];
    f_ = [2 1 3 4];
    p_=patch('Faces',f_,'Vertices',v_,'FaceColor','red');
    p_.EdgeColor='none';
    p_.FaceColor='white';
    xlabel ('x-dimension')
    ylabel ('y-dimension')
    title ('Pressure and velocity field around the square')
    colorbar
    figure()
    plot(residual_max)  %residual plot
    xlabel ('Iteration number')
    ylabel ('Maximum value of residual')
    title ('Convergence plot')
    */
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

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    fprintf(stdout, "Program ended at %g seconds.\n",cpu_time_used);
}
