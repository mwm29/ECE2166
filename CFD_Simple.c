#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
//#include <mpi.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

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

double findMaxMat(double **Mat, int x, int y) {
	double max = Mat[0][0];
	for(int i = 0; i < x; i++) {
		for(int j = 0; j < y; j++) {
			if(Mat[i][j] > max) {
				max = Mat[i][j];
			}
		}
	}
	return max;
}

void setBoundary(double **Matu, double **Matv, int nx, int ny) {
	int i,j;
	for(i = (nx/2 - nx/5 + 2)-1; i < (nx/2 + nx/5); i++) {
		for(j = (ny/2 - ny/5 + 2)-1; j < (ny/2 + ny/5); j++){
			Matu[i][j] = 0.0;
			Matv[i][j] = 0.0;
		}
	}
}

void matCopy(double **MatA, double **MatB, int nx, int ny) {
	int i,j;
	for(int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			MatB[i][j] = MatA[i][j];
		}
	}
}

void  matMultVec(double **Mat, double *Vec, int dim, double *res) {
	int i, j;
	for(i = 0; i < dim; i++) {
		for(j = 0; j < dim; j++) {
			res[i] += Mat[i][j] * Vec[j];
		}
	}
}

int main ( int argc, char **argv ){
	//Loop variables
	int i,j;
	
    //space variables
    int nx = 10;           //number of columns
    int ny = 10;           //number of rows
    int nxy=nx*ny;
	
    //discretization variables
    double dx = 1;       //x-grid size
    double dy = 1;       //y-grid size
	double dxi = 1/dx;
	double dyi = 1/dy;
	
    //fluid density & viscosity
    double density = 1000;
    double viscosity = 1;
    double ratio = density/viscosity;       //ratio of density and dynamic viscosity
	double invRatio = 1/ratio;
	
    //size and dimension of pressure and velocity
    double **p = (double **)malloc(nx*sizeof(double *));
	double **u = (double **)malloc((nx+1)*sizeof(double *));
	double **v = (double **)malloc((nx+1)*sizeof(double *));
	double *residual = (double *)calloc(nxy,sizeof(double));
    double *dp = (double *)calloc(nxy,sizeof(double));
	for(i = 0; i < nx; i++) {
		p[i] = (double *)calloc(ny, sizeof(double));
	}
	for(i = 0; i < (nx+1); i++) {
		u[i] = (double *)calloc((ny+1),sizeof(double));
		v[i] = (double *)calloc((ny+2),sizeof(double));
	}
    for(i = 0; i < (nx+1); i++){
		for(j = 0; j < (ny+1); j++) {
			u[i][j] = 0.1;
		}
    }
	
    //temporary variables
    double **u1 = (double **)malloc((nx+1)*sizeof(double *));
	double **v1 = (double **)malloc((nx+1)*sizeof(double *));
	double **dp1 = (double **)malloc(nx*sizeof(double *));
	double **residual1 = (double **)malloc(nx*sizeof(double *));
	for(i = 0; i < (nx+1); i++) {
		u1[i] = (double *)calloc((ny+1),sizeof(double));
		v1[i] = (double *)calloc((ny+2),sizeof(double));
	}
	matCopy(u, u1, nx+1, ny+1);
	for(i = 0; i < nx; i++) {
		dp1[i] = (double *)calloc(ny,sizeof(double));
		residual1[i] = (double *)calloc(ny,sizeof(double));
	}	

    //timestep value, relaxation factor, number of iterations (ENTER)
    double dt = 1;
    double relaxation_factor = 0.5;
    int total_iterations = 200;
    double *residual_max = (double *)calloc(total_iterations,sizeof(double));

    //check CFL criteria
	double CFL_x = findMaxMat(u, nx+1, ny+1)*dt*dxi;
	double CFL_y = findMaxMat(u, nx+1, ny+1)*dt*dyi;

    //calculate sparse matrix (AUTO)
    double J_a = 2*(dxi*dxi+dyi*dyi);
    double J_b = -dyi*dyi;
    double J_c = -dxi*dxi;
	
    double **J = (double **)malloc(nxy*sizeof(double *));
	double **invJ = (double **)malloc(nxy*sizeof(double *));
	for(i = 0; i < nxy; i++) {
		J[i] = (double *)calloc(nxy,sizeof(double));
		invJ[i] = (double *)calloc(nxy,sizeof(double));
	}
	
    for(i = 0; i < nxy-1; i++){
        J[i][i+1] = J_b/J_a;
    }

    for(i = 1; i < nxy; i++){
        J[i][i-1] = J_b/J_a;
    }

    for(i = 0; i < nxy-nx; i++){
        J[i][i+nx] = J_c/J_a;
    }

    for(i = 0; i < nxy-nx; i++){
        J[i+nx][i] = J_c/J_a;
    }

    for(i = 0; i < ny-1; i++){
        J[(i+1)*nx][(i+1)*nx-1]=0;
		J[(i+1)*nx-1][(i+1)*nx]=0;
    }

    for(i = 0; i < nxy; i++){
		double sum = 0;
		for(j = 0; j < nxy; j++) {
			sum += J[i][j];
		}
		J[i][i] = -sum;
    }
	gsl_matrix *mat = gsl_matrix_calloc(nxy,nxy);
	double *J_inv = calloc(nxy*nxy,sizeof(double));
	for(i = 0; i < nxy; i++) {
		for(j = 0; j < nxy; j++) {
			gsl_matrix_set(mat,i,j,J[i][j]);
		}
	}
	gsl_mat_inv(mat, &J_inv, nxy);
	for(i = 0; i < nxy; i++) {
		for(j = 0; j < nxy; j++) {
			invJ[i][j] = J_inv[i*nxy+j];
		}
	}
	
    //calculate velocity field using Navier-Stokes equations
    for(j = 1; j < ny; j++){
        for(i = 1; i < nx; i++){
			double a1=-((u[j][i+1])*(u[j][i+1])-(u[j][i-1])*(u[j][i-1]))*(0.5*dxi)-(u[j+1][i]*(v[j+1][i+1]+v[j+1][i])-u[j-1][i]*(v[j][i+1]+v[j][i]))*(0.25*dyi);
			double a3=(u[j][i+1]-2*u[j][i]+u[j][i-1])*(dxi*dxi);    
			double a4=(u[j+1][i]-2*u[j][i]+u[j-1][i])*(dyi*dyi);
			double A=a1+(a3+a4)*invRatio;
			u1[j][i]=u[j][i]+dt*(A-(p[j][i]-p[j][i-1])*dxi);
        }
    }
	
	
    for(j = 1; j < ny; j++){
        for(i = 1; i <= nx; i++){
			double b1=-((((v[j+1][i])*(v[j+1][i])-(v[j-1][i])*(v[j-1][i]))*(0.5*dyi))-((v[j][i+1]*(u[j-1][i]+u[j][i]))-(v[j][i-1]*(u[j-1][i-1]+u[j][i-1])))*(0.25*dxi));
			double b3=(v[j+1][i]-2*v[j][i]+v[j-1][i])*(dyi*dyi);
			double b4=(v[j][i+1]-2*v[j][i]+v[j][i-1])*(dxi*dxi);
			double B=b1+(b3+b4)*invRatio;
			v1[j][i]=v[j][i]+dt*(B-(p[j][i-1]-p[j-1][i-1])*dyi);
        }
    }
	
    //apply boundary conditions
    // flow around square
	setBoundary(u1,v1,nx+1,ny+1);
	
    //iterate for pressure and velocity corrections
	int iteration;
    for(iteration = 0; iteration < total_iterations; iteration++){

        for(j = 0; j < ny; j++){
            for(i = 0; i < nx; i++){
                residual1[j][i]=(u1[j][i+1]-u1[j][i]+v1[j+1][i]-v1[j][i])/(-J_a*dt);
            }
        }
		
		for(j = 0; j < ny; j++) {
			for(i = 0; i < nx; i++) {
				residual[nx*j+i] = residual1[j][i];
			}
		}
		
        //Calculate changes in pressure field
        matMultVec(invJ, residual, nxy, dp);
		
		for(i = 0; i < nx; i++) {
			for(j = 0; j < ny; j++) {
				dp1[i][j] = dp[nx*i+j];
			}
		}

        for(j = 1; j < ny; j++){
            for(i = 1; i < nx; i++){
                u1[j][i]=u1[j][i]+relaxation_factor*(dp1[j][i-1]-dp1[j][i])*dt*dxi;
            }
        }

        for(j = 1; j < ny; j++){
            for(i = 1; i < nx; i++){
                v1[j][i]=v1[j][i]+relaxation_factor*(dp1[j-1][i-1]-dp1[j][i-1])*dt*dyi;
            }
        }

        for(j = 0; j < ny; j++){
            for(i = 0; i < nx; i++){
                p[j][i] = p[j][i] + relaxation_factor * dp1[j][i];
            }
        }                                   

        //are these even used?
        matCopy(u1, u, nx, ny);
		matCopy(v1, v, nx, ny);

        residual_max[iteration] = findMaxMat(residual1, nx, ny);

        if(residual_max[iteration] < 0.0001){
            break;
        }
    }//iteration ends
	
	for(i = 0; i < nx+1; i++) {
		for(j = 0; j < ny+1; j++) {
			printf("%lf ",u[i][j]);
		}
		printf("\n");
	}
	
	//free variables
	for(i = 0; i < nx; i++) {
		free(p[i]);
		free(dp1[i]);
		free(residual1[i]);
	}
	for(i = 0; i < (nx+1); i++) {
		free(u[i]);
		free(v[i]);
		free(u1[i]);
		free(v1[i]);
	}
	for(i = 0; i < nxy; i++) {
		free(J[i]);
	}
	free(p);
	free(u);
	free(v);
	free(residual);
    free(dp);
    free(u1); 
	free(v1);
	free(dp1);
	free(residual1);
	free(residual_max);
	free(J);
	
	return 0;
}
