#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
//#include <mpi.h>

double findMaxMat(double *Mat, int x, int y) {
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

void setBoundary(double *Matu, double *Matv, int nx, int ny) {
	int i,j;
	for(i = (nx/2 - nx/5); j < (nx/2 + nx/5)+2; i++) {
		for(j = (ny/2 - ny/5); j < (ny/2 - ny/5 + 2); j++){
			Matu[i][j] = 0.0;
			Matv[i][j] = 0.0;
		}
	}
}

void matCopy(double MatA, double MatB, int nx, int ny) {
	int i,j;
	
	for(int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			MatB[i][j] = MatA[i][j];
		}
	}
}

int main ( int argc, char **argv ){
    //space variables (ENTER)
    int nx = 50;           //number of columns
    int ny = 50;           //number of rows
    //nx = 5;ny = 5;
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

    //size and dimension of pressure and velocity (AUTO)
    double *p = calloc(nxy,sizeof(double));
    double *u = calloc((nx+1)*(ny+1),sizeof(double));
    for(int i = 0; i < (nx+1)*(ny+1); i++){
        u[i] = 0.1;
    }
    double *v = calloc((nx+2)*(ny+1),sizeof(double));
    double *residual = calloc(nxy,sizeof(double));
    double *dp = calloc(nxy,sizeof(double));

    //temporary variables
    double *u1 = calloc((nx+1)*(ny+1),sizeof(double));
    for(int i = 0; i < (nx+1)*(ny+1); i++){
        u1[i] = 0.1;
    }
    double *v1 = calloc((nx+2)*(ny+1),sizeof(double));
    double *dp1 = calloc(nxy,sizeof(double));
    double *residual1 = calloc(nxy,sizeof(double));

    //timestep value, relaxation factor, number of iterations (ENTER)
    double dt = 1;
    double relaxation_factor = 0.5;
    int total_iterations = 200;
    double *residual_max = calloc(total_iterations,sizeof(double));

    //check CFL criteria
	double CFL_x = findMax(u, nx+1, ny+1) * dt/dx;
	double CFL_y = findMax(u, nx+1, ny+1)*dt/dy;

    //calculate sparse matrix (AUTO)
    double J_a = 2*(dxi*dxi+dyi*dyi);
    double J_b = -dyi;
    double J_c = -dxi;
    double *J = malloc(nxy*nxy*sizeof(double)); //check sizing here
    for(int i = 0; i < nxy*5; i++){
        J_row[i] = -1;
        J_col[i] = -1;
    }
    //J=spalloc(Nx*Ny,Nx*Ny,(Nx-2)*(Ny-2)*4+Nx*Ny);
	
	int i;
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
        J[i*nx+1][i*nx]=0;
		J[i*nx][i*nx+1]=0;
    }

	int j;
	int sum;
    for(i = 0; i < nxy; i++){
		sum = 0;
      for(j = 0; j < nxy; j++) {
		  sum += J[i][j];
	  }
	  J[i][i] = -sum;
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
	setBoundary(u1,v1,nx,ny);

    //iterate for pressure and velocity corrections
	int iteration;
    for(iteration = 0; iteration < total_iterations; iteration++){

        for(j = 0; j < ny; j++){
            for(i = 0; i < nx; i++){
                residual1[j][i]=(u1[j][i+1]-u1[j][i]+v1[j+1][i]-v1[j][i])/(-J_a*dt);
            }
        }

        //Calculate changes in pressure field
        //dp=J\residual;                                              
		
        for(j = 1; j < ny; j++){
            for(i = 1; i < nx; i++){
                u1[j][i]=u1[j][i]+relaxation_factor*(dp[j][i-1]-dp1[j][i])*dt*dxi;
            }
        }

        for(j = 1; j < ny; j++){
            for(int i = 1; i <= nx; i++){
                v1[j][i]=v1[j][i]+relaxation_factor*(dp1[j-1][i-1]-dp1[j][i-1])*dt*dyi;
            }
        }

        for(j = 0; j < ny; j++){
            for(i = 0; i < nx; i++){
                p[j*nx+i] = p[j*nx+i] + relaxation_factor * dp1[j*nx+i];
            }
        }
		
		
        //p = p + relaxation_factor*dp1;                                    

        //are these even used?
        matCopy(u1, u, nx, ny);
		matCopy(v1, v, nx, ny);

        residual_max[iteration] = findMaxMat(residuall, nx, ny);                  %output maximum value of residual

        if(residual_max[iteration] < 0.0001){
            break;
        }
    }//iteration ends

	//free variables
	
	//export data to file
	
	return 0;
}
