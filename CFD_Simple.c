#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
//#include <mpi.h>

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


int main ( int argc, char **argv ){
    //space variables (ENTER)
    int nx = 50;           //number of columns
    int ny = 50;           //number of rows
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
    //residual_max = zeros(total_iterations,1);
    double *residual_max = calloc(total_iterations,sizeof(double));

    //check CFL criteria (CHECK!)
    //CFL_x = max(max(u))*dt/dx;
    //CFL_y = max(max(u))*dt/dy;

    //calculate sparse matrix (AUTO)
    double J_a = 2*(1/(dx*dx)+1/(dy*dy));
    double J_b = -1/(dy*dy);
    double J_c = -1/(dx*dx);
    int *J_row = malloc(((nx-2)*(ny-2)*4+nx*ny)*sizeof(int));
    int *J_col = malloc(((nx-2)*(ny-2)*4+nx*ny)*sizeof(int));
    double *J = malloc(((nx-2)*(ny-2)*4+nx*ny)*sizeof(double));
    for(int i = 0; i < ((nx-2)*(ny-2)*4+nx*ny); i++){
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
    for(int i = 0; i < nxy; i++){
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
    for(int i = 0; i < ny-1; i++){
        int index = sparse_find(J_row, J_col, J_num, i*nx+1, i*nx);
        // if(index == -1){
        //     fprintf(stderr,"not found1\n");
        // }
        J_num--;
        J_row[index] = J_row[J_num];
        J_col[index] = J_col[J_num];
        J[index] = J[J_num];

        index = sparse_find(J_row, J_col, J_num, i*nx, i*nx+1);
        // if(index == -1){
        //     fprintf(stderr,"not found2\n");
        // }
        J_num--;
        J_row[index] = J_row[J_num];
        J_col[index] = J_col[J_num];
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
        J[J_num] = sparse_row_sum(J_row, J, J_num, i);
        J_num++;
    }
    //%spy(J)   %for checking sparse matrix

    //calculate velocity field using Navier-Stokes equations
    for(int j = 1; j < ny; j++){
        for(int i = 1; i < nx; i++){
            double a1 = -(u[j*(nx+1)+i+1]*u[j*(nx+1)+i+1] - u[j*(nx+1)+i-1]*u[j*(nx+1)+i-1])/(2*dx) - (u[(j+1)*(nx+1)+i]*(v[(j+1)*(nx+2)+i+1]+v[(j+1)*(nx+2)+i]) - u[j*(nx+1)+i]*(v[j*(nx+2)+i+1]+v[j*(nx+2)+i]))/(4*dy);
            double a3 = (u[j*(nx+1)+i+1]-2*u[j*(nx+1)+i]+u[j*(nx+1)+i-1])/(dx*dx);
            double a4 = (u[(j+1)*(nx+1)+i]-2*u[j*(nx+1)+i]+u[(j-1)*(nx+1)+i])/(dy*dy);

            double A = a1+(a3+a4)/ratio;

            u1[j*(nx+1)+1] = u[j*(nx+1)+1] + dt*(A-(p[j*(nx+1)+i]-p[j*(nx+1)+i-1])/(dx));
        }
    }

    for(int j = 1; j < ny; j++){
        for(int i = 1; i <= nx; i++){
            double b1 = -(v[(j+1)*(nx+2)+i]*v[(j+1)*(nx+2)+i] - v[(j-1)*(nx+2)+i]*v[(j-1)*(nx+2)+i])/(2*dx) - ((v[j*(nx+2)+i+1]*(u[(j-1)*(nx+1)+i]+u[j*(nx+1)+i]))-(v[j*(nx+2)+i-1]*(u[(j-1)*(nx+1)+i-1]+u[j*(nx+1)+i-1])))/(4*dx);
            double b3 = (v[(j+1)*(nx+2)+i]-2*v[j*(nx+2)+i]+v[(j-1)*(nx+2)+i])/(dy*dy);
            double b4 = (v[j*(nx+2)+i+1]-2*v[j*(nx+2)+i]+v[j*(nx+2)+i-1])/(dx*dx);

            double B = b1+(b3+b4)/ratio;

            v1[j*(nx+1)+i] = v1[j*(nx+1)+i] + dt*(B-(p[j*(nx+1)+i-1]-p[(j-1)*(nx+1)+i-1])/dy);
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




    //iterate for pressure and velocity corrections
    //for iteration=1:total_iterations           % Iteration Loop
    for(int iteration = 0; iteration < total_iterations; iteration++){

        for(int j = 0; j < ny; j++){
            for(int i = 0; i < nx; i++){
                // index order here?
                residual1[j*nx+i] = (u1[j*(nx+1)+i+1] - u1[j*(nx+1)+i] + v1[(j+1)*(nx+1)+i] - v1[j*(nx+1)+i])/(-J_a*dt);
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
        dp=J\residual;                                              %changes in pressure field

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
                p[j*nx+i] = p[j*nx+i] + relaxation_factor * dp1[j*nx+i];
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
}
