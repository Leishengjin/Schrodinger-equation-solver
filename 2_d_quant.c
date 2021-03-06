//Time evolution of the dynamics of a 2-D quantum particle in a harmonic oscillator.
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>

#define N 128 //Size of the spatial grid (Note that N variable is reserved, so this is later saved into another variable n.).

static double pi;
static double buf[N][2];   // for fft
static double buf2[N][2];  // for fft_2D

int main(nfile,filenames)
int nfile; char *filenames[];
{
  int i,j,k,n,it;
  double delta_t,xk,ccr,cci; //Initialise the vectors that are going to be used.
  double cfou[N+1][2],psi[N][N][2],pot,cpot[N][N][2],cmom[N][2]; //Initialise the matrices that are to be used.
  double xx1,xx2,xmag,xnorm,cmomr,cmomi;
  void fft(),fft_2D();


  FILE * temp = fopen("data.temp", "w"); //Initialise the file that the states will be recorded into.

  n=N; //number of space points
  pi=4.*atan((double)1.0); // define pi

  delta_t = 0.1; //size of the time step

//Fourier coefficients (?)
  for(i=0;i<n+1;i++){
    cfou[i][0]=cos(2.*pi*i/n);
    cfou[i][1]=sin(2.*pi*i/n);
  }

//Define the initial state and the potential; also the time evolution operators
//Loop is for all space points (0:n)
  for(i=0;i<n;i++){ for(j=0;j<n;j++){


    psi[i][j][0]=exp(-(i-0.75*n)*(i-0.75*n)/36. -(j-n/2)*(j-n/2)/36. ) +
                 exp(-(i-n/2)*(i-n/2)/36.       -(j-0.75*n)*(j-0.75*n)/36. ); //Initial state is two Gaussians with exchange symmetry
    psi[i][j][1]=0.; //Initial state has no imaginary part (0 momentum)

    pot = 0.001*((i-n/2)*(i-n/2)+(j-n/2)*(j-n/2)); //Potential is harmonic (quadratic)
    if(abs(i-j)<2){pot += 0.1;} // delta function repulsion

    cpot[i][j][0] =  cos(delta_t*pot)/(n*n); //Real part of the time evolution operator due to the potential
    cpot[i][j][1] = -sin(delta_t*pot)/(n*n); //Imaginarty part of the time evolution operator due to the potential

    } // close j-loop

    xk = 2.*sin(i*pi/n); //Derivative operator (momentum) in k-space
    cmom[i][0] =  cos(delta_t*xk*xk); //Real part of the time evolution operator due to the kinetic energy
    cmom[i][1] = -sin(delta_t*xk*xk); //Imaginary part of the time evolution operator due to the kinetic energy
  } // close i-loop

//Evolve the initial state using the time evolution operator defined above.
//Loop is for all time points (1000 points since 0:0.1:100).
  for(it=0;it<2000;it++){

    for(i=0;i<n;i++) for(j=0;j<n;j++){
      ccr = psi[i][j][0]*cpot[i][j][0] - psi[i][j][1]*cpot[i][j][1]; //First, act on the state by the potential part
                                                                     // of the evolution operator (real part)
      cci = psi[i][j][0]*cpot[i][j][1] + psi[i][j][1]*cpot[i][j][0]; //(imaginary part)
      psi[i][j][0] = ccr; //update the state (real)
      psi[i][j][1] = cci; //(imaginary)
    }

    fft_2D(n,cfou,psi,1); //Evaluate the Fourier transform to pass into momentum space and
                         // to be able to use the linear kinetic evolution operator.

    for(i=0;i<n;i++) for(j=0;j<n;j++){
      cmomr = cmom[i][0]*cmom[j][0] - cmom[i][1]*cmom[j][1]; // real part of exp(ikx^2)*exp(iky^2)
      cmomi = cmom[i][0]*cmom[j][1] + cmom[i][1]*cmom[j][0]; // imag part of exp(ikx^2)*exp(iky^2)
      ccr = psi[i][j][0]*cmomr - psi[i][j][1]*cmomi; //Act on the state by the kinetic part (real)
      cci = psi[i][j][0]*cmomi + psi[i][j][1]*cmomr; //(imaginary)
      psi[i][j][0] = ccr; //update real part
      psi[i][j][1] = cci; //and the imaginary part
    }

    fft_2D(n,cfou,psi,-1); //Evaluate the inverse Fourier transform to go back into position space.

    if(it%1 == 0){

// calculate xx1 = <x1>  and xx2 = <x2> given x1>x2
       xx1 = 0.; //initialise
       xx2 = 0.;
       xnorm = 0.;
       for(i=0;i<n;i++) for(j<=i;j<n;j++){
         xmag = psi[i][j][0]*psi[i][j][0] + psi[i][j][1]*psi[i][j][1]; 
         xx1 += i*xmag; //(+= means sum)
         xx2 += j*xmag; 
         xnorm += xmag; //sum all values in the |psi| matrix to get the norm
       }
       xx1 = xx1/xnorm; //normalise
       xx2 = xx2/xnorm;
       fprintf(temp,"%d %lf %lf\n",it,xx1,xx2);
//     fflush(temp);
     } // closes if

  } // closes it loop

}
