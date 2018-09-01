//Time evolution of the dynamics of 3 quantum particles in a harmonic oscillator.
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>

#define N 128 //Size of the spatial grid (Note that N variable is reserved, so this is later saved into another variable n.).

static double pi;
static double buf[N][2];   // for fft
static double buf3[N][2];  // for fft_3D 
static double psi[N][N][N][2],cpot[N][N][N][2]; 

int main(nfile,filenames)
int nfile; char *filenames[];
{
  int i,j,k,n,it;
  double delta_t,xk,ccr,cci,n01,n02,n03; //Initialise the vectors that are going to be used.
  double cfou[N+1][2],pot,cmom[N][2]; //Initialise the matrices that are to be used.
  double xx1,xx2, xx3,xmag,xnorm,cmomr,cmomi;
  void fft(),fft_3D();


  FILE * temp = fopen("data.temp", "w"); //Initialise the file that the states will be recorded into.

  n=N; //number of space points
  pi=4.*atan((double)1.0); // define pi

  delta_t = 0.01; //size of the time step

  //Fourier coefficients (?)
  for(i=0;i<n+1;i++){
    cfou[i][0]=cos(2.*pi*i/n);
    cfou[i][1]=sin(2.*pi*i/n);
  }
  
  // Define the positions of the three particles: "a displaced particle and two particles near the origin"
  n01=0.75*n;
  n02=0.5*n +5;
  n03=0.5*n -5;

  //Define the initial state and the potential; also the time evolution operators
  //Loop is for all space points (0:n)
  for(i=0;i<n;i++){ for(j=0;j<n;j++){ for(k=0;k<n;k++){

      //Initial state is Gaussians on the axes with exhange symmetry
      psi[i][j][k][0]=exp(-(i-n01)*(i-n01)/36. -(j-n02)*(j-n02)/36. -(k-n03)*(k-n03)/36 ) +
                      exp(-(i-n01)*(i-n01)/36. -(k-n02)*(k-n02)/36. -(j-n03)*(j-n03)/36 ) +
                      exp(-(j-n01)*(j-n01)/36. -(i-n02)*(i-n02)/36. -(k-n03)*(k-n03)/36 ) +
                      exp(-(j-n01)*(j-n01)/36. -(k-n02)*(k-n02)/36. -(i-n03)*(i-n03)/36 ) +
                      exp(-(k-n01)*(k-n01)/36. -(i-n02)*(i-n02)/36. -(j-n03)*(j-n03)/36 ) +
                      exp(-(k-n01)*(k-n01)/36. -(j-n02)*(j-n02)/36. -(i-n03)*(i-n03)/36 );
      psi[i][j][k][1]=0.; //Initial state has no imaginary part (or 0 momentum)

      pot = 0.001*((i-n/2)*(i-n/2)+(j-n/2)*(j-n/2)+(k-n/2)*(k-n/2)); //External potential is harmonic (quadratic)
      if(abs(i-j)<2){pot += 1.0;} // delta function repulsion 
      if(abs(i-k)<2){pot += 1.0;}
      if(abs(j-k)<2){pot += 1.0;}

      cpot[i][j][k][0] =  cos(delta_t*pot)/(n*n*n); //Real part of the time evolution operator due to the potential
      cpot[i][j][k][1] = -sin(delta_t*pot)/(n*n*n); //Imaginary part of the time evolution operator due to the potential

     } // close k-loop
    } // close j-loop

    xk = 2.*sin(i*pi/n); //Derivative operator (momentum) in k-space
    cmom[i][0] =  cos(delta_t*xk*xk); //Real part of the time evolution operator due to the kinetic energy
    cmom[i][1] = -sin(delta_t*xk*xk); //Imaginary part of the time evolution operator due to the kinetic energy
  } // close i-loop

//Evolve the initial state using the time evolution operator defined above.
//Loop is for all time points (1000 points since 0:0.1:100).
  for(it=0;it<200000;it++){

    for(i=0;i<n;i++) for(j=0;j<n;j++) for(k=0;k<n;k++){
      ccr = psi[i][j][k][0]*cpot[i][j][k][0] - psi[i][j][k][1]*cpot[i][j][k][1]; //First, act on the state by the potential part
                                                                     // of the evolution operator (real part)
      cci = psi[i][j][k][0]*cpot[i][j][k][1] + psi[i][j][k][1]*cpot[i][j][k][0]; //(imaginary part)
      psi[i][j][k][0] = ccr; //update the state (real)
      psi[i][j][k][1] = cci; //(imaginary)
    }

    fft_3D(n,cfou,psi,1); //Evaluate the Fourier transform to pass into momentum space and
                         // to be able to use the linear kinetic evolution operator.

    for(i=0;i<n;i++) for(j=0;j<n;j++) for(k=0;k<n;k++){
      cmomr = cmom[i][0]*cmom[j][0]*cmom[k][0] - cmom[i][1]*cmom[j][1]*cmom[k][0] - 
              cmom[i][1]*cmom[j][0]*cmom[k][1] - cmom[i][0]*cmom[j][1]*cmom[k][1]; // real part of exp(ikx^2)*exp(iky^2)
      cmomi = cmom[i][0]*cmom[j][1]*cmom[k][0] + cmom[i][1]*cmom[j][0]*cmom[k][0] +
              cmom[i][0]*cmom[j][0]*cmom[k][1] + cmom[i][1]*cmom[j][1]*cmom[k][1]; // imag part of exp(ikx^2)*exp(iky^2)
      ccr = psi[i][j][k][0]*cmomr - psi[i][j][k][1]*cmomi; //Act on the state by the kinetic part (real)
      cci = psi[i][j][k][0]*cmomi + psi[i][j][k][1]*cmomr; //(imaginary)
      psi[i][j][k][0] = ccr; //update real part
      psi[i][j][k][1] = cci; //and the imaginary part
    }

    fft_3D(n,cfou,psi,-1); //Evaluate the inverse Fourier transform to go back into position space.

// %10 less data
     if(it%100 == 0){
    // calculate xx1 = <x1>, xx2 =<x2>  and xx3 = <x3> 
       xx1 = 0.; //initialise
       xx2 = 0.;
       xx3 = 0.;
       xnorm = 0.;
       for(i=0;i<n;i++) for(j=i;j<n;j++) for(k=j; k<n; k++){
         xmag = psi[i][j][k][0]*psi[i][j][k][0] + psi[i][j][k][1]*psi[i][j][k][1]; 
         xx1 += i*xmag; //(+= means sum)
         xx2 += j*xmag; 
         xx3 += k*xmag; 
         xnorm += xmag; //sum all values in the |psi| matrix to get the norm
       }
       xx1 = xx1/xnorm; //normalise
       xx2 = xx2/xnorm;
       xx3 = xx3/xnorm;
       fprintf(temp,"%d %lf %lf %lf\n",it,xx1,xx2,xx3);
       fflush(temp);
     } // closes if

  } // closes it loop

}//closes the program



// Define the 3D Fourier transform
void fft_3D(n,cfou,vect3D,mode) int n; double cfou[][2]; double vect3D[n][n][n][2]; int mode;
{
  int i,j,k;
  void fft();

  for(j=0;j<n;j++) for(k=0;k<n;k++){  // FT the x direction
    for(i=0;i<n;i++){  //carry the i'th x dir.
      buf3[i][0] = vect3D[i][j][k][0];
      buf3[i][1] = vect3D[i][j][k][1];
    }

    // do the FT
    fft(n,cfou,buf3,mode);

    for(i=0;i<n;i++){  //carry the i'th x dir. back
      vect3D[i][j][k][0] = buf3[i][0];
      vect3D[i][j][k][1] = buf3[i][1];
    }
  }

  for(j=0;j<n;j++) for(k=0;k<n;k++){  // FT the y direction
    for(i=0;i<n;i++){  //carry the i'th y dir. 
      buf3[i][0] = vect3D[j][i][k][0];
      buf3[i][1] = vect3D[j][i][k][1];
    }

    // do the FT
    fft(n,cfou,buf3,mode);

    for(i=0;i<n;i++){  //carry the i'th y dir. back
      vect3D[j][i][k][0] = buf3[i][0];
      vect3D[j][i][k][1] = buf3[i][1];
    }
  }
 
  for(j=0;j<n;j++) for(k=0;k<n;k++){  // FT the z direction
    for(i=0;i<n;i++){  //carry the i'th z dir.
      buf3[i][0] = vect3D[k][j][i][0];
      buf3[i][1] = vect3D[k][j][i][1];
    }

    // do the FT
    fft(n,cfou,buf3,mode);

    for(i=0;i<n;i++){  //carry the i'th z dir. back
      vect3D[k][j][i][0] = buf3[i][0];
      vect3D[k][j][i][1] = buf3[i][1];
    }
  }
}//closes the fft_3D program


//The fast Fourier transform algorithm

void fft(n,cfou,vect,mode) int n; double cfou[][2]; double vect[][2]; int mode;
{
  int i,j,k,nn,J,J0,K,L,L0,M,length,index,index_0,index_incr;
  int no2;
  double cfr,cfi,buf1r,buf1i,buf2r,buf2i;

  for(i=0;i<n;i++){
    for(j=0;j<2;j++){
      buf[i][j]=vect[i][j];
    }
  }

  if(mode == 1){
    index_0=0;
    index_incr=n;
  }
  else{
    index_0=n;
    index_incr= -n;
  }
  length=1;
  no2=n/2;

  for(nn=no2;nn>0;nn=nn/2){  //   nn = no of ft to be computed at this step
    index_incr/=2;
    index=index_0;
    J0=0;
    L0=0;
    for(k=0;k<length;k++){   //   length = length of the ft at this step
      cfr=cfou[index][0];
      cfi=cfou[index][1];
      for(j=0;j<nn;j++){     //   go over all the ft's at this step
        J=j+J0;
        K=J+no2;
        L=j+L0;
        M=L+nn;

        buf1r=buf[L][0];
        buf1i=buf[L][1];
        buf2r=buf[M][0];
        buf2i=buf[M][1];
        vect[J][0]=buf1r + cfr*buf2r - cfi*buf2i;
        vect[J][1]=buf1i + cfi*buf2r + cfr*buf2i;
        vect[K][0]=buf1r - cfr*buf2r + cfi*buf2i;
        vect[K][1]=buf1i - cfi*buf2r - cfr*buf2i;
      }   // for - j

      index+=index_incr;
      J0=J0+nn;
      L0=L0+nn+nn;
    }    //  for - k

    length=length+length;

    if(nn != 1){
      for(i=0;i<n;i++){
        for(j=0;j<2;j++){
          buf[i][j]=vect[i][j];
        }
      }
    }

  }   //  for - nn
  return;
}//closes the fft program