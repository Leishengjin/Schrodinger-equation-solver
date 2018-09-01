//Time evolution of the dynamics of a 1-D quantum particle in a harmonic oscillator.
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#define N 128 //Size of the spatial grid (Note that N variable is reserved, so this is later saved into another variable n.).

static double pi;
static double buf[N][2];

int main(nfile,filenames)
int nfile; char *filenames[];
{
  int i,j,k,n,it;
  double delta_t,xk,ccr,cci; //Initialise the parameters that are going to be used.
  double cfou[N+1][2],psi[N][2],pot[N],cpot[N][2],cmom[N][2]; //Initialise the matrices that are to be used.
  void fft();

  FILE * temp = fopen("data.temp", "w"); //Initialise the file that the states will be recorded into.
  FILE * gnuplotPipe = popen ("gnuplot -persistent", "w"); //Pipe to GNUplot for plotting.

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
  for(i=0;i<n;i++){
    psi[i][0]=exp(-(i-n/2)*(i-n/2)/40. ); //Real part of the initial state is a Gaussian
    psi[i][1]=0.; //Initial state has no imaginary part
    pot[i] = 0.00005*(i-n/2)*(i-n/2); //Potential is harmonic (quadratic)
    cpot[i][0] =  cos(delta_t*pot[i])/n; //Real part of the time evolution operator due to the potential
    cpot[i][1] = -sin(delta_t*pot[i])/n; //Imaginarty part of the time evolution operator due to the potential
    xk = 2.*sin(i*pi/n); //Derivative operator (momentum) in k-space
    cmom[i][0] =  cos(delta_t*xk*xk); //Real part of the time evolution operator due to the kinetic energy
    cmom[i][1] = -sin(delta_t*xk*xk); //Imaginary part of the time evolution operator due to the kinetic energy
  }

//Evolve the initial state using the time evolution operator defined above.
//Loop is for all time points (1000 points since 0:0.1:100).
  for(it=0;it<500;it++){
    for(i=0;i<n;i++){
      ccr = psi[i][0]*cpot[i][0] - psi[i][1]*cpot[i][1]; //First, act on the state by the potential part of the evolution operator (real part)
      cci = psi[i][0]*cpot[i][1] + psi[i][1]*cpot[i][0]; //(imaginary part)
      psi[i][0] = ccr; //update the state (real)
      psi[i][1] = cci; //(imaginary)
    }

    fft(n,cfou,psi,1); //Take the Fourier transform to pass into momentum space and to be able to use the linear kinetic evolution operator.

    for(i=0;i<n;i++){
      ccr = psi[i][0]*cmom[i][0] - psi[i][1]*cmom[i][1]; //Act on the state by the kinetic part (real)
      cci = psi[i][0]*cmom[i][1] + psi[i][1]*cmom[i][0]; //(imaginary)
      psi[i][0] = ccr; //update real part
      psi[i][1] = cci; //and the imaginary part
    }

    fft(n,cfou,psi,-1); //Take the inverse Fourier transform to pass back into position space.

    rewind(temp); //????
    for(i=0;i<n;i++){
      fprintf(temp,"%d %lf \n",i,psi[i][0]*psi[i][0]+psi[i][1]*psi[i][1]); //Print the norm square of the evolved state.
    }
    fflush(temp); // Save into the temp file ?????

    sleep(0.1);

    fprintf(gnuplotPipe, "%s \n", "plot \"data.temp\" w l "); //Plot the state at each time increment.
    

  }
   
}


//I don't know what is going on after this point...


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
//        System.out.println("nn="+nn+" length="+length+" k="+k+
//                           " j="+j+" J="+J+" K="+K+
//                           " L="+L+" M="+M+
//                           " cfr= "+cfr+" cfi= "+cfi);
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
}
