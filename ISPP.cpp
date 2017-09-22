#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<complex.h>

// Function prototypes go here!


int main() {

// Here I'm going to declare my variables!

// i will be an integer
int i;
// d will be a double
double wl, er, ei;
double omega, c, pi, hbar, eta, kappa;
double SPP_WL, SPP_MOM, SPP_PL;
double complex epsilon_metal, N_SPP;

c = 299792458.;
pi = 3.14159265359;
hbar = 1.0545718e-34;
// My file variable goes here!
FILE *fp;

// Open the file DIEL/Ag_JC.txt for reading!
fp = fopen("DIEL/Ag_JC.txt","r");

// Read one data point from file stored in variable fp
for (i=0; i<1778; i++) {
  // first column is wavelength
  fscanf(fp,"%lf",&wl);
  // second column is epsilon_real
  fscanf(fp,"%lf",&er);
  // third column is epsilon_imaginary
  fscanf(fp,"%lf",&ei);

  // Define complex dielectric function of the metal at the current wavelength
  epsilon_metal = er + I*ei;

  // Calculate the complex SPP index at the current wavelength
  N_SPP = csqrt( epsilon_metal/(epsilon_metal + 1) );

  // get real part of N_SPP and store it to eta
  eta = creal(N_SPP);

  // get imaginary part of N_SPP and store it to kappa
  kappa = cimag(N_SPP);

  // get omega from wavelength
  omega = 2*pi*c/(wl*1e-9);
  // calculate SPP wavelength
  SPP_WL = 2*pi*c/(omega*eta);
  // Calculate SPP momentum
  SPP_MOM = hbar*omega*eta/c;
  // Calculate SPP propagation length
  SPP_PL = c/(2*kappa*omega);

  // Print wavelength, epsilon_real, epsilon_imaginary
  printf("  %e  %f  %f  %f  %e  %e  %e\n",hbar*omega, wl, er, ei, SPP_WL, SPP_MOM, SPP_PL); 

}

// Main code goes here!
printf(" Hello World! \n");


}


// Functions go here!
