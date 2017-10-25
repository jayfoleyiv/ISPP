#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include</usr/include/complex.h>

// Function prototypes go here!


int main() {

// Here I'm going to declare my variables!

// i will be an integer
int i;
// d will be a double
// Variables for material 1
double wl, er1, ei1;
double omega, c, pi, hbar, eta1, kappa1;
double SPP_WL1, SPP_MOM1, SPP_PL1;
double complex epsilon_metal1, N_SPP1;

//=======
// Variables for material 2
//>>>>>>> b14423256b13e566b977e03495aaea05f8f21103
double er2, ei2, eta2, kappa2;
double SPP_WL2, SPP_MOM2, SPP_PL2;
double complex epsilon_metal2, N_SPP2;

//ISPP variables - only defined for material 2
double N2, k2, a, b, alpha, beta;



//=======
// ISPP variables - only defined for material 2

double angle;
//>>>>>>> b14423256b13e566b977e03495aaea05f8f21103

c = 299792458.;
pi = 3.14159265359;
hbar = 6.5821e-16;
// My file variable goes here!
FILE *fp1, *fp2, *fpw;
char filename[1000];

// Open the file DIEL/Ag_JC.txt for reading!
fp1 = fopen("DIEL/Ag_JC.txt","r");
fp2 = fopen("DIEL/Au_JC.txt","r");

printf(" PLEASE ENTER INCIDENT ANGLE!\n");
scanf("%lf",&angle);
printf(" YOU ENTERED %f AS AN INCIDENT ANGLE\n",angle);
// convert angle into radians
angle = angle * pi/180.;

printf(" PLEASE TELL ME WHERE I SHOULD WRITE THE DATA!\n");
scanf("%s",filename);
printf(" DATAWILLBE WRITTEN TO A FILE CALLED %s\n", filename);

fpw = fopen(filename, "w");

// Read one data point from file stored in variable fp
for (i=0; i<1778; i++) {

  // Reading from file 1
  // first column is wavelength
  fscanf(fp1,"%lf",&wl);
  // second column is epsilon_real
  fscanf(fp1,"%lf",&er1);
  // third column is epsilon_imaginary
  fscanf(fp1,"%lf",&ei1);

  // Reading from file 2
  // first column is wavelength
  fscanf(fp2,"%lf",&wl);
  // second column is epsilon_real
  fscanf(fp2,"%lf",&er2);
  // third column is epsilon_imaginary
  fscanf(fp2,"%lf",&ei2);

  // Define complex dielectric function of the metal 1 at the current wavelength
  epsilon_metal1 = er1 + I*ei1;

  // Define complex dielectric function of the metal 2 at the current wavelength
  epsilon_metal2 = er2 + I*ei2;

  // Calculate the complex SPP index at the current wavelength for metal 1 
  N_SPP1 = csqrt( epsilon_metal1/(epsilon_metal1 + 1) );

  // Calculate the complex SPP index at the current wavelength for metal 2
  N_SPP2 = csqrt( epsilon_metal2/(epsilon_metal2 + 1) );

  // get real part of N_SPP1 and store it to eta
  eta1 = creal(N_SPP1);

  // get imaginary part of N_SPP1 and store it to kappa
  kappa1 = cimag(N_SPP1);

  // get real part of N_SPP2 and store it to eta
  eta2 = creal(N_SPP2);

  // get imaginary part of N_SPP2 and store it to kappa
  kappa2 = cimag(N_SPP2);
//>>>>>>> b14423256b13e566b977e03495aaea05f8f21103

  // get omega from wavelength
  omega = 2*pi*c/(wl*1e-9);
  // calculate SPP wavelength

//=======
  SPP_WL1 = 2*pi*c/(omega*eta1);
  // Calculate SPP momentum
  SPP_MOM1 = hbar*omega*eta1/c;
  // Calculate SPP propagation length
  SPP_PL1 = c/(2*kappa1*omega);

  // get omega from wavelength
  omega = 2*pi*c/(wl*1e-9);
  // calculate SPP wavelength
//>>>>>>> b14423256b13e566b977e03495aaea05f8f21103
  SPP_WL2 = 2*pi*c/(omega*eta2);
  // Calculate SPP momentum
  SPP_MOM2 = hbar*omega*eta2/c;
  // Calculate SPP propagation length
  SPP_PL2 = c/(2*kappa2*omega);


// Calculate N2

alpha = eta1*sin(angle);

beta = kappa1*sin(angle);

a = (alpha*alpha)+(beta*beta)+(eta2*eta2)-(kappa2*kappa2);

b = ((kappa2-beta)*(kappa2-beta))+(eta2-alpha)*(eta2-alpha)*((kappa2+beta)*(kappa2+beta)+((eta2+alpha)*(eta2+alpha)));

N2 =((1./sqrt(2.))*sqrt(a+sqrt(b)));








  // Print wavelength, epsilon_real, epsilon_imaginary
  fprintf(fpw,"  %f %f %f\n",angle*180./pi, (omega*1e-6/c)* N2, hbar*omega); 

//=======

  // Print wavelength, epsilon_real, epsilon_imaginary
  //fprintf(fpw,"  %e  %f  %f  %f  %e  %e  %e\n",hbar*omega, wl, er, ei, SPP_WL, SPP_MOM, SPP_PL); 
//>>>>>>> b14423256b13e566b977e03495aaea05f8f21103




}

}
// Functions go here!
