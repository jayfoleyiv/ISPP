#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include</usr/include/complex.h>

// Function prototypes go here!
int *VEC_INT(int dim);
double *VEC_DOUBLE(int dim);
char *VEC_CHAR(int dim);
double complex *VEC_CDOUBLE(int dim);
void CMatMult2x2(int Aidx, double complex *A, int Bidx, double complex *B, int Cidx, double complex *C);
void TransferMatrix(double thetaI, double k0, double complex *rind, double *d,
double complex *cosL, double *beta, double *alpha, double complex *m11, double complex *m21);
double complex MaxwellGarnett(double complex epsilonM, double volumefraction, double complex epsilonI);
int Nlayer;
int polflag;
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

double er2, ei2, eta2, kappa2;
double SPP_WL2, SPP_MOM2, SPP_PL2;
double complex epsilon_metal2, N_SPP2;

//ISPP variables - only defined for material 2
double N2, k2, a, b, alpha, beta;
double angle;

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

  double complex N_SPP_Eff = MaxwellGarnett(N_SPP2, sin(angle), N_SPP1);

  // get real part of N_SPP1 and store it to eta
  eta1 = creal(N_SPP1);

  // get imaginary part of N_SPP1 and store it to kappa
  kappa1 = cimag(N_SPP1);

  // get real part of N_SPP2 and store it to eta
  eta2 = creal(N_SPP2);

  // get imaginary part of N_SPP2 and store it to kappa
  kappa2 = cimag(N_SPP2);


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

b = ((kappa2-beta)*(kappa2-beta)+(eta2-alpha)*(eta2-alpha))*((kappa2+beta)*(kappa2+beta)+((eta2+alpha)*(eta2+alpha)));

N2 =((1./sqrt(2.))*sqrt(a+sqrt(b)));








  // Print wavelength, epsilon_real, epsilon_imaginary
  fprintf(fpw,"  %f %f %f %f\n",angle*180./pi, (omega*1e-6/c)* N2, (omega*1e-6/c)*creal(N_SPP_Eff), hbar*omega); 
  fprintf(fpw,"  %f %f %f %f\n",angle*180./pi, (omega*1e-6/c)* N2, (omega*1e-6/c)* creal(N_SPP_Eff), hbar*omega); 


}

}

// Functions go here!
double complex MaxwellGarnett(double complex epsilonM, double volumefraction, double complex epsilonI) {

  double complex eps_spp1, eps_spp2, eps_eff, ri_eff;

  eps_spp2 = epsilonM*epsilonM;
  eps_spp1 = epsilonI*epsilonI;

  eps_eff = eps_spp2*(2*volumefraction*(eps_spp1 - eps_spp2) + eps_spp1 + 2*eps_spp2)/(2*eps_spp2 + eps_spp1 + volumefraction*(eps_spp2 - eps_spp1));
  ri_eff = csqrt(eps_eff);

  return ri_eff;
  
  eps_spp2= epsilonM*epsilonM;
  eps_spp1= epsilonI*epsilonI;

  eps_eff= eps_spp2*(2*volumefraction*(eps_spp1 - eps_spp2) + eps_spp1 + 2*eps_spp2)/(2*eps_spp2 + eps_spp1 + volumefraction*(eps_spp2-eps_spp1));
  ri_eff = csqrt(eps_eff);

  return ri_eff;

}

void TransferMatrix(double thetaI, double k0, double complex *rind, double *d, 
double complex *cosL, double *beta, double *alpha, double complex *m11, double complex *m21) { 
  
  int i, j, k, indx;
  double complex *kz, *phiL, *D, *Dinv, *Pl, *phil;
  double complex *EM, ctheta, tmp, *tmp2, *tmp3, c0, c1, ci, kx;

  kz   = VEC_CDOUBLE(Nlayer);
  phil = VEC_CDOUBLE(Nlayer);
  D    = VEC_CDOUBLE(4*Nlayer);
  Dinv = VEC_CDOUBLE(4*Nlayer);
  Pl   = VEC_CDOUBLE(4*Nlayer);
  EM   = VEC_CDOUBLE(4);
  tmp2 = VEC_CDOUBLE(4); 
  tmp3 = VEC_CDOUBLE(4);

  c0 = 0. + I*0.;
  c1 = 1. + I*0.;
  ci = 0. + I*1.;

  //  x-component of incident wavevector...
  //  should be in dielectric material, so the imaginary 
  //  component should be 0.
  kx = k0*rind[0]*sin(thetaI);

  //  Now get the z-components of the wavevector in each layer
  for (i=0; i<Nlayer; i++) {
     kz[i] = (rind[i]*k0)*(rind[i]*k0) - kx*kx;
     kz[i] = csqrt(kz[i]);
     // Want to make sure the square root returns the positive imaginary branch
     if (cimag(kz[i])<0.)  {
        kz[i] = creal(kz[i]) - cimag(kz[i]);
     }
   }



   //  Calculate the P matrix
   for (i=1; i<Nlayer-1; i++) {
     phil[i]=kz[i]*d[i];

     //  Upper left (diagonal 1)
     Pl[i*4] = cexp(-ci*phil[i]);  
     //  upper right (off diagonal 1)
     Pl[i*4+1] = c0;
     //  lower left (off diagonal 2)
     Pl[i*4+2] = c0;
     //  lower right (diagonal 2)
     Pl[i*4+3] = cexp(ci*phil[i]);

   }

 
   //  Calculate the D and Dinv matrices
   for (i=0; i<Nlayer; i++) {
     ctheta = kz[i]/(rind[i]*k0);
     //  p-polarized incident waves
     if (polflag==1) {  

       //  Upper left (diagonal 1)
       D[i*4] = ctheta;
       // upper right
       D[i*4+1] = ctheta;
       // lower left
       D[i*4+2] = rind[i];
       // lower right
       D[i*4+3] = -rind[i];

     } 
     //  s-polarized incident waves
     if (polflag==2) {

       // upper left
       D[i*4] = 1;
       // upper right
       D[i*4+1] = 1;
       // lower left
       D[i*4+2] = rind[i]*ctheta;
       // lower right
       D[i*4+3] = -1*rind[i]*ctheta;

     }
     //  Now compute inverse
     //  Compute determinant of each D matrix
     tmp = D[i*4]*D[i*4+3]-D[i*4+1]*D[i*4+2];
     tmp = 1./tmp;

     //printf("  tmp is %12.10f  %12.10f\n",creal(tmp),cimag(tmp));    
     Dinv[i*4]=tmp*D[i*4+3];
     Dinv[i*4+1]=-1*tmp*D[i*4+1];
     Dinv[i*4+2]=-1*tmp*D[i*4+2];
     Dinv[i*4+3]=tmp*D[i*4];
 
   }


   // Initial EM matrix
   EM[0] = c1;
   EM[1] = c0;
   EM[2] = c0;
   EM[3] = c1;
   for (i=Nlayer-2; i>0; i--) {
     CMatMult2x2(i, Pl  , i,  Dinv, 0, tmp2);
     CMatMult2x2(i, D   , 0, tmp2,  0, tmp3);
     CMatMult2x2(0, tmp3, 0, EM  ,  0, tmp2); 

     for (j=0; j<2; j++) {
       for (k=0; k<2; k++) {
          EM[2*j+k] = tmp2[2*j+k];
       }
     }
   }
   CMatMult2x2(0, EM  , Nlayer-1, D   , 0, tmp2);
   CMatMult2x2(0, Dinv, 0, tmp2, 0, EM); 


   //  Finally, collect all the quantities we wish 
   //  to have available after this function is called
   *m11 = EM[0*2+0];  //  
   *m21 = EM[1*2+0];
   *beta = creal(kx);
   *alpha = cimag(kx);
   *cosL = ctheta;

   free(kz);  
   free(phil);
   free(D);
   free(Dinv);
   free(Pl);
   free(EM);
   free(tmp2);
   free(tmp3);

}
void CMatMult2x2(int Aidx, double complex *A, int Bidx, double complex *B, int Cidx, double complex *C) {
     int i, j, k, m, n;

     double complex sum;

     for (k=0; k<2; k++) {
       for (i=0; i<2; i++) {

          sum = 0. + 0*I;
          for (j=0; j<2; j++) {

            m = 2*i + j;
            n = 2*j + k;           
            sum += A[Aidx*4+m]*B[Bidx*4+n];

          }
          
          C[Cidx*4 + (2*i+k)] = sum;
        }
      }

}

// Functions
int *VEC_INT(int dim){
  int *v,i;
  v = (int *)malloc(dim*sizeof(int));
  if (v==NULL) {
     printf("\n\nVEC_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0;
  return v;
}
double *VEC_DOUBLE(int dim){
  int i;
  double *v;
  v = (double *)malloc(dim*sizeof(double));
  if (v==NULL) {
     printf("\n\nVEC_DOUBLE: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0;
  return v;
}

double complex *VEC_CDOUBLE(int dim) {
  int i;
  double complex *v;
  v = (double complex *)malloc(dim*sizeof(double complex));
  if (v==NULL) {
     printf("\n\nVEC_CDOUBLE:  Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0 + I*0.;
  return v;
}

char *VEC_CHAR(int dim){
  char *v;
  v = (char *)malloc(dim*sizeof(char));
  if (v==NULL) {
     printf("\n\nVEC_CHAR: Memory allocation error\n\n");
     exit(0);
  }
  return v;
}
