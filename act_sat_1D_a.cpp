// Exact Computation of Auto-Chemotaxis
// for a one cell system

// 1D case

 #include <iostream.h>
 #include <math.h>
 #include <fstream.h>
 #include <stdio.h>


 #define MBIG 1000000000
 #define MSEED 161803398
 #define MZ 0
 #define FAC (1.0/MBIG)

float ran3(long *idum)
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;k++)
                        for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}
 

 double p[20002]; // position of cell at ith timestep and sample number k
 
 double sum_a[20002];
 double sum_ab[20002];

 int i,j,k,ith,jinit,s,tau,ti,x;
 double sum,tp,t,ts,a,b,c,pi,r,N;
 double Dp,D,grad,Beta,alpha,v,phi_init;
 double x1,x2,wy,lambda,frac;
 double a_conc,kappa;
 int sample,sample_max,i_max;
 

 long idum=-78;


 
 int main () {

 pi = 3.1415265359; 

 D = 1;

 Dp = 1;

 Beta = 1;

 kappa = 1;

 ts = 0.1;
 lambda = 0;

 ith = 30;

 sample_max = 1000;
 i_max = 10000;

  

 div_t divresult_a;
 
 ofstream posOut("act1a.txt");

 for (k=5; k<=5; k++) {

 if (k==1) {frac = 0;}
 if (k==2) {frac = 0.25;}
 if (k==3) {frac = 0.50;}
 if (k==4) {frac = 0.75;}
 if (k==5) {frac = 1.0;}
 if (k==6) {frac = 1.25;}
 if (k==7) {frac = 1.50;}

 alpha = frac * 2 * Dp;


 // Clear the arrays

 for (i=0; i<=i_max; i++) {
  p[i] = 0;
  sum_a[i] = 0;
  sum_ab[i] = 0; 
 }


 for (sample=1; sample<=sample_max; sample++) {
  
 p[0] = 0; // Initial Condition -- Cell at the origin

 divresult_a = div (sample,1);
 // if (divresult_a.rem == 0) {cout<<sample<<endl;}

 for (i=1; i<=i_max; i++) {
 
 t = i*ts;
 
 // The chemical profile at t = 0 is that of a delta function

 phi_init = (exp(-lambda*t)/sqrt(4*pi*Dp*t))*exp(-pow(p[i],2)/(4*Dp*t));


 // Compute the gradient at the current cell position

 sum = 0;

 if (t > 0) {

 jinit = 0;
// if (i>=ith) {jinit=i-ith;}

 for (j=jinit; j<=(i-1); j++) {
 tp = j*ts;
 a = pow(4*Dp*(t-tp),-1.5);
 b = p[i] - p[j];
 c = exp( ( -(b*b)/ (4*Dp*(t-tp)) ) - (lambda*(t-tp)) );
 sum = sum + (a*b*c*ts);
 }
 }

 grad = sum * (-2*Beta / pow(pi,0.5));

 grad = grad - (p[i]/(2*Dp*t))*phi_init;


 // Compute the absolute concentration at the current cell position

 sum = 0;
 
 jinit = 0;
// if (i>=ith) {jinit=i-ith;}

 for (j=jinit; j<=(i-1); j++) {
 tp = j*ts;
 a = pow(4*pi*Dp*(t-tp),-0.5);
 b = p[i] - p[j];
 c = exp( ( -(b*b)/ (4*Dp*(t-tp)) ) - (lambda*(t-tp)) );
 sum = sum + a*c*ts;
 } 

 a_conc = Beta * sum + phi_init;

 

 // Generate a Gaussian random number

 do { x1 = 2.0 * ran3(&idum) - 1.0;
      x2 = 2.0 * ran3(&idum) - 1.0;
      wy = x1 * x1 + x2 * x2;
    } while ( wy >= 1.0 );

 wy = sqrt( (-2.0 * log( wy ) ) / wy );
 r = x1 * wy;
 

 v = (pow((2*D/ts),0.5) * r) + kappa * alpha * (grad / a_conc);
 
 p[i+1] = p[i] + (v*ts); 


 // Calculate the Average

 sum_a[i] = sum_a[i] + p [i];
 sum_ab[i] = sum_ab[i] + pow(p [i],2);
 
 } // end of timeloop 


 } // end of sample loop 


 

 
 // Output results to a file
 

 for (i=0; i<=i_max; i++) { 

 divresult_a = div (i,1);
 if (divresult_a.rem == 0) {
 posOut<<i*ts<<" "<<(sum_ab[i]/(sample_max)) - pow((sum_a[i]/(sample_max)),2)<<endl;
 }

 }

 
 } // end of k loop

 return 0;
 }
