// Exact Computation of Auto-Chemotaxis
// for a one cell system in 2D

// using Simpson's 3/8 rule for numerical integration

 #include <iostream.h>
 #include <math.h>
 #include <fstream.h>
 #include <string>
 using std::string; 



 #define IM1 2147483563
 #define IM2 2147483399
 #define AM (1.0/IM1)
 #define IMM1 (IM1-1)
 #define IA1 40014
 #define IA2 40692
 #define IQ1 53668
 #define IQ2 52774
 #define IR1 12211
 #define IR2 3791
 #define NTAB 32
 #define NDIV (1+IMM1/NTAB)
 #define EPS 1.2e-7
 #define RNMX (1.0-EPS)

 float ran2(long *idum) {

  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
       k = (*idum)/IQ1;
       *idum=IA1*(*idum-k*IQ1)-k*IR1;
       if (*idum < 0) *idum += IM1;
       if (j < NTAB) iv[j] = *idum;
      }
      iy=iv[0];
    }
   k=(*idum)/IQ1;
   *idum=IA1*(*idum-k*IQ1)-k*IR1;
   if (*idum < 0) *idum += IM1;
   k=idum2/IQ2;
   idum2=IA2*(idum2-k*IQ2)-k*IR2;
   if (idum2 < 0) idum2 += IM2;
   j = iy/NDIV;
   iy = iv[j]-idum2;
   iv[j] = *idum;
   if (iy < 1) iy += IMM1;
   if ((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
  }

 double px[20001] [1]; // position of cells at ith timestep and sample number k
 double py[20001] [1];

 double gradx [1];
 double grady [1];

 double vx[20002];
 double vy[20002];

 double rg[20002];
 double gr[20002];

 int i,j,k,h,celln,celln1,ith,jinit,kr,st,st1;
 double sum,tp,t,ts,a,b1,b2,c,pi,r1,r2;
 double Dp,D,Beta,alpha,alpha1;
 double x1,x2,wy,lambda,d,e,sc1,sc2;
 double sum_a,sum_ab,ss,dist,x1a,x2a,y1,y2,mm;
 double sum_b,sum1,phi,a1,c1,x,y,w,px_s,py_s,dist_a;

 int sample,sample_max,i_max,k1,k2,k3;
 int cell_max; 

 long idum=-586215;

 
 int main() {

 pi = 3.1415265359; 

 D = 0.01;
 Dp = 1;
 Beta = -1;
 alpha = 6.24 * 2;
 alpha1 = 0; 
 w = 10;
 px_s = 5;
 py_s = 5;

 i_max = 10000;

 ofstream posOut("act1.txt");


 for (kr=1; kr<=1; kr++) {

 for (i=0; i<=i_max+2; i++) {rg[i]=0; gr[i]=0;}

 ts = 0.1;
 cell_max = 0;
 
 lambda = 0.1;
 ith = 300;

 sample_max = 100;
 
 ss = 0.01;
 
 div_t divresult_a;
 

 for (sample=0; sample<=sample_max; sample++) {

 for (celln=0; celln<=cell_max; celln++) {
  
 px[0] [celln] = 0;
 py[0] [celln] = 0;

 }

 //p[0] [0] = 0;
 //p[0] [1] = 0;

 divresult_a = div (sample,1);
 if (divresult_a.rem == 0) {cout<<sample<<endl;}


 for (i=0; i<=i_max; i++) {
 
 t = i*ts;

 // Compute the gradient at the current cell position

 for (celln=0; celln<=cell_max; celln++) {

 sum = 0;
 sum1 = 0;
 
 if (t > 0) {

 jinit = 0;
 if (i>=ith) {jinit=i-ith;}
  
 for (j=jinit; j<=(i-1); j++) {
 
 sum_a = 0;
 sum_b = 0;

 tp = j*ts;
 
 a = pow(4*Dp*(t-tp),-2);
 
 for (celln1=0; celln1<=cell_max; celln1++) { 

  b1 = px[i] [celln] - px[j] [celln1];
  b2 = py[i] [celln] - py[j] [celln1];  

  c = exp( ( -(b1*b1+b2*b2)/ (4*Dp*(t-tp)) ) - (lambda*(t-tp)) );
  
  c1 = 3;

  if ((j==jinit)||(j==(i-1))) {c1 = 1;}

  sum_a = sum_a + (a*b1*c*c1);
  sum_b = sum_b + (a*b2*c*c1); 
 } // end of celln1
 
 sum = sum + (sum_a*0.375*ts);
 sum1 = sum1 + (sum_b*0.375*ts); 
 } // end of j

 } // end of if stat

 gradx [celln] = sum * (-2*Beta / pow(pi,1));
 grady [celln] = sum1 * (-2*Beta / pow(pi,1));

 } // end of celln loop

 for (celln=0; celln<=cell_max; celln++) { 

 // Generate a Gaussian random number

 do { x1 = 2.0 * ran2(&idum) - 1.0;
      x2 = 2.0 * ran2(&idum) - 1.0;
      wy = x1 * x1 + x2 * x2;
    } while ( wy >= 1.0 );

 wy = sqrt( (-2.0 * log( wy ) ) / wy );

 r1 = x1 * wy;
 r2 = x2 * wy;

 vx[i] = (pow(((2*D)/ts),0.5) * r1) + (alpha * gradx [celln]);
 vy[i] = (pow(((2*D)/ts),0.5) * r2) + (alpha * grady [celln]); 

 // Add positive chemotaxis to a fixed chemo-attractant source

 vx[i] = vx[i] + alpha1 * (2/w) * (px_s - px [i] [celln]) * exp( - ( pow(px_s - px [i] [celln],2) + pow(py_s - py [i] [celln],2) ) / w );
 vy[i] = vy[i] + alpha1 * (2/w) * (py_s - py [i] [celln]) * exp( - ( pow(px_s - px [i] [celln],2) + pow(py_s - py [i] [celln],2) ) / w );
  
 px [i+1] [celln] = px [i] [celln] + (vx[i]*ts); 
 py [i+1] [celln] = py [i] [celln] + (vy[i]*ts);
 
 // posOut<<px [i] [celln]<<" "<<py [i] [celln]<<endl; 

 } // end of celln loop 

 // dist_a = sqrt(pow(px [i+1] [celln] - px_s,2)+pow(py [i+1] [celln] -py_s,2));

 // if (dist_a <= 0.1) { rg[i+1] = rg[i+1] + 1; }
    
 rg[i+1] = rg[i+1] + px [i+1] [celln] * px [i+1] [celln] + py [i+1] [celln] * py [i+1] [celln] ;
 gr[i+1] = gr[i+1] + sqrt(pow(px [i+1] [celln],2) + pow(py [i+1] [celln],2)); 

 } // end of i loop

  
 } // end of sample loop

 } // end of kr loop
 
 for (i=0; i<=i_max; i++) {posOut<<i*ts<<" "<<(rg[i]/(sample_max+1))-(pow(gr[i]/(sample_max+1),2))<<endl;} 

 return 0;
 }


