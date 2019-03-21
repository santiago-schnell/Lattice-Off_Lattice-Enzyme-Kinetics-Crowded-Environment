// Exact Computation of Auto-Chemotaxis
// for a one cell system

// 4 chemicals

 #include <iostream.h>
 #include <math.h>
 #include <fstream.h>

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

 double p[300002] [5]; // position of cells at ith timestep and sample number k

 double grad_1 [5];
 double grad_2 [5];
 double grad_3 [5];
 double grad_4 [5];

 double rg[300002];
 double rg1[300002];
 double rg2[300002];
 double rg3[300002];
 double rg4[300002]; 

 double gr[5];

 int i,j,k,celln,celln1,ith,jinit,kr;
 double sum,tp,t,ts,ts_a,a,b,c,pi,r,v1,v2;
 double Dp_1,Dp_2,Dp_3,Dp_4,D,Beta,alpha_1,alpha_2,alpha_3,alpha_4,v;
 double x1,x2,wy,lambda_1,lambda_2,lambda_3,lambda_4,e;
 double sum_a,sum_ab,ss,dist;
 double sum_b,sum1,phi,a1,b1,c1,epsilon_p,epsilon_c;

 int sample,sample_max,i_max;
 int cell_max; 

 long idum = -7856;

 
 int main () {

 pi = 3.1415265359; 

 D = 0.1; //0.001;

 cout<<"Dp_1"<<endl;
 //cin>>Dp_1;
 
 Dp_1 = 0.07;

 cout<<"Dp_2"<<endl;
 //cin>>Dp_2;

 Dp_2 = 0.198;

 cout<<"Dp_3"<<endl;
 //cin>>Dp_3;
 
 Dp_3 = 0.253;

 cout<<"Dp_4"<<endl;
 //cin>>Dp_4;

 Dp_4 = 0.120;

 Beta = -1;

 //alpha = 0;
 
 ofstream posOut("act6.txt");


 for (kr=1; kr<=1; kr++) {
 
 cout<<"alpha_1"<<endl;
 cin>>alpha_1;
 
 // alpha_1 = 0.245;

 cout<<"alpha_2"<<endl;
 //cin>>alpha_2;

 alpha_2 = 0.0;

 cout<<"alpha_3"<<endl;
 //cin>>alpha_3;

 alpha_3 = 0.4;
 
 cout<<"alpha_4"<<endl;
 //cin>>alpha_4;

 alpha_4 = 0.577;

 for (i=0; i<=i_max+2; i++) {rg[i]=0; rg1[i]=0; rg2[i]=0; rg3[i]=0; gr[i]=0;}

 ts = 0.1;
 

 cell_max = 0;

 cout<<"lambda_1"<<endl; 
 //cin>>lambda_1;

 lambda_1 = 0.032;
 
 cout<<"lambda_2"<<endl;
 //cin>>lambda_2;

 lambda_2 = 0.014;

 cout<<"lambda_3"<<endl; 
 //cin>>lambda_3;

 lambda_3 = 0.068;
 
 cout<<"lambda_4"<<endl;
 //cin>>lambda_4;

 lambda_4 = 0.029;

 ith = 3000;

 sample_max = 0;
 i_max = 3000;
 ss = 0.01;
 
 cout<<((alpha_1*0.936237)/(4*pow(Dp_1,1.5)*sqrt(lambda_1))) - (alpha_2/(4*pow(Dp_2,1.5)*sqrt(lambda_2))) + ((alpha_3*0.907162)/(4*pow(Dp_3,1.5)*sqrt(lambda_3))) - ((alpha_4*0.939294)/(4*pow(Dp_4,1.5)*sqrt(lambda_4)));
 
 div_t divresult_a;
 
 for (sample=0; sample<=sample_max; sample++) {
  
 p[0] [0] = 0;


 divresult_a = div (sample,1);
 if (divresult_a.rem == 0) {cout<<sample<<endl;}

 for (i=0; i<=i_max; i++) {
 
 t = i*ts;

 // Compute the gradient at the current cell position

 // Chemical A 

 sum_a = 0;

 
 if (t > 0) {

 jinit = 0;
 if (i>=ith) {jinit=i-ith;}
  
 for (j=jinit; j<=(i-1); j++) {
 
 tp = j*ts;
 
 a = pow(4*Dp_1*(t-tp),-1.5);
 
 b = p[i] [0] - p[j] [0];
 c = exp( ( -(b*b)/ (4*Dp_1*(t-tp)) ) - (lambda_1*(t-tp)) );

 sum_a = sum_a + (a*b*c)*ts;

 }

 // Chemical A is a chemo-repellent

 grad_1 [0] = sum_a * (-2*Beta / pow(pi,0.5));

 
 // Chemical B 

 sum_a = 0;

 jinit = 0;
 if (i>=ith) {jinit=i-ith;}
  
 for (j=jinit; j<=(i-1); j++) {
 
 tp = j*ts;
 
 a = pow(4*Dp_2*(t-tp),-1.5);
 
 b = p[i] [0] - p[j] [0];
 c = exp( ( -(b*b)/ (4*Dp_2*(t-tp)) ) - (lambda_2*(t-tp)) );

 sum_a = sum_a + (a*b*c)*ts;

 }

 // Chemical B is a chemo-attractant 

 grad_2 [0] = sum_a * (2*Beta / pow(pi,0.5));


 // Chemical C 

 sum_a = 0;

 

 jinit = 0;
 if (i>=ith) {jinit=i-ith;}
  
 for (j=jinit; j<=(i-1); j++) {
 
 tp = j*ts;
 
 a = pow(4*Dp_3*(t-tp),-1.5);
 
 b = p[i] [0] - p[j] [0];
 c = exp( ( -(b*b)/ (4*Dp_3*(t-tp)) ) - (lambda_3*(t-tp)) );

 sum_a = sum_a + (a*b*c)*ts;

 }

 // Chemical C is a chemo-repellent

 grad_3 [0] = sum_a * (-2*Beta / pow(pi,0.5));

 
 // Chemical D 

 sum_a = 0;

 jinit = 0;
 if (i>=ith) {jinit=i-ith;}
  
 for (j=jinit; j<=(i-1); j++) {
 
 tp = j*ts;
 
 a = pow(4*Dp_4*(t-tp),-1.5);
 
 b = p[i] [0] - p[j] [0];
 c = exp( ( -(b*b)/ (4*Dp_4*(t-tp)) ) - (lambda_4*(t-tp)) );

 sum_a = sum_a + (a*b*c)*ts;

 }

 // Chemical D is a chemo-attractant 

 grad_4 [0] = sum_a * (2*Beta / pow(pi,0.5));


 }


 // Generate a Gaussian random number

 do { x1 = 2.0 * ran2(&idum) - 1.0;
      x2 = 2.0 * ran2(&idum) - 1.0;
      wy = x1 * x1 + x2 * x2;
    } while ( wy >= 1.0 );

 wy = sqrt( (-2.0 * log( wy ) ) / wy );
 r = x1 * wy;

 if (i==100) {
  v = (pow(((2*D)/ts),0.5) * r) + (alpha_1 * grad_1[0]) + (alpha_2 * grad_2[0]) + (alpha_3 * grad_3[0]) + (alpha_4 * grad_4[0]);
 }
 else
 {
  v = (alpha_1 * grad_1[0]) + (alpha_2 * grad_2[0]) + (alpha_3 * grad_3[0]) + (alpha_4 * grad_4[0]);
 }

 p[i+1] [0] = p[i] [0] + (v*ts); 
 
 // v1 = sqrt(pow((alpha * grad [0]),2));
 
 // v2 = sqrt(pow((pow(((2*D)/ts),0.5) * r),2));
 

 //if (v1 > v2) {rg[i] = 1;} else {rg[i]=0;} 

 

 //if ( abs(p [i+1] [0] - p[i+1] [1]) <= 5) { rg[i+1] = rg[i+1] + 1; }
 
 //if ( ((p [i+1] [0]) <= 0.5)&&((p [i+1] [0]) >= 0.2) ) { rg1[i+1] = rg1[i+1] + 1; }
 //if ( ((p [i+1] [1]) <= -0.2)&&((p [i+1] [1]) >= -0.5) ) {rg2[i+1] = rg2[i+1] + 1; }
 


 rg[i+1] = rg[i+1] + pow(p [i+1] [0],2);
 rg1[i+1] = rg1[i+1] + p [i+1] [0];

 
 //rg1[i+1] = rg1[i+1] + p [i+1] [1];
 //rg2[i+1] = rg2[i+1] + p [i+1] [2];
 //rg3[i+1] = rg3[i+1] + p [i+1] [3];
 //rg4[i+1] = rg4[i+1] + p [i+1] [4];


 //sum1 = 1;

 //for (celln=0; celln<=cell_max; celln++) {
  
 //sum1 = sum1 * p[i+1] [celln];

 //}


 //gr[i+1] = gr[i+1] + sum1;
 
 } 

 }
 
 
 
 posOut<<endl;
 
 for (i=0; i<=i_max; i++) {

 divresult_a = div (i,5);
 if (divresult_a.rem == 0) {

 // posOut<<i*ts<<" "<<((rg[i]/(sample_max+1)))-pow(((rg1[i]/(sample_max+1))),2)<<endl;
 posOut<<i*ts<<" "<<p[i][0]<<endl;

 }

 }

 }
  

 return 0;
 }

 
 // Output results to a file
 
 //for (i=0; i<=i_max; i++) { 

 //sum_a = 0;
 //sum_ab = 0;

 //for (sample=0; sample<=sample_max; sample++) {
 
 //for (celln=0; celln<=cell_max; celln++) {

 // Calculate the Variance

 //sum_a = sum_a + p [i] [celln] [sample];
 //sum_ab = sum_ab + pow((p [i] [celln] [sample]),2);
 //posOut<<p [i] [sample]<<" ";
  
 //}
 //}

 //divresult_a = div (i,1);
 //if (divresult_a.rem == 0) {
 //posOut<<i*ts<<" "<<(sum_ab/((sample_max+1)*(cell_max+1))) - pow((sum_a/((sample_max+1)*(cell_max+1))),2)<<endl;
 //}
 
 
 
 
 //} 



 // Calculate Cell density

 //for (celln=0; celln<=cell_max; celln++) { 

 //for (j=90; j<=110; j++) {

 //d = (j-100) * ss;
 //e = ((j+1)-100) * ss;

 //if ((p [i] [celln] [sample] > d)&&(p [i] [celln] [sample] <= e)) {
 
 //dist = (e - p [i] [celln] [sample])/ss;
 
 //cd [i] [j+1] = cd [i] [j+1] + 1 - dist;
 //cd [i] [j] = cd [i] [j] + dist;
 //}
 //} 

 //}

 //}
 
 //}
 
 // Output Cell Density 

 //for (i=0; i<=i_max; i++) { 
 //for (j=90; j<=110; j++) {
 //posOut<<cd [i] [j]/(sample_max+1)<<" "; 
 //}
 //posOut<<endl;
 //}

 //}