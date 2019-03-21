// Exact Computation of Auto-Chemotaxis
// for a one cell system in 2D

// using Simpson's 3/8 rule for numerical integration

 #include <iostream.h>
 #include <math.h>
 #include <fstream.h>
 #include <string>
 using std::string; 

 #include <stdio.h>

 #include <stdlib.h>



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


 double px[10001] [1]; // position of cells at ith timestep and sample number k
 double py[10001] [1];

 double gradx [1];
 double grady [1];

 double rg[10002];
 double rg1[10002];
 double rg2[10002];
 double rg3[10002];
 double rg4[10002]; 
 double rg5[10002];
 double rg6[10002];
 double rg7[10002];
 double rg8[10002];
 double rg9[10002];
 double rg10[10002];
 double rg11[10002];
 double rg12[10002];
 double rg13[10002];
 double rg14[10002]; 
 double rg15[10002];
 double rg16[10002];
 double rg17[10002];
 double rg18[10002];
 double rg19[10002];
 double rg20[10002];


 double vx[10002];
 double vy[10002];

 double gr[10002];
 double gr1[10002];

 int i,j,k,h,celln,celln1,ith,jinit,kr,st,st1;
 double sum,tp,t,ts,a,b1,b2,c,pi,r1,r2;
 double Dp,D,Beta,alpha;
 double x1,x2,wy,lambda,d,e,sc1,sc2;
 double sum_a,sum_ab,ss,dist,x1a,x2a,y1,y2,mm;
 double sum_b,sum1,phi,a1,c1,x,y;

 int sample,sample_max,i_max,k1,k2,k3;
 int cell_max; 

 string s;
 char  buffer[200];

 long idum=-586215;

 
 int main () {

 pi = 3.1415265359; 

 D = 0.01;
 Dp = 1;
 Beta = -1;

 const char *ptr1 = 0;
 
 ofstream posOut("act1.txt");


 for (kr=1; kr<=1; kr++) {

 cin>>alpha;
 cin>>st;
 cin>>st1;

 i_max = 600;

 for (i=0; i<=i_max+2; i++) {rg[i]=0; rg1[i]=0; rg2[i]=0; gr[i]=0;}

 ts = 0.1;
 cell_max = 0;
 
 lambda = 0.1;
 ith = 300;

 sample_max = 10000;
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
 //if (i>=ith) {jinit=i-ith;}
  
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
 }
 
 sum = sum + (sum_a*0.375*ts);
 sum1 = sum1 + (sum_b*0.375*ts); 
 }

 }


 gradx [celln] = sum * (-2*Beta / pow(pi,1));
 grady [celln] = sum1 * (-2*Beta / pow(pi,1));


// phi = exp(-(p[i] [celln]*p[i] [celln])/(1+4*Dp*t)) * (exp(-lambda*t)/sqrt(1+4*Dp*t)) + Beta * sum1;
 
// grad [celln] = grad [celln] / phi;

 }

 for (celln=0; celln<=cell_max; celln++) { 

 // Generate a Gaussian random number

 do { x1 = 2.0 * ran3(&idum) - 1.0;
      x2 = 2.0 * ran3(&idum) - 1.0;
      wy = x1 * x1 + x2 * x2;
    } while ( wy >= 1.0 );

 wy = sqrt( (-2.0 * log( wy ) ) / wy );

 r1 = x1 * wy;
 r2 = x2 * wy;

 vx[i] = (pow(((2*D)/ts),0.5) * r1) + (alpha * gradx [celln]);
 vy[i] = (pow(((2*D)/ts),0.5) * r2) + (alpha * grady [celln]); 

 
 px[i+1] [celln] = px[i] [celln] + (vx[i]*ts); 
 py[i+1] [celln] = py[i] [celln] + (vy[i]*ts);
 

 } 

  rg[i] = rg[i] + vx [i] * vx [i];
  rg1[i] = rg1[i] + vx [i] * vx [i-st-st1];
  rg2[i] = rg2[i] + vx [i] * vx [i-2*st-st1];
  rg3[i] = rg3[i] + vx [i] * vx [i-3*st-st1];
  rg4[i] = rg4[i] + vx [i] * vx [i-4*st-st1];
  rg5[i] = rg5[i] + vx [i] * vx [i-5*st-st1];
  rg6[i] = rg6[i] + vx [i] * vx [i-6*st-st1];
  rg7[i] = rg7[i] + vx [i] * vx [i-7*st-st1];
  rg8[i] = rg8[i] + vx [i] * vx [i-8*st-st1];
  rg9[i] = rg9[i] + vx [i] * vx [i-9*st-st1];
  rg10[i] = rg10[i] + vx [i] * vx [i-10*st-st1];
  rg11[i] = rg11[i] + vx [i] * vx [i-11*st-st1];
  rg12[i] = rg12[i] + vx [i] * vx [i-12*st-st1];
  rg13[i] = rg13[i] + vx [i] * vx [i-13*st-st1];
  rg14[i] = rg14[i] + vx [i] * vx [i-14*st-st1];
  rg15[i] = rg15[i] + vx [i] * vx [i-15*st-st1];
  rg16[i] = rg16[i] + vx [i] * vx [i-16*st-st1];
  rg17[i] = rg17[i] + vx [i] * vx [i-17*st-st1];
  rg18[i] = rg18[i] + vx [i] * vx [i-18*st-st1];
  rg19[i] = rg19[i] + vx [i] * vx [i-19*st-st1];
  rg20[i] = rg20[i] + vx [i] * vx [i-20*st-st1];

  gr[i] = gr[i] + vx [i];
  
  //if (i>=10) {

  //rg[i+1] = rg[i+1] + pow(px[i+1] [0],5);

  //rg1[i+1] = rg1[i+1] + px[i+1] [0] * px[i] [0];
  //rg2[i+1] = rg2[i+1] + px[i+1] [0] * px[i-1] [0];
  //rg3[i+1] = rg3[i+1] + px[i+1] [0] * px[i-2] [0];
  //rg4[i+1] = rg4[i+1] + px[i+1] [0] * px[i-3] [0];
  //rg5[i+1] = rg5[i+1] + px[i+1] [0] * px[i-4] [0];
  //rg6[i+1] = rg6[i+1] + px[i+1] [0] * px[i-5] [0];
  //rg7[i+1] = rg7[i+1] + px[i+1] [0] * px[i-6] [0];
  //rg8[i+1] = rg8[i+1] + px[i+1] [0] * px[i-7] [0];
  //rg9[i+1] = rg9[i+1] + px[i+1] [0] * px[i-8] [0];
 
  //gr[i+1] = gr[i+1] + px[i+1] [0];
  
  //}
 
 }

  
 } 


 mm = (rg[i_max]/(sample_max+1))-pow(gr[i_max]/(sample_max+1),2);

 // Calculate auto-correlation
  
  posOut<<0<<" "<<((rg[i_max]/(sample_max+1))-pow(gr[i_max]/(sample_max+1),2))/mm<<endl;

  posOut<<st+st1<<" "<<((rg1[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<2*st+st1<<" "<<((rg2[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-2*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<3*st+st1<<" "<<((rg3[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-3*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<4*st+st1<<" "<<((rg4[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-4*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<5*st+st1<<" "<<((rg5[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-5*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<6*st+st1<<" "<<((rg6[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-6*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<7*st+st1<<" "<<((rg7[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-7*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<8*st+st1<<" "<<((rg8[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-8*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<9*st+st1<<" "<<((rg9[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-9*st-st1]/(sample_max+1)))/mm<<endl;
 
  posOut<<10*st+st1<<" "<<((rg10[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-10*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<11*st+st1<<" "<<((rg11[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-11*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<12*st+st1<<" "<<((rg12[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-12*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<13*st+st1<<" "<<((rg13[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-13*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<14*st+st1<<" "<<((rg14[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-14*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<15*st+st1<<" "<<((rg15[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-15*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<16*st+st1<<" "<<((rg16[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-16*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<17*st+st1<<" "<<((rg17[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-17*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<18*st+st1<<" "<<((rg18[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-18*st-st1]/(sample_max+1)))/mm<<endl;
 
  posOut<<19*st+st1<<" "<<((rg19[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-19*st-st1]/(sample_max+1)))/mm<<endl;

  posOut<<20*st+st1<<" "<<((rg20[i_max]/(sample_max+1))-(gr[i_max]/(sample_max+1))*(gr[i_max-20*st-st1]/(sample_max+1)))/mm<<endl;

   
  // for (i=0; i<=i_max; i++) {posOut<<i<<" "<<pow(gr[i]/(sample_max+1),2)+pow(gr1[i]/(sample_max+1),2)<<endl;}
 
 i = 0;

 if (i==1) {


 // Create movie of chemical concentration plot

 // loop over number of frames

 for (k3=1; k3<=100; k3++) {

 i = k3*30;

 cout<<"Movie frame ... "<<k3<<endl;

 h += sprintf( buffer, "%d", k3 );


 s = "test";
 s = s + buffer;
 s = s + ".grd";
 
 ptr1= s.data ( );

 cout<<ptr1<<endl;


 ofstream posOut(ptr1);

 sum = 0;

 for (k2=0; k2<=100; k2++) {
 for (k1=0; k1<=100; k1++) {

 sum = sum + 1;

 }
 }

 posOut<<"{TITL Grid Vers 2 3}"<<endl;
 posOut<<"{CART -90.0 45.0 0 0 METR 0}"<<endl;
 posOut<<"{DPAL "<<sum<<endl;

 x1a = -50;
 x2a = 200;

 y1 = -160;
 y2 = 40; 

 sc1 = (x2a - x1a)*0.01;
 sc2 = (y2 - y1)*0.01; 

 for (k2=0; k2<=100; k2++) {
 for (k1=0; k1<=100; k1++) {

 x = x1a + (x2a - x1a)*0.01*k1;
 y = y1 + (y2 - y1)*0.01*k2;

 //posOut<<x<<" "<<y<<endl;

 sum_a = 0; 

 for (j=0; j<=(i-1); j++) {
 
 t = i*ts;
 tp = j*ts;
 
 a = pow(4*pi*Dp*(t-tp),-1);
  
 b1 = x - px[j] [0];
 b2 = y - py[j] [0];  

 c = exp( ( -(b1*b1+b2*b2)/ (4*Dp*(t-tp)) ) - (lambda*(t-tp)) );
  
 //c1 = 3;

 //if ((j==0)||(j==(i_max-1))) {c1 = 1;}

 sum_a = sum_a + (a*c*ts); 
 
 } 

 posOut<<"( "<<x<<","<<y<<" ) "<<log(sum_a+1e-20)<<endl;
  
 }
 }
 
 posOut<<"}"<<endl;
 posOut<<"{ENDF}"<<endl;

 posOut.close(); 

 } // end of loop over movie frames


 }
 }

 return 0;
 }


