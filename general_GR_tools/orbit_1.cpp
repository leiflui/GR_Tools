#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
using std::ifstream;
#define PI 3.14159265358979323846
#define PI 3.14159265358979323846
#define N 6



// #define a 0.0     //define spin 
// double a2 = a*a;
double L, E, kappa; //save global constants of motions
#define t0 0        //always start at t=0
// double RHor=1+sqrt(1-a2)+0.000001; //Define BH Horizon
#define H -1  //0 for photons and -1 for particles


//setup reading initial conditions from script
ifstream inFile;
#include <iostream>
#include <fstream>
double data[8];

//Define parameters globally
double r0;
double theta0;
double phi0;
double tdot0;
double rdot0;
double thetadot0;
double phidot0;
double a;
double a2;

void initial(double *y0, double *ydot0, double *nm)
{ 

  //read initial conditions from IP.txt and save in data. The initial condition describing local, isotropic emission is shown in the thesis.
  std::ifstream input("IP.txt"); 

  for (int i = 0; i < 7; i++) {
      input >> data[i];                 //save data 
      std::cout<< data[i]<<std::endl;
      }

  //rename the data, so we can remember, what is what
  r0=data[0];
  theta0=data[1];
  phi0=data[2];
  tdot0=data[3];
  rdot0=data[4];
  thetadot0=data[5];
  phidot0=data[6];
  a = data[7];
  a2 = a*a;
  // save the initial coordinates in y0
  y0[0]=r0;                      
  y0[1]=theta0; 
  y0[2]=phi0; 
  y0[3]=t0;
   
  //Define some constants 
  double r2=r0*r0;               
  double r32=r0*sqrt(r0);         
  double costheta2=cos(theta0)*cos(theta0);
  double sintheta2=sin(theta0)*sin(theta0);
  double sum=r2+a2*costheta2;
  double delta=r2-2.0*r0+a2;

  // Calculate the line element squared ds^2
  double uaua=-(1-2.0*r0/sum)*tdot0*tdot0-4.0*a*r0*sintheta2*tdot0*phidot0/sum+sum*rdot0*rdot0/delta+sum*thetadot0*thetadot0+(r2+a2+2.0*r0*a2*sintheta2/sum)*sintheta2*phidot0*phidot0; //L  

  //normalise the velocities with the line element ds
  double norm=sqrt(fabs(uaua));
   std::cout<<"Norm: "<<norm<<std::endl;
  if (H == -1){
    rdot0=rdot0/norm;         
    thetadot0=thetadot0/norm;  
    phidot0=phidot0/norm;     
    tdot0=tdot0/norm;
  }
  // std::cout<<tdot0<<"\t"<<phidot0<<std::endl;
  // put the four momentum in y0, so we have the 6 equations of motions needed to evolve the geodesics.                                                        
  y0[4]= rdot0*sum/delta;    //pr       
  y0[5]= thetadot0*sum;      //p_theta  

//set up four velocity in ydot0
  ydot0[0] = rdot0; // rdot 0         
  ydot0[1] = thetadot0; // theta_dot
  ydot0[2] = phidot0; // phi_dot
  ydot0[3] = tdot0; //tdot 0
  nm[0] = norm; //norm
  //calculate the the 3 constant of motion.
  E = 2.0*a*r0*sintheta2*phidot0/sum-(-1.0+2.0*r0/sum)*tdot0;  //energy 
  L=(sum*delta*phidot0-2.0*a*r0*E)*sintheta2/(sum-2.0*r0);     //angular momentum     
  kappa=y0[5]*y0[5]+a2*sintheta2*(E*E+H)+L*L/sintheta2;    // Kappa (Qarters constant)
 
} 


// setup the 6 parameters for the Fourth order Runge-Kutta integrator
void geodesic(double *y, double *dydlamda)  
{
  double r, theta, phi, t, pr, ptheta;
  // write out the vector y0, so it is easier to remember what parameter is what
  r=y[0];
  theta= y[1];
  phi = y[2];
  t = y[3];
  pr = y[4];
  ptheta = y[5];

  //define some parameters
  double r2=r*r;
  double twor=2.0*r;
  double cos2=cos(theta)*cos(theta);
  double sin2=sin(theta)*sin(theta);
  double sum=r2+a2*cos2;
  double delta=r2-twor+a2;
  double sd=sum*delta;

  //Calculate the four velocity and two momentum derivatives, which is used to evolve the geodesics (see thesis)
  dydlamda[0] = pr*delta/sum;     //rdot
  dydlamda[1] = ptheta/sum;       //theta dot
  dydlamda[2] = (twor*a*E+(sum-twor)*L/sin2)/sd;   //phi dot
  dydlamda[3] = E+(twor*(r2+a2)*E-twor*a*L)/sd;    //t dot
  dydlamda[4] = ((r-1.0)*((r2+a2)*H-kappa)+r*delta*H+twor*(r2+a2)*E*E-2.0*a*E*L)/sd-2.0*pr*pr*(r-1.0)/sum;  //p_r dot
  dydlamda[5] = sin(theta)*cos(theta)*(L*L/(sin2*sin2)-a2*(E*E+H))/sum;  //p_theta dot 
}



// Evolve the 6 parameters with a standard fourth order Runge-Kutta integrator. Here the parameters are found in the book Numerical Recipes in C by William Press.
void rkck (double *y, double *dydx, double h, double *yout, double *yerr) 
{
  // Constants from  Cash and Karp in the book Numerical Recipes in C by William Press.
  static const double b21 = 0.2, b31 = 3.0/40.0, b32 = 9.0/40.0, b41 = 0.3, b42 = -0.9,
    b43 = 1.2, b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0,
    b54 = 35.0/27.0, b61 = 1631.0/55296.0, b62 = 175.0/512.0,
    b63 = 575.0/13824.0, b64 = 44275.0/110592.0, b65 = 253.0/4096.0,
    c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0,
    dc1 = c1-2825.0/27648.0, dc3 = c3-18575.0/48384.0,
    dc4 = c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6 = c6-0.25; //matrix of coefficient 

  
  int i;
  int n=N;
  double ak2[6], ak3[6], ak4[6], ak5[6], ak6[6], ytemp[6];

  for (i=0;i<n;i++) ytemp[i]=y[i]+b21*h*dydx[i];    //First step
  geodesic(ytemp, ak2); //update integration and get derivatives in ak2, which is used in following step
  
  
  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);   //second step
  geodesic(ytemp, ak3);
 
  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);  //3
  geodesic(ytemp, ak4);

  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);  // 4 
  geodesic(ytemp, ak5);  

  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);  //fith step can be used to adjust the step-size so we can setup adaptive step-size
  geodesic(ytemp, ak6); 

  for (i=0;i<n;i++) yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]); //accumulate increments with proper weights (using fourth order).  output geodeisic:
 
  for (i=0;i<n;i++) yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]); // calculate/estimate the error between 5th and 4th order to see if it is below a threshold and smaller stepsize are needed
}



// Define functions max and min to find which of two parameters are larger/smaller. This will be used later to define adaptive step-sizes.
double MAX(double x, double y)
{
    if (x >= y)
      return x;
    else 
      return y;
}

double MIN(double x, double y)
{
    if (x <= y)
      return x;
    else
      return y;
}


// Give the following error message if the system cannot update the null geodesics
void nrerror(const string error_text)
{
  cerr << "Numerical Recipes run-time error..." << endl;
  cerr << error_text << endl;
  cerr << "...now exiting to system..." << endl;
  exit(1);
}



// This give us adaptive size for the integrator: if error small make h (step) larger and too large make h smaller
void rkqs (double *y, double *dydx, double &x, double htry, double eps, double *yscal, double &hdid, double &hnext)  
{
  const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
  
  int i;

  double errmax, h, htemp, xnew;

  int n=N;
  h=htry;
  double yerr[6], ytemp[6];  


  for (;;)   
    {   
      rkck(y, dydx, h, ytemp, yerr);

      errmax=0.0;
      for(i=0;i<n;i++) errmax = MAX(errmax, fabs(yerr[i]/yscal[i]));  //errmax is the error tolotance and fabs(yerr[i]/yscal[i]) is the relative error, which ever is bigger will be errmax

      errmax/=eps;  
      if(errmax<=1.0) break; //if errmax/eps<1 (less than treshhold), then we have obtained right h value. step succeded compute next step

      // Shrink step-size 
      htemp=SAFETY*h*pow(errmax, PSHRNK); // Shrinks stepsize h with by a factor (1/err)^(1/5)
      h=(h>=0.0 ? MAX(htemp,0.1*h) : MIN(htemp, 0.1*h)); //cannot vary stepsize to rapidly
    
      xnew=x+h;
      if(xnew==x) nrerror("stepsize underflow in rkqs");  //error if stepsize gets too low
    }

  if (errmax>ERRCON) hnext=SAFETY*h*pow(errmax, PGROW); //if errmax > 2e-4 (treshhold) then use higher h
  else hnext=5.0*h; 

  x+=(hdid=h);    //  new stepsize increment x + the old stepsize  (x= h_old + h_new)

  for (i=0;i<n;i++) y[i]=ytemp[i];  //save y-values and move on to next integration step

}



int main() //our main program
{ 

  int n = N;
  double y[6], dydx[6], yscal[6], y01[7]; //make a 7th element in array
  double nm[1];

  double ylaststep[6], ynextstep[6];
  double xlaststep=0, xnextstep=0;

  //Open two files, which can be used to save data.
  std::ofstream datafile;
  datafile.open("geodesic_e7.dat");  
 
   
  //Read the initial conditions t0 adaptive step size in
  double htry=0.5, eps=1e-11, hdid=0.0, hnext=0.0; 
  int i=0, q=0;
  const double TINY=1.0e-3;
  double x=0.0;
 
  //Setup the initial integration, by running 
  initial(y, dydx, nm);




   //Setting up stopping condition: only integrate until we hit the reflection surface or until photons are 1000 R_g away, where the observer is located. This is important in order to calculate the how long time the reflection surface lag behind the direct corona continuum
   //Stopping conditions thick disk
  //while (z1>zmodel&&y[0]>RHor&&y[0]<1000)  
  
  //Stopping conditions thin disk
  while (y[0]<2000&&y[0]>6)
//while y[1]>PI/2
    {

       //prepare to calculate next step and save the current position and step as ylaststep and xlaststep  
      for (i=0;i<n; i++) ylaststep[i]=y[i];
      // find the new derivatives
      xlaststep=x; 
      geodesic(y, dydx);
        

      //update y and see how big of an increment.
      for (i=0; i<n; i++) yscal[i]=fabs(y[i])+fabs(dydx[i]*htry)+TINY;

      y01[0]=y[0];
      y01[1]=y[1];
      y01[2]=y[2];
      y01[3]=L;
      y01[4]=y[4];
      y01[5]=y[5];
      y01[6]=y[3];
      
      // find the next step size
      rkqs(y, dydx, x, htry, eps, yscal, hdid, hnext);
      q++;
       //save data before and after reflection surface.
      datafile<<y[0]<<" "<<y[1]<<" "<<y[2]<<" "<<y[3]<<" "<<dydx[0]<<" "<<dydx[1]<<" "<<dydx[2]<<" "<<dydx[3]<<" "<<nm[0]<<endl;
      
      if(y[0]<6) //if r < horizon radius (particle are absorped at event horizon)
    {
      y[0]=100001;  
      break;
    }
     
      //else if (q>1500000) //give a limit to how many steps we want to integrate
      else if (q>2*150000) 
    {
      break;
    }
      
      
      //else if (y[3]>1000000) //give a limit to how many steps we want to integrate
      else if (y[3]>2*100000)
    {
      break;
    }
      
      
      htry=hnext/5; //use the next stepsize to when calculating the new integration (adaptive step-size)

    }

  
  datafile.close();
  return 0;
}

