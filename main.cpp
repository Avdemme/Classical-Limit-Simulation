#include "stdio.h"
#include "stdlib.h"
#include <fstream>
#include <iostream>
#include "math.h"
#include "string.h"

// random number function

double rand_uniform() {
  return (double)rand() / (double)RAND_MAX;
}

//initializes main function
using namespace std;
int main ( int argc, char** argv ) {

 
  ofstream output;

  double efield;

 
  
  
  // choosing the strength of the electric field
 
  
  
   cout << "electric field strength ( i recommend 500-2000):" << endl;
  
   cin >> efield;
  
  
   //enter number of particles

  int pnum;
   cout << "number of particles (I would keep it under 50-400 depending on the strength of the electric field ):" << endl;
   cin >> pnum;
  

  double sqpnum = sqrt(pnum);
  int isqpnum = (int)floor(sqpnum + 0.1);
  // cout << isqpnum << endl;
   
  // initializing position array into a square of particles that increases in density as number increases
 
   double pos[isqpnum][isqpnum][2];
  
  
  for ( int i = 0; i< sqpnum; i++ ){
    for ( int j = 0; j < sqpnum; j++){

       pos[i][j][0] =j * (1 /sqpnum);
       pos[i][j][1] =j * (1 /sqpnum);
    }
  }
   
  
  // initializing the center of mass
  
  double com[2];
  double comx = 0;
  double comy = 0;
  
  
  for ( int i = 0; i < sqpnum ; i++ ){
    for ( int j = 0; j < sqpnum ; j++ ){
      
       comx = comx + ((pos[i][j][0]) / pnum);
       comy = comy +((pos[i][j][1]) / pnum);
   
    }
  } 
  
  com[0] = comx;
  com[1] = comy;

  //setting up the motion of the com solo to create a classical baseline

   double continue1 = 0;
   double composi[2] = {comx , comy};
   double composf[2] = {comx , comy};
   double comveli[2] = {0,0};
   double comvelf[2] = {0,0};
   double comacci[2] = {efield,0};
   double comaccf[2] = {efield,0};
   double dt = 0.001;
   double gen1 = 0;

   while (continue1 == 0){
     composf[0] = composi[0] + 0.5 * comveli[0] * dt + 0.25 * comacci[0] * pow(dt,2);
     composf[1] = composi[1] + 0.5 * comveli[1] * dt + 0.25 * comacci[1] * pow(dt,2);
     
     if (composf[0] > 1000 ) {
       continue1 = 1;
     }
     if (composf[1] > 1000 ) {
       continue1 = 1;
     }
     comvelf[0] = comveli[0] + comacci[0] * dt;
     comvelf[1] = comveli[1] + comacci[1] * dt;
     

     comveli[0] = comvelf[0];
     comveli[1] = comvelf[1];
     composi[0] = composf[0];
     composi[1] = composi[1];
     gen1 = gen1+1;
   }

   double classicalposition[2] = {composf[0], composf[1] };
   /* ERROR CHECK
    cout << classicalposition[0] << "\t" << classicalposition[1] << endl;
    cout << efield << "\t" << gen1 << endl;
   */

  


   // Setting up the motion of the system with quantum motion
   double continue2 = 0;
   double partposi[isqpnum][isqpnum][2];
   double partposf[isqpnum][isqpnum][2];

   double partveli[isqpnum][isqpnum][2];
   double partvelf[isqpnum][isqpnum][2];

   for (int i = 0; i < isqpnum ; i++ ){
     for (int j = 0; j < isqpnum ; j++ ){
       for (int k =0; k < 2 ; k++ ){
	 partposi[i][j][k] =  pos[i][j][k];
	 partposf[i][j][k] =  pos[i][j][k];
	 partveli[i][j][k] =  0;
	 partvelf[i][j][k] =  0;
       }
     }
   }

   double dt2 = 0.001;
   double gen2 = 0;
   double comqx = 0;
   double comqy = 0;
   
   // creating the overall motion (including quantum fluctuations) for the N particle system

   while (continue2 == 0){

     // compute the position of the com
     for ( int i = 0; i < sqpnum ; i++ ){
       for ( int j = 0; j < sqpnum ; j++ ){
      
	 comqx = comqx + ((partposi[i][j][0]) / pnum);
	 comqy = comqy +((partposi[i][j][1]) / pnum);
   
       }
     } 
     // adding our end conditions
   
     if (comqx > 1000 ) {
       continue2 = 1;
     }
     if (comqy > 1000 ) {
       continue2 = 1;
     }

     // creating the quantum fluctuations
     double quantfluct[isqpnum][isqpnum][2];
     double sign1 = 1;
     double sign2 = 1;
     double fluctmag = 0;

     for ( int i = 0; i < isqpnum; i++){
       for ( int j = 0; j < isqpnum; j++ ) {
	 for ( int k = 0; k < 2 ; k++ ) {
	   sign1 = rand_uniform();
	   if ( sign1 > 0.5 ) {
	     sign2 = 1;
	   }
	   else if (sign1 < 0.5 ) {
	     sign2 = -1;
	   }

	   fluctmag = rand_uniform();

	   quantfluct[i][j][k] = sign2 * fluctmag;
	 }
       }
     }

     //setting up the electric field
     double Ef[isqpnum][isqpnum][2];
       for (int i = 0; i < isqpnum; i++){
	 for (int j = 0; j < isqpnum; j++ ){
	   for (int k = 0; k < 2; k++){
	     Ef[i][j][0] = efield;
	     Ef[i][j][1] = 0;
	   }
	 }
       }


     // setting up the interparticle forces
     double intforce[isqpnum][isqpnum][isqpnum][isqpnum][2][2];
     double rad=0;
     double const1 = -0.5;
     for( int i = 0; i < isqpnum ; i++ ) {
       for ( int j = 0; j < isqpnum; j++ ) {
	 for ( int k = 0; k<2; k++ ) {
	   for (int l = 0; l < isqpnum; l++){
	     for (int m = 0; m < isqpnum; m++){
	       for (int n = 0; n < isqpnum; n++){
		 if ( partposi[i][j][k] != partposi[l][m][n]){
		   if ( k == 0 && n == 0 ){
		     rad = partposi[i][j][k] - partposi[l][m][n];
		     if ( rad < 5.0 && rad > 0.9 ) {

		       intforce[i][j][l][m][k][n] = const1 / pow(rad,2.);
		       
		     }
		   }
		   if ( k==1 && n==1 ){
		     rad = partposi[i][j][k] - partposi[l][m][n];
		     if ( rad < 5.0 && rad > 0.9 ) {
		       
		       intforce[i][j][l][m][k][n] = const1 / pow(rad,2.);
		       
		     }
		   }
		 }
	       }
	     }
	   }
	 }
       }
     }

     // creating an array for the total interparticle force on each particle
     double totintforce[isqpnum][isqpnum][2];
     double fsign;
     double fsignind;
    
     
     for( int i = 0; i < isqpnum ; i++ ) {
       for ( int j = 0; j < isqpnum; j++ ) {
	 for ( int k = 0; k<2; k++ ) {
	   for (int l = 0; l < isqpnum; l++){
	     for (int m = 0; m < isqpnum; m++){
	       for (int n = 0; n < isqpnum; n++){
		 
		 fsignind = rand_uniform();
		 if (fsignind > 0.5 ) {
		   fsign = 1;
		 }
		 if (fsignind < 0.5 ) {
		   fsign = -1;
		 }
		 
		 totintforce[i][j][k] = fsign * (totintforce[i][j][k] + intforce[i][j][l][m][k][n]);  
	       } 
	     }
	   }
	 }
       }
     }
     
     // creating a total force array
     double totforce[isqpnum][isqpnum][2];
     for(int i = 0; i < isqpnum ; i++){
       for(int j=0; j < isqpnum ; j++){
	 for (int k = 0; k < 2 ; k++ ) {
	   totforce[i][j][k] = totintforce[i][j][k] + Ef[i][j][k];
	 }
       }
     }

     // we can finally integrate now thank god

     for ( int i=0; i < isqpnum; i++ ){
       for ( int j = 0; j < isqpnum; j++ ){
	 for ( int k = 0; k<2; k++ ){
	   partposf[i][j][k] = partposi[i][j][k]+ quantfluct[i][j][k] + 0.5 * partveli[i][j][k] * dt + 0.25 * totforce[i][j][k] * pow(dt,2.);
	   partvelf[i][j][k] = partveli[i][j][k] +  totforce[i][j][k] * dt;
	   partposi[i][j][k] = partposf[i][j][k];
	   partveli[i][j][k] = partvelf[i][j][k];
	 }
       }
     }
   }
   // finally returning the final value of the com as a test
   // cout << comqx << "\t" << comqy << endl;
   // determining the distance error
   double distance;
   distance = sqrt(pow((comqx - classicalposition[0]),2) + pow((comqy - classicalposition[1]),2));
   cout <<efield << "\t" << pnum << "\t" << distance << endl;

   output.open ( "output.txt", ios::app);
   output << pnum << "\t" << distance << endl;
   output.close();

   
  
   

  return 1;
}
