// MC-sphere.c, adapted by D. Yllanes from the code provided in
// DJ Amit and V Martin-Mayor, "Field Theory, the Renormalization Group
//         and Critical Phenomena" (third ed., World Scientific, Singapore, 2005)
// Calculates the volume of the sphere in D dimensions by a Monte Carlo method
// The PRNG is a simple congruential generator, for simplicity
// Compile as:
// 	gcc MC-sphere.c -o sphere -lm 
// Run as
// ./sphere seed D N
// D: the space dimension , N: the number of Monte Carlo data to be generated
#include<stdio.h> 
#include<stdlib.h>  
#include<math.h>    

#define DMAX 100 // Maximum space dimension

//Definitions for the PRNG
#define FNORM (5.4210108624275218e-20)// max double such that RAND_MAX*FNORM<1
#define FRAND (randcong64()*FNORM) //A random number 0<= r <1
unsigned long long zseed; //Our random number
unsigned long long randcong64(void); //The PRNG
void init_seed(void); //Initialise seed


int main(int argc, char **argv)
{  
  double r[DMAX];    //We declare the type and the name of variables.
  double sum,sum2,modulus;
  unsigned int seed,d,N;
  int i,j;
  //we now read the arguments
  switch (argc){
    case 3:
      sscanf(argv[1],"%u",&d);
      sscanf(argv[2],"%u",&N);
      break;
    default:
      fprintf(stderr,"Usage: %s dimension N\n", argv[0]);
      exit(1);
  }

  if(d>DMAX)
    { //We haven't allocated enough memory for vector r!
      fprintf(stderr,"Recompile with DMAX larger than %d\n",d);
      exit(1);
    }

  init_seed();  //Initialize the pseudo-random sequence

  sum=sum2=0; //Put to zero the accumulation variables
  for(j=0;j<N;j++) //Monte Carlo starts!
    {
      //We now extract the random vector and calculate its modulus
      modulus=0;
      for(i=0;i<d;i++)
        {
                    r[i]=2.*FRAND-1.;
                    modulus+=r[i]*r[i];
                }
      if(modulus<1.)
                {//If the modulus is smaller than one, we are inside the sphere
		 //In order to estimate errors, we need to accumulate the sum 
		 //of our variable and of its square (in this case both are 
		 //the same)
                    sum++; //
                    sum2++;//
                }
    }
  //Divide by N, and get mean values
  sum/=N; 
  sum2/=N; 
  //Print out
  printf("Volume= %.8g +- %.8g\n",sum*pow(2.0,d), pow(2.0,d)*sqrt((sum2-sum*sum)/(N-1))); 
  printf("Exact = %.8g\n",2*pow(M_PI, 0.5*d)/d/tgamma(0.5*d)); 
  return 0; //end
}

//The 64-bit congruential pseudorandom generator
unsigned long long randcong64(void)
{
  return zseed=zseed*3202034522624059733LLU+1;
}

//We initialise the PRNG by reading from /dev/urandom
void init_seed(void)
{
  FILE *frandom;

  frandom=fopen("/dev/urandom","r");
  fread(&zseed,(size_t) 8,(size_t) 1,frandom);
  fclose(frandom);
}


