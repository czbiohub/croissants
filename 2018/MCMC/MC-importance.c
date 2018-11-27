//Integrate c*x^(c-1)*cos(x) for 0<x<1
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
  double  sum_b, sum2_b;
  double  sum_i, sum2_i;
  unsigned int N;
  int j;
  double c;
  //we now read the arguments
  switch (argc){
    case 3:
      sscanf(argv[1],"%lf",&c);
      sscanf(argv[2],"%u",&N);
      break;
    default:
      fprintf(stderr,"Usage: %s c N\n", argv[0]);
      exit(1);
  }


  init_seed();  //Initialize the pseudo-random sequence

  sum_b=sum2_b=0;
  sum_i=sum2_i=0;
  double random;
  double f,y;
  for(j=0;j<N;j++) //Monte Carlo starts!
    {
      random = FRAND;
      //Brute force: we just evaluate f(random) and average it
      f = c*pow(random,c-1)*cos(random);
      sum_b += f;
      sum2_b += f*f;

      //Importance: we first generate a random variable r with 
      //probability density p(x)=c*x^(c-1), we then evaluate cos(r) 

      //F(y) = int_0^y dx  c*x^(c-1) = y^c
      //F^-1(random) is distributed with p(x) 
      y = pow(random, 1/c);
      f = cos(y);

      sum_i += f;
      sum2_i += f*f;
    }
  //Divide by N, and get mean values
  sum_b/=N; 
  sum2_b/=N; 
  sum_i/=N; 
  sum2_i/=N; 
  //Print out
  printf("Brute force= %.8g +- %.8g\n",sum_b, sqrt((sum2_b-sum_b*sum_b)/(N-1))); 
  printf("Importance= %.8g +- %.8g\n",sum_i, sqrt((sum2_i-sum_i*sum_i)/(N-1))); 
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


