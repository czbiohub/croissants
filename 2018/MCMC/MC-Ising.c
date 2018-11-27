//Sample program for the simulation of the Ising model in 2D
// adapted by D. Yllanes from the code provided in
// DJ Amit and V Martin-Mayor, "Field Theory, the Renormalization Group
//         and Critical Phenomena" (third ed., World Scientific, Singapore, 2005)

//We implement Heat-bath dynamics, single-cluster and Swendsen-Wang updates.
//When compiling, define the system size L
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h> 

//Global definitions

//We implement both the single-spin dynamics (heat bath) 
//and two versions of cluster dynamics
//#define HEATBATH        //Uncomment the dynamical rules you want to use
//#define SWENDSENWANG
//#define WOLFF

//The macro UPDATE should contain the name of the updating routine.
//It is set according to the previous choice. We need to nest several
//#if to make sure that only one updating rule is choosen.

#ifdef HEATBATH 
#define UPDATE Heat_bath() 
#else
#ifdef SWENDSENWANG
#define UPDATE SwendsenWang()
#else
#define UPDATE Wolff()
#define SINGLECLUSTER //To print-out information related to the update
#endif
#endif

//define L when compiling
#define Lx   L        //Lattice size in the X-direction
#define Ly   L        //Lattice size in the Y-direction

#define V (Lx*Ly)    //Lattice volume    
#define NOBS  3      //Number of observables measured during the simulation
#define NRES  3      //Number of partial averages printed on the screen

#define MAXIT 10000  //Maximum length of the data block

//Definitions for the PRNG
#define TWO64MINUS1 18446744073709551615ULL //2^64 - 1 
#define TWO64 18446744073709551616.   //2^64 
#define TWO63  9223372036854775808ULL  //2^63 
#define FNORM (5.4210108624275218e-20)// max double such that RAND_MAX*FNORM<1
#define FRANDOM (randcong64()*FNORM) //A random number 0<= r <1
#define MYRANDOM randcong64()
unsigned long long zseed; //Our random number
unsigned long long randcong64(void); //The PRNG
void init_seed(void); //Initialise seed


typedef struct {unsigned char x,y;} COOR; //Coordinates for the cluster method
                                          //We are limited to Lx=Ly=256
                                          //If you want more, set type to int.

typedef struct   //Data_structure  describing the simulation
{
 unsigned int itmax, // Number of measurements per block
    mesfr,           // frecuency of measurements
    nbin,            // number of data block to be generated
    itcut,           // next block to be generated                     
    nterm,           // number of thermalization steps (without measurements)
    flag,            // starting configuration: 0(random), 1(cold),2(backup)
    lx,              // lattice dimension along X  (just to check!) 
    ly;              // lattice dimension along Y  (just to check!)
   double beta       // The adimensional combination beta J
    ;   
} s_data;

#define NDAT_INT 8 //Number of integer-type data in structre of type data
#define NDATA_DOUBLE 1 //Number of double-type data  structre of type data


//Declaration of global variables 

//Variables for input-output:

s_data data; //The data describing the simulation
char directory[1024];     //The directory where output shold be written
double vdata[NOBS][MAXIT];    //data saved to disk for detailed analysis
double obs[NOBS];             //Measurements taken in the current configuration
double results[NRES];         //Partial average

//Variables for simulation:

char u[V]; //the spin is 1 or -1, you save memory using the smallest type
unsigned long long prob[5]; // Heat-bath probability table. Local field: -4,-2,0,2,4
unsigned long long  p;        // [1-exp(- 2 beta J)]*2^64
char my_spin,my_new_spin; //Variables for the cluster method
double Nspin,Nspin2,Ncluster; //cluster statistic
double averageNspin,averageNspin2,averageNcluster; 
int plus_x[Lx]; //Addressing: nearest neighbor in the positive x-direction,etc
int plus_y[Ly];
int minus_x[Lx];
int minus_y[Ly];
int fc_plus_y[Ly]; //Addresing for the cluster method, in function flip_cluster
int fc_minus_y[Ly];

//Variables for measurements:

double cos_x[Lx]; // cos(2 pi x/Lx), etc. Needed for second-momentum
double cos_y[Ly]; // correlation-length
double sin_x[Lx];
double sin_y[Ly];

//Function declarations

void Heat_bath(void);        //Heat-bath update
void Wolff(void);            //Single-cluster dynamics
void SwendsenWang(void);     //All clusters dynamics
void flip_cluster(COOR,int); //The core of the cluster methods
void read_data(void);        //Read starting data
void back_up(void);          //Back up
void read_back_up(void);     //Restart the simulation
void measurements(void);         //Take measurements
void write_measurements(int);    //Write measurements to disk
void check_file_output(int); //
void show_results(int);     //Show partial averages
void print_and_exit(char *format, ...); //Print out an error message and exits
void Init(void);//Initialize spins, probabilities, etc. 
void Init_Boltz(double); //Initialize the probability table

int main(void)
{
  int it,iblo,imes,i;
  
  read_data(); //Read all control variables
  Init(); //Initialize (or read backup) everithing:spin, addressing,PRNG,...

  for(it=0;it<data.nterm;it++) //Do nterm steps for thermalization
     UPDATE;  //The selected dynamical rule


  for(iblo=data.itcut;iblo<data.nbin;iblo++) //Next block of measurement
    {

      check_file_output(iblo);// Output data-dile should not exist!
      for(i=0;i<NRES;i++)//Put to zero accumulation variables
        results[i]=0;
#ifdef SINGLECLUSTER  //We put to zero accumulators for cluster statistics
      averageNspin=0; 
#endif
#ifdef SWENDSENWANG
      averageNcluster=0;
      averageNspin2=0;
#endif
      for(imes=0;imes<data.itmax;imes++)
        {
          for(it=0;it<data.mesfr;it++)//Monte Carlo loop without measurements
            UPDATE;

          measurements(); //Take your measurements
          for(i=0;i<NOBS;i++) 
            vdata[i][imes]=obs[i]; //Prepare to write

          //Accumulate for partial average.
          results[0]+=obs[0];        //The (minus) energy density
          results[1]+=obs[1]*obs[1]; //M^2
          results[2]+=obs[2];        //F
        }
      for(i=0;i<NRES;i++)
        results[i]/= data.itmax; //Normalize partial averages

#ifdef SINGLECLUSTER  //We normalize accumulators for cluster statistics
      averageNspin/=data.itmax*data.mesfr; 
#endif
#ifdef SWENDSENWANG
      averageNcluster/=data.itmax*data.mesfr; 
      averageNspin2/=data.itmax*data.mesfr; 
#endif

      write_measurements(iblo); //Write measurements to disk
      data.itcut=iblo+1; //The backup knows which blocks should be done next.
      back_up();          //Save crucial information
      show_results(iblo);     //Print to screen  partial averages
    }

  return 0; // We're done!
}

void Heat_bath(void)
{
  int j,x,y,neigh_px,neigh_py,neigh_mx,neigh_my;
  int h;

  j=0; //We update the lattice sequentially, j=x+Lx*y
  for(y=0;y<Ly;y++)
    {
      neigh_py=plus_y[y]; //It is better to use local memory
      neigh_my=minus_y[y];
    for(x=0;x<Lx;x++)
      {
        neigh_px=plus_x[x];
        neigh_mx=minus_x[x];

        //The probability table is indexed according to h=(4+local-field)/2
        h=(4+u[j+neigh_px]+u[j+neigh_mx]+u[j+neigh_py]+u[j+neigh_my])>>1;
        if(MYRANDOM<prob[h]) 
          u[j]=1;
        else
          u[j]=-1;

        j++;
      }
    }
}

void Wolff(void)
{
  //We pick a random site, build the cluster to which it belongs and flip it
  int site;
  COOR position;
  position.x=FRANDOM*Lx;  //The coordinates of the random site
  position.y=FRANDOM*Ly;
  site=position.x+Lx*position.y;
  my_spin=u[site];        //The current value of the spins in the cluster
  my_new_spin=-my_spin;   //The new value of the spins in the cluster
  Nspin=0; //We count the number of spins in the cluster
  flip_cluster(position,site);
  averageNspin+=Nspin;    //We average the number of spins in cluster
}

void SwendsenWang(void)
{
  //We go through the lattice. An spin not already asigned to
  //a cluster should be 1 or -1, hence spin+2 should be 3 or 1.
  //In binary 3 is 11 while 1 is 01. Hence the bitwise "and" of spin+2 with 3
  //is larger than zero.
  //A spin=2  has been asigned to a cluster that will end-up being -1
  //hence  spin+2=4, in binary 100. Its bitwise "and" with 3 (binary (011)) 
  //is zero, hence we know that the spin already belongs to a cluster.
  //A spin=6  has been asigned to a cluster that will end-up being +1
  //hence  spin+2=8, in binary 1000. Its "and" with 3 (binary (0011)) is zero.

  int site;
  COOR position; 
  site=0;
  Ncluster=0;
  Nspin2=0; //We average the square of the numper of spins in each cluster
  for(position.y=0;position.y<Ly;position.y++)
   for(position.x=0;position.x<Lx;position.x++)
     {
       if((u[site]+2)&3)
         { //This spin still does not belong to a cluster
           my_spin=u[site];
           if(MYRANDOM<TWO63)
             my_new_spin=2;
           else
             my_new_spin=6;
           Ncluster++; //We start to build a new cluster
           Nspin=0;
           flip_cluster(position,site);
           Nspin2+=Nspin*Nspin;//Average square of number of spins in cluster
         }
       site++;
     }
  averageNcluster+=Ncluster;
  averageNspin2+=Nspin2;

  for(site=0;site<V;site++) //We undo the trick!
    {
      if(u[site]==2)
        u[site]=-1;
      else
        u[site]=1;
    }
}


//This is a recursive rutine trying to extend the cluster as far as a possible
//in each direction of growth
void flip_cluster(COOR position,int site)
{
  int site_new,shift,shift_y;

  u[site]=my_new_spin; //From now on, this site already belongs to a cluster

  Nspin++;             //There is a new spin in cluster 

  shift=plus_x[position.x]; //We try to grow along positive X
  site_new=site+shift; 
  if(u[site_new]==my_spin) //Spin sign in  forward X direction correct?

    if(MYRANDOM<p)         //OK, the cluster growths
      {
        position.x+=shift; //The new coordinates of the cluster grower
        flip_cluster(position,site_new);
        position.x-=shift; //The cluster grower is back again
      }

  shift=minus_x[position.x]; //We try to grow along negative X
  site_new=site+shift; 
  if(u[site_new]==my_spin) //Spin sign in backward X direction correct?
    if(MYRANDOM<p)         //OK, the cluster growths
      {
        position.x+=shift; //The new coordinates of the cluster grower
        flip_cluster(position,site_new);
        position.x-=shift;
      }


  shift=plus_y[position.y]; //We try to grow along positive Y
  shift_y=fc_plus_y[ position.y]; //To move the cordinate y
  site_new=site+shift; 
  if(u[site_new]==my_spin) //Spin sign in  forward Y direction correct?
    if(MYRANDOM<p)         //OK, the cluster growths
      {
        position.y+=shift_y; //The new coordinates of the cluster grower
        flip_cluster(position,site_new);
        position.y-=shift_y;
      }

  shift=minus_y[position.y]; //We try to grow along negative Y
  shift_y=fc_minus_y[ position.y];
  site_new=site+shift; 
  if(u[site_new]==my_spin) //Spin sign in backward Y direction correct?
    if(MYRANDOM<p)         //OK, the cluster growths
      {
        position.y+=shift_y; //The new coordinates of the cluster grower
        flip_cluster(position,site_new);
        position.y-=shift_y;
      }
}

void measurements(void)
{
  char spin;

  int site,neigh_px,neigh_py,x,y;
  int Mx[Lx],My[Ly];
  double m00,m01r,m01i,m10r,m10i;
  double ener;

  //We put to zero accumulation variables
  ener=0;
  m00=m01r=m01i=m10r=m10i=0;
  for(x=0;x<Lx;x++)
    Mx[x]=0;
  for(y=0;y<Ly;y++)
    My[y]=0;

  //Take measurements!
  site=0;
  for(y=0;y<Ly;y++)
    {
      neigh_py=plus_y[y]; //It is better to use local memory
    for(x=0;x<Lx;x++)
      {
        neigh_px=plus_x[x];
        spin=u[site];
        ener+=spin*(u[site+neigh_px]+u[site+neigh_py]); //(minus) energy
        Mx[x]+=spin; //To save multiplications for the Fourier transform
        My[y]+=spin; 
        m00+=spin;   //The magnetization;
        site++;
      }
    }
  obs[0]=ener/((double) V);
  obs[1]=m00/((double) V);

  //For the Fourier transform, we take k along the longest
  //lattice direction!
  if(Lx==Ly)
    {

      for(x=0;x<Lx;x++) //Take Fourier transform for k=(2 pi /Lx,0)
        {
          m10r+=Mx[x]*cos_x[x];
          m10i+=Mx[x]*sin_x[x];
        }
      m10r/=((double) V);
      m10i/=((double) V);
      for(y=0;y<Ly;y++) //Take Fourier transform for k=(0,2 pi /Ly)
        {
          m01r+=My[y]*cos_y[y];
          m01i+=My[y]*sin_y[y];
        }
      m01r/=((double) V);
      m01i/=((double) V);
      obs[2]=0.5*(m01r*m01r+m01i*m01i+m10r*m10r+m10i*m10i);
    }
  else 
    {
      if(Lx<Ly)
        {
          for(y=0;y<Ly;y++) //Take Fourier transform for k=(0,2 pi /Ly)
            {
              m01r+=My[y]*cos_y[y];
              m01i+=My[y]*sin_y[y];
            }
          m01r/=((double) V);
          m01i/=((double) V);
          obs[2]=m01r*m01r+m01i*m01i;
        }
      else //Lx>Ly
        {
          for(x=0;x<Lx;x++) //Take Fourier transform for k=(2 pi /Lx,0)
            {
              m10r+=Mx[x]*cos_x[x];
              m10i+=Mx[x]*sin_x[x];
            }
          m10r/=((double) V);
          m10i/=((double) V);
          obs[2]=m10r*m10r+m10i*m10i;
        }
    }
}

//Shows partial averages on the screen, to check that everithing is smooth
void show_results(int blo)
{
  char name[1024];
  FILE *Foutput;
  sprintf(name,"%s/partial_averages.dat",directory);

  if(NULL==(Foutput=fopen(name,"a")))
    print_and_exit("I could not open %s\n",name);

  printf("Block %d: e=%.8g, M^2=%8g, F=%8g",blo,results[0],
         results[1],results[2]);
  fprintf(Foutput,"%d %.8g %.8g %.8g",blo,results[0],results[1],results[2]);
#ifdef SINGLECLUSTER
  printf("\nAverage-cluster=%.8g",averageNspin);
  fprintf(Foutput," %.8g",averageNspin);
#endif
#ifdef SWENDSENWANG
  printf("\nAverage-cluster^2=%.8g, average # clusters=%.8g",
         averageNspin2,averageNcluster);
  fprintf(Foutput," %.8g %.8g",averageNspin2,averageNcluster);
#endif
  printf("\n");
  fprintf(Foutput,"\n");

  fflush(stdout); //Print to screen now!!
  fflush(Foutput);//Print to disk now!!
  
}

//Writes to disk all measurements taken during block blo
void write_measurements(int blo)
{
  int i;
  char name[1024];
  FILE *Foutput;
  sprintf(name,"%s/OUT%03d.DAT",directory,blo);
  if(NULL==(Foutput=fopen(name,"wb")))
    print_and_exit("I could not open %s\n\a",name);

  fwrite(&data,sizeof(data),1,Foutput);
  for(i=0;i<NOBS;i++)    
    if(! fwrite(&vdata[i][0],sizeof(double)*data.itmax,1,Foutput))
      print_and_exit("Error writting observable %d in block %d\n\a",i,blo);
  fclose(Foutput);

}

//Read the variables that determine the run from input file
void read_data(void)
{

  char coso[1024];
  int j;
  unsigned int * ptdata_int;
  double * ptdata_real;
  FILE *Finput;
  Finput=fopen("input","r");
  if (Finput==NULL)
    print_and_exit("I could not open 'input'.\n");

  fgets(coso,1024,Finput);
  sscanf(coso,"%s",directory);

  for (j=0,ptdata_int=&data.itmax;j<NDAT_INT;j++)
    {
        fgets(coso,1024,Finput);
        sscanf(coso,"%d",ptdata_int++);
    }
  for (j=0,ptdata_real=&data.beta;j<NDATA_DOUBLE;j++)
    {
        fgets(coso,1024,Finput);
        sscanf(coso,"%lf",ptdata_real++);
    }

  fclose(Finput);
  //Print to screen what has been read!
  printf("Writing directory %s\n",directory);
  printf("itmax  %u \n",data.itmax);
  printf("mesfr  %u \n",data.mesfr);
  printf("nbin   %u \n",data.nbin);
  printf("itcut  %u \n",data.itcut);
  printf("nterm  %u \n",data.nterm);
  printf("flag   %u \n",data.flag);
  printf("lx     %u \n",data.lx);
  printf("ly     %u \n",data.ly);
  printf("beta   %lf\n",data.beta);

  if(data.itmax>MAXIT)
    print_and_exit("Recompile with MAXIT>%u",data.itmax);
  if(data.lx!=Lx)
    print_and_exit("Lx should be %u, rather than %u\n\a",Lx,data.lx);
  if(data.ly!=Ly)
    print_and_exit("Ly should be %u, rather than %u\n\a",Ly,data.ly);
}


//Writes all information needed to continue the run
void back_up(void)
{
  FILE *Fconfig;
  char name[1024],name_dollar[1024],name_old[1024];

  //We keep in disk two backup files: conf and conf.old.
  //Everithing is first written in conf.$$$.
  // Should everithing work fine,
  //We put the old conf in conf.old and conf.$$$ in conf

  sprintf(name_dollar,"%s/conf.$$$",directory);
  sprintf(name,"%s/conf",directory);
  sprintf(name_old,"%s/conf.old",directory);


  if(NULL==(Fconfig=fopen(name_dollar,"wb")))
    print_and_exit("I could not open %s\n\a",name_dollar);
  fwrite(&data,sizeof(data),1,Fconfig); //data structure

  //Now write the data for the PRNG
  fwrite(&zseed,sizeof(zseed),1,Fconfig);
  //Finally, the configuration
  if(!fwrite(&u,sizeof(u),1,Fconfig))
    {
      fclose(Fconfig);
      print_and_exit("Error writing backup file\n");
    }
  fclose(Fconfig);  
  remove(name_old);
  rename(name,name_old);
  rename(name_dollar,name);
}


//Reads all information needed to continue the run
void read_back_up(void)
{
  FILE *Fconfig;
  char name[1024];
  s_data datab;

  sprintf(name,"%s/conf",directory);
  
  if(NULL==(Fconfig=fopen(name,"rb")))
    print_and_exit("I could not open %s\n\a",name);

  fread(&datab,sizeof(datab),1,Fconfig); //data structure

  if(data.lx!=Lx)
    print_and_exit("Error in backup file\n!"
                   "Lx should be %u, rather than %u\n\a",Lx,datab.lx);
  if(data.ly!=Ly)
      print_and_exit("Error in backup file\n!"
                     "Ly should be %u, rather than %u\n\a",Ly,datab.ly);

  data.itcut=datab.itcut; //Next data-block to be generated

  //Now read the data for the PRNG
  fread(&zseed,sizeof(zseed),1,Fconfig);
  //Finally, the configuration
  if(!fread(&u,sizeof(u),1,Fconfig))
    {
      fclose(Fconfig);
      print_and_exit("Error reading backup file\n");
    }
  fclose(Fconfig);  
  printf("backup file read!\n");
}


//This function initialize (or read from backup)  spins, probabilities, etc. 
void Init(void)
{
  int x,y,i;
  char name[1024];
  FILE * Fcheck;

  //First take care of addresing
  for(x=0;x<Lx;x++)
    {
      plus_x[x]=1;

      minus_x[x]=-1;
    }
  //Periodic boundary conditions
  minus_x[0]=Lx-1;
  plus_x[Lx-1]=-minus_x[0];

  for(y=0;y<Ly;y++)
    {
      fc_plus_y[y]=1;
      fc_minus_y[y]=-1;
      plus_y[y]=Lx;
      minus_y[y]=-Lx;
    }
  //Periodic boundary conditions
  fc_minus_y[0]=Ly-1;
  fc_plus_y[Ly-1]=-fc_minus_y[0];
  minus_y[0]=Lx*(Ly-1);
  plus_y[Ly-1]=-minus_y[0];

  //Now table for Fourier Transform
  for(x=0;x<Lx;x++)
    {
      cos_x[x]=cos(2.*M_PI*x/Lx);
      sin_x[x]=sin(2.*M_PI*x/Lx);
    }
  for(y=0;y<Ly;y++)
    {
      cos_y[y]=cos(2.*M_PI*y/Ly);
      sin_y[y]=sin(2.*M_PI*y/Ly);
    }
  
  //Initialize tables for simulation
  Init_Boltz(data.beta);

  //The spin configuration depends on the flag:
  if(data.flag==2)
    {
      read_back_up(); //We read the spin configuration and the state of PRNG
    }
  else
    {
      //If this is a new simulation there should be not a backup file
      sprintf(name,"%s/conf.dat",directory);
      if(NULL==(Fcheck=fopen(name,"rb")))
        { //OK, there is not backup file
          ;
        }
      else
        print_and_exit("File %s already exists!\n\a",name);


      //First, random-number generator
      init_seed();
      if(data.flag==1)
        {
          printf("\nCold start!\n\n");
          for(i=0;i<V;i++)
            u[i]=1;
        }
      else
        {
          printf("\nHot start!\n\n");
          for(i=0;i<V;i++)
            {
              if(MYRANDOM<TWO63)
                u[i]=1;
              else
                u[i]=-1;
            }
        }
    }
}

void Init_Boltz(double beta)
{
  double h,pippo;
  int i;
  for(i=0;i<5;i++) //Table for SD
    {
      h=beta*(-4+2*i);
      pippo=TWO64*(exp(h)/(exp(h)+exp(-h)));
      prob[i]=(unsigned long long) pippo;
    }
  pippo=TWO64*(1.-exp(-2.*beta));
  p=(unsigned long long) pippo;// Conditional Link occupation prob.
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


void check_file_output(int blo)
{

  char name[1024];
  FILE *Fcheck;

  sprintf(name,"%s/OUT%03d.dat",directory,blo);
  if(NULL==(Fcheck=fopen(name,"rb")))
    { //OK, there is not backup file
      ;
    }
  else
    print_and_exit("There already is a conf file!\n\a");

}

// Function prepared by Luis Antonio Fernandez.
// It is used as C printf function. Therefore, the difficult part is to handle
// a variable number of arguments. This is what the va_list is for. 
void print_and_exit(char *format, ...)
{
  va_list list;

  va_start(list,format);
  vprintf(format,list);
  va_end(list);
  exit(1);
}
