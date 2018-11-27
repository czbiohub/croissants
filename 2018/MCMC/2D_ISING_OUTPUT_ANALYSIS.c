//Analysis program for the output of the simulation of the Ising model in 2D
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h> //Needed for the function print_and_exit()

//Global definitions
#define MAXIT 10000  //Maximum-lenght of the data file
#define MAXBETAS 21  //Maximum number of extrapolated beta values
#define MAXBLO  200  //Maximum number of Jack-Knife blocks
#define NBIN    50  //Number of bins in energy histograms

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



//Variables. Global vectors can be assigned a bigger length.
double v_read_e[MAXIT]; //We read here data file. One vector per observable
double v_read_m[MAXIT];
double v_read_f[MAXIT];


//Variables for jack-knife. For each observable we also need correlation 
//with energy density
double expo[MAXBETAS][MAXBLO];   //Boltzmann factor
double ener[MAXBETAS][MAXBLO];  
double ener2[MAXBETAS][MAXBLO];
double M2[MAXBETAS][MAXBLO];
double M2ener[MAXBETAS][MAXBLO];
double M4[MAXBETAS][MAXBLO];
double M4ener[MAXBETAS][MAXBLO];
double F[MAXBETAS][MAXBLO];
double Fener[MAXBETAS][MAXBLO];

//Extrapolated beta values
double v_beta[MAXBETAS];

double volume; //Lattice volume
double geometric; //2 sin(pi/L), needed for correlation-lenght

//Function declarations
void check_compatibility(s_data,s_data,char *);//data files are miscible?
void print_out_sdata(s_data);//Print-out a data structure.
void print_and_exit(char *format, ...); //Print out an error message and exits
void make_sums(int,int,int,int,double,double); //Makes the JK blocks
void calculate_standard(double *,double *,double *,double *,double *,double *,
                        int,int); //Calculate errors for standard obervable
//The following two functions calculate errors for <M^4>/<M^2>^2 and 
//the correlation-length, respectively.
void calculate_M4overM22(double *,double *,double *,double *,int,int);
void calculate_xi(double *,double *,double *,double *,int,int);
//Normalize histogram and calculate error
void calculate_Histog(double *,double *,int,double);

int main(int argc, char **argv)
{
 
  s_data data,datab;
  char name_list[1024], name_input[1024];
  FILE *Flist, *Finput, *Foutput;
  int nfile;        //Total number of files
  int nfilediscard; //number of data file to be discarded for thermalization
  int nblo;          //Number of JK-blocks
  int lblo;          //Lenght (in files!) of data blocks
  int myblo;         //Block to which the just read file  belongs.
  int nbeta;         //Number of extrapolated beta-values
  double dbeta;      //Maximum distance for extrapolation.
  //Variables for analysis
  double e,e2,e_max,e_min;
  int ifile,it,ibeta;
  double Obs[MAXBETAS],DObs[MAXBETAS],error[MAXBETAS],errorD[MAXBETAS];
  
  //To control if  re-weighting is sensible, check the energy histograms:
  double Histog[MAXBLO][NBIN],Histog_Plus[MAXBLO][NBIN];
  double Histog_Minus[MAXBLO][NBIN];
  double errHist[NBIN],errHist_Plus[NBIN],errHist_Minus[NBIN];
  double lbin,e_ini,pippo_Plus,pippo_Minus; 
  int iblo;
  int mybin; //Histogram-bin to which an energy belongs.
  double L;//A silly variable


  if(argc!=7)
    print_and_exit("Usage: %s name_list nfile nfilediscard nblo dbeta nbeta\n",
                   argv[0]);
  //Read input parameters
  sscanf(argv[1],"%s",name_list);
  sscanf(argv[2],"%d",&nfile);
  sscanf(argv[3],"%d",&nfilediscard);
  sscanf(argv[4],"%d",&nblo);
  sscanf(argv[5],"%lf",&dbeta);
  sscanf(argv[6],"%d",&nbeta);
  //Check of compatibility
  
  if(!(nbeta&1)) //We want an odd number of beta values, in order
    nbeta++;     //to have simulation's beta included in the list.
  if(nbeta>MAXBETAS)
    print_and_exit("Compile again with MAXBETAS>%d\n",nbeta);


  if(nblo+1>MAXBLO) //We put the mean value in the block indexed as "nblo"
    {
      printf("Compile again with MAXBLO>%d\n",nblo+1);
    }

  //Adjustments to make (nfile-nfilediscard) divisible by nblo;
  lblo=(nfile-nfilediscard)/nblo;
  if(lblo<1)
    {
      printf("nblo should be smaller than %d\n\a",nfile-nfilediscard);
    }
  nfilediscard=nfile-lblo*nblo;
  printf("\nWe discard %d files, and keep %d block of %d files each\n",
         nfilediscard,nblo, lblo);

  //First reading to calculate average energy,
  //and to check for compatibility of all data
  if(NULL==(Flist=fopen(name_list,"r")))
    print_and_exit("I could not open file %s\n",name_list);

  e=e2=0; //Accumulators are put to zero
  e_max=-2000; //A ridiculous value should be assigned first. In this way 
  e_min=2000;  //we  know that the first read data will be larger than e_max
  //The first time, the structure data is read.
  for(ifile=0;ifile<nfilediscard;ifile++)
    {
      fscanf(Flist,"%s",name_input);
      if(NULL==(Finput=fopen(name_input,"r")))
        print_and_exit("I could not open file %s\n",name_input);
      fclose(Finput);
    }
  

  for(ifile=nfilediscard;ifile<nfile;ifile++)
    {
      fscanf(Flist,"%s",name_input);

      if(ifile==nfilediscard)  //The first time we read a file, we keep the structure data.
        {

          if(NULL==(Finput=fopen(name_input,"r")))
            print_and_exit("I could not open file %s\n",name_input);

          fread(&data,sizeof(data),1,Finput);
          if(data.itmax>MAXIT) //Compatibility check 
            print_and_exit("Recompile with MAXIT>%d\n",data.itmax);

          fclose(Finput);      

          //We can now assign the beta values:
          if(nbeta>1)
            for(ibeta=0;ibeta<nbeta;ibeta++)
              v_beta[ibeta]=data.beta+dbeta*(-1+ibeta/((nbeta-1.)/2.));
          else
            v_beta[0]=data.beta;
          //The volume is also known, by now:
          volume=(double) (data.lx*data.ly);
          //geometric=2 sin (pi/L), L:maximum lattice dimension
          L=data.lx;
          if(data.ly>L)
            L=data.ly;
          geometric=2.*sin(M_PI/L);

        }
      if(NULL==(Finput=fopen(name_input,"r")))
        print_and_exit("I could not open file %s\n",name_input);
      fread(&datab,sizeof(data),1,Finput);      
      check_compatibility(data,datab,name_input);
      if(data.itmax!=fread(&v_read_e,sizeof(double),data.itmax,Finput))
        print_and_exit("Error reading energies from %s\n",name_input);
      fclose(Finput);      
      for(it=0;it<data.itmax;it++)
        {
          e+=v_read_e[it];
          e2+=v_read_e[it]*v_read_e[it];
          e_max=(v_read_e[it]>e_max)?v_read_e[it]:e_max;
          e_min=(v_read_e[it]<e_min)?v_read_e[it]:e_min;
        }
    }
  fclose(Flist);  
  e/=data.itmax*(nfile-nfilediscard);
  e2/=data.itmax*(nfile-nfilediscard);
  printf("\n<e>=%.8g, var_e=%.8g, e_max=%.8g e_min=%.8g\n",
         e,sqrt(e2-e*e),e_max,e_min);
  printf("Maximum safe extrapolation: 2/(vol*sqrt(var_e))=%.8g\n",
         2./(volume*sqrt(e2-e*e)));
  printf("You are asking for %.4g times that much.\n\n",
         dbeta*volume*sqrt(e2-e*e)/2);

  //Setting variables for histogramms
  lbin=(e_max-e_min)/(NBIN-1.); //The size of the histogram-bins
  e_ini=e_min-0.5*lbin; //First bin covers  e_min-0.5*lbin < e < e_min+0.5*lbin
                        //Last bin goes from e_max-0.5*lbin to e_max+0.5*lbin
  printf("\nCheck: e_max is assigned bin %d, which should equal %d\n\n",
         (int)((e_max-e_ini)/lbin),NBIN-1);

  //Starting Jack-Knife! We first put accumulators to zero:
  for(ibeta=0;ibeta<nbeta;ibeta++)
    for(myblo=0;myblo<=nblo;myblo++)
      {
        expo[ibeta][myblo]=0;  
        ener[ibeta][myblo]=0;  
        ener2[ibeta][myblo]=0;
        M2[ibeta][myblo]=0;
        M2ener[ibeta][myblo]=0;
        M4[ibeta][myblo]=0;
        M4ener[ibeta][myblo]=0;
        F[ibeta][myblo]=0;
        Fener[ibeta][myblo]=0;
      }
  for(myblo=0;myblo<=nblo;myblo++)
    for(mybin=0;mybin<NBIN;mybin++)
      {
        Histog[myblo][mybin]=0;
        Histog_Plus[myblo][mybin]=0;
        Histog_Minus[myblo][mybin]=0;
      }

  //Now we read again, forming Jack-Knife blocks
  if(NULL==(Flist=fopen(name_list,"r")))
     print_and_exit("I could not open file %s\n",name_list);

  for(ifile=0;ifile<nfilediscard;ifile++)
    fscanf(Flist,"%s",name_input);


  for(ifile=nfilediscard;ifile<nfile;ifile++)
    {
      
      myblo=(ifile-nfilediscard)/lblo; //This file belongs to block "myblo"
      if((ifile-nfilediscard)%(10*lblo)==0)
        printf("File number %d. Now treating block number %d\n",ifile,myblo);

      fscanf(Flist,"%s",name_input);
      if(NULL==(Finput=fopen(name_input,"r")))
        print_and_exit("I could not open file %s\n",name_input);

      fread(&datab,sizeof(data),1,Finput);      
      if(data.itmax!=fread(&v_read_e,sizeof(double),data.itmax,Finput))
        print_and_exit("Error reading energies from %s\n",name_input);

      if(data.itmax!=fread(&v_read_m,sizeof(double),data.itmax,Finput))
        print_and_exit("Error reading magnetizations from %s\n",name_input);


      if(data.itmax!=fread(&v_read_f,sizeof(double),data.itmax,Finput))
        print_and_exit("Error reading f from %s\n",name_input);

      fclose(Finput); 
      //Add to corresponding JK-blocks
      make_sums(nblo,myblo,data.itmax,nbeta,data.beta,e);

      //Now the histogram. In the reweighting we multiply by
      //exp(Vdbeta(e-<e>) rather than exp(Vdbeta) to avoid an overflow
      for(it=0;it<data.itmax;it++)
        {       
          mybin=(int)((v_read_e[it]-e_ini)/lbin);
          pippo_Plus=exp(volume*dbeta*(v_read_e[it]-e));
          pippo_Minus=exp(-volume*dbeta*(v_read_e[it]-e));;
          for(iblo=0;iblo<=nblo;iblo++)
            if(iblo!=myblo)
              {
                Histog[iblo][mybin]++;
                Histog_Plus[iblo][mybin]+=pippo_Plus;
                Histog_Minus[iblo][mybin]+=pippo_Minus;
              }
        }
    }
  fclose(Flist);  

  //Calculate normalizating factors for histogram, obtain errors and 
  //print-out to disk

  calculate_Histog(&Histog[0][0],errHist,nblo,lbin);
  calculate_Histog(&Histog_Plus[0][0],errHist_Plus,nblo,lbin);
  calculate_Histog(&Histog_Minus[0][0],errHist_Minus,nblo,lbin);
   Foutput=fopen("Histogram.dat","wt");
  for(mybin=0;mybin<NBIN-1;mybin++)
    fprintf(Foutput,"%.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",e_min+mybin*lbin,
            Histog[nblo][mybin],errHist[mybin],
            Histog_Plus[nblo][mybin],errHist_Plus[mybin],
            Histog_Minus[nblo][mybin],errHist_Minus[mybin]);
  fclose(Foutput);

  //Print out averages and errors
  //First the observables and their derivatives
  //The energy:
  calculate_standard(Obs,DObs,error,errorD,&ener[0][0],&ener2[0][0],nbeta,nblo);
  Foutput=fopen("Energy.dat","wt");
  for(ibeta=0;ibeta<nbeta;ibeta++)
    fprintf(Foutput,"%.8g %.8g %.8g %.8g %.8g\n",
            v_beta[ibeta],Obs[ibeta],error[ibeta],DObs[ibeta],errorD[ibeta]);
  fclose(Foutput);

  //The magnetization:
  calculate_standard(Obs,DObs,error,errorD,&M2[0][0],&M2ener[0][0],nbeta,nblo);
  Foutput=fopen("M2.dat","wt");
  for(ibeta=0;ibeta<nbeta;ibeta++)
    fprintf(Foutput,"%.8g %.8g %.8g %.8g %.8g\n",
            v_beta[ibeta],Obs[ibeta],error[ibeta],DObs[ibeta],errorD[ibeta]);
  fclose(Foutput);

  //Now, adimensional functions of observables:
  // 1) Binder-cumulant like <M^4>/<M^2>^2:
  calculate_M4overM22(Obs,DObs,error,errorD,nbeta,nblo);
  Foutput=fopen("M4overM22.dat","wt");
  for(ibeta=0;ibeta<nbeta;ibeta++)
    fprintf(Foutput,"%.8g %.8g %.8g %.8g %.8g\n",
            v_beta[ibeta],Obs[ibeta],error[ibeta],DObs[ibeta],errorD[ibeta]);
  fclose(Foutput);

  // 2) Correlation-length
  calculate_xi(Obs,DObs,error,errorD,nbeta,nblo);
  Foutput=fopen("xi.dat","wt");
  for(ibeta=0;ibeta<nbeta;ibeta++)
    fprintf(Foutput,"%.8g %.8g %.8g %.8g %.8g\n",
            v_beta[ibeta],Obs[ibeta],error[ibeta],DObs[ibeta],errorD[ibeta]);
  fclose(Foutput);

  return 0;
}


void check_compatibility(s_data dataA,s_data dataB,char * name)
{
  int error;
  error=0;
  if(dataA.itmax!=dataB.itmax)
    error=1;
  if(dataA.mesfr!=dataB.mesfr)
    error=1;
  if((dataA.lx!=dataB.lx)||(dataA.ly!=dataB.ly))
    error=1;
  if(fabs(dataA.beta-dataB.beta)>1e-12)
    error=1;

  if(error)
    {
      printf("Incompatibility in file %s\n",name);
      printf("First data structure read:\n");
      print_out_sdata(dataA);
      printf("Second data structure read:\n");
      print_out_sdata(dataB);
      exit(1);
    }
}

void print_out_sdata(s_data mydata)
{
  printf("itmax  %u \n",mydata.itmax);
  printf("mesfr  %u \n",mydata.mesfr);
  printf("nbin   %u \n",mydata.nbin);
  printf("itcut  %u \n",mydata.itcut);
  printf("nterm  %u \n",mydata.nterm);
  printf("flag   %u \n",mydata.flag);
  printf("lx     %u \n",mydata.lx);
  printf("ly     %u \n",mydata.ly);
  printf("beta   %lf\n",mydata.beta);
}
 
//Makes the JK blocks. First make the sums, then  distribute them.
void make_sums(int nblo,int iblo,int maxit,int nbeta,double beta,double e_mean)
{
  int it,jblo,ibeta;
  double myener,expo_factor,m2,m4,f;
  double v_expo[MAXBETAS],v_e[MAXBETAS],v_ee[MAXBETAS],v_m2[MAXBETAS];
  double v_m2e[MAXBETAS],v_m4[MAXBETAS],v_m4e[MAXBETAS],v_f[MAXBETAS];
  double v_fe[MAXBETAS];

  //First put to zero accumulators
  for(ibeta=0;ibeta<nbeta;ibeta++)
    {
      v_expo[ibeta]=0; v_e[ibeta]=0; v_ee[ibeta]=0;v_m2[ibeta]=0;
      v_m2e[ibeta]=0; v_m4[ibeta]=0; v_m4e[ibeta]=0; v_f[ibeta]=0;
      v_fe[ibeta]=0;
    }
  for(it=0;it<maxit;it++)
    {
      myener=v_read_e[it];
      m2=v_read_m[it]*v_read_m[it];
      m4=m2*m2;
      f=v_read_f[it];
      if(nbeta==1)
        {
          v_expo[0]+=1;
          v_e[0]+=myener;
          v_ee[0]+=myener*myener;
          v_m2[0]+=m2;
          v_m2e[0]+=m2*myener;
          v_m4[0]+=m4;
          v_m4e[0]+=m4*myener;
          v_f[0]+=f;
          v_fe[0]+=f*myener;
        }
      else
        {
          for(ibeta=0;ibeta<nbeta;ibeta++)
            {
              expo_factor=exp(volume*(v_beta[ibeta]-beta)*(myener-e_mean));
              v_expo[ibeta]+=expo_factor;              
              v_e[ibeta]+=expo_factor*myener;          
              v_ee[ibeta]+=expo_factor*myener*myener;      
              v_m2[ibeta]+=expo_factor*m2;             
              v_m2e[ibeta]+=expo_factor*m2*myener;       
              v_m4[ibeta]+=expo_factor*m4;             
              v_m4e[ibeta]+=expo_factor*m4*myener;       
              v_f[ibeta]+=expo_factor*f;               
              v_fe[ibeta]+=expo_factor*f*myener;         
            }
        }
    }
  //We are ready for meaking the blocks:
  for(ibeta=0;ibeta<nbeta;ibeta++)
    for(jblo=0;jblo<=nblo;jblo++)
      {
        if(jblo!=iblo)
          {//iblo never is nblo, hence there it is the average of all data
            expo[ibeta][jblo]+=v_expo[ibeta];
            ener[ibeta][jblo]+=v_e[ibeta];
            ener2[ibeta][jblo]+=v_ee[ibeta];
            M2[ibeta][jblo]+=v_m2[ibeta];
            M2ener[ibeta][jblo]+=v_m2e[ibeta];
            M4[ibeta][jblo]+=v_m4[ibeta];
            M4ener[ibeta][jblo]+=v_m4e[ibeta];
            F[ibeta][jblo]+=v_f[ibeta];
            Fener[ibeta][jblo]+=v_fe[ibeta];
          }
      }
}

//Gives out the vector obs containing the average of the calculated
//operator for each beta, and DObs wich contains the beta derivative
//Vectors  error and errorD are error-estimates for both
//data is indexed as data[ibeta][iblo]=data[ibeta*MAXBLO+iblo]
//We need to pass data and dataE to make this work for a generic
//observable.
void calculate_standard(double *Obs,double *DObs,double *error,
                        double *errorD,double *data,double *dataE,
                        int nbeta,int nblo)
{
  int iblo,ibeta;
  double sum,sum2,sumderiv,sumderiv2;
  double *vdat,*vdate;
  double temp,tempD,tempe;

  for(ibeta=0;ibeta<nbeta;ibeta++)
    {
      sum=sum2=sumderiv=sumderiv2=0;
      vdat=data+ibeta*MAXBLO; //vdat points to the first JK-block at ibeta
      vdate=dataE+ibeta*MAXBLO;

      temp=vdat[nblo]/expo[ibeta][nblo];
      tempe=ener[ibeta][nblo]/expo[ibeta][nblo];

      Obs[ibeta]=temp;
      DObs[ibeta]=volume*(vdate[nblo]/expo[ibeta][nblo]-temp*tempe);
      for(iblo=0;iblo<nblo;iblo++)
        {
          temp=vdat[iblo]/expo[ibeta][iblo];
          tempe=ener[ibeta][iblo]/expo[ibeta][iblo];
          tempD=volume*(vdate[iblo]/expo[ibeta][iblo]-temp*tempe);
          sum+=temp;
          sum2+=temp*temp;
          sumderiv+=tempD;
          sumderiv2+=tempD*tempD;
        }
      sum/=nblo;
      sum2/=nblo;
      sumderiv/=nblo;
      sumderiv2/=nblo;
      error[ibeta]=sqrt((nblo-1.)*(sum2-sum*sum));
      errorD[ibeta]=sqrt((nblo-1.)*(sumderiv2-sumderiv*sumderiv));
    }
}

//As the previous function, but only  for <M^4>/<M^2>^2
//To keep formulae simple, we use logarithmic derivation
void  calculate_M4overM22(double *Obs,double *DObs,double *error,
                        double *errorD,int nbeta,int nblo)
{
  int iblo,ibeta;
  double sum,sum2,sumderiv,sumderiv2;
  double temp,tempD,tempe,derivM2,derivM4,m4,m2,m4e,m2e;

  for(ibeta=0;ibeta<nbeta;ibeta++)
    {
      sum=sum2=sumderiv=sumderiv2=0;
      m4=M4[ibeta][nblo]/expo[ibeta][nblo];
      m4e=M4ener[ibeta][nblo]/expo[ibeta][nblo];
      m2=M2[ibeta][nblo]/expo[ibeta][nblo];
      m2e=M2ener[ibeta][nblo]/expo[ibeta][nblo];
      tempe=ener[ibeta][nblo]/expo[ibeta][nblo];
      derivM4=volume*(m4e-m4*tempe);
      derivM2=volume*(m2e-m2*tempe);
      temp=m4/(m2*m2);
      tempD=temp*(derivM4/m4 - 2.*derivM2/m2);
      Obs[ibeta]=temp;
      DObs[ibeta]=tempD;
      for(iblo=0;iblo<nblo;iblo++)
        {
          m4=M4[ibeta][iblo]/expo[ibeta][iblo];
          m4e=M4ener[ibeta][iblo]/expo[ibeta][iblo];
          m2=M2[ibeta][iblo]/expo[ibeta][iblo];
          m2e=M2ener[ibeta][iblo]/expo[ibeta][iblo];
          tempe=ener[ibeta][iblo]/expo[ibeta][iblo];
          derivM4=volume*(m4e-m4*tempe);
          derivM2=volume*(m2e-m2*tempe);
          temp=m4/(m2*m2);
          tempD=temp*(derivM4/m4 - 2*derivM2/m2);
          sum+=temp;
          sum2+=temp*temp;
          sumderiv+=tempD;
          sumderiv2+=tempD*tempD;
        }
      sum/=nblo;
      sum2/=nblo;
      sumderiv/=nblo;
      sumderiv2/=nblo;
      error[ibeta]=sqrt((nblo-1.)*(sum2-sum*sum));
      errorD[ibeta]=sqrt((nblo-1.)*(sumderiv2-sumderiv*sumderiv));
    }
}


//As the previous function, for correlation-length
//To keep formulae simple, we use logarithmic derivation
void  calculate_xi(double *Obs,double *DObs,double *error,
                        double *errorD,int nbeta,int nblo)
{
  int iblo,ibeta;
  double sum,sum2,sumderiv,sumderiv2;
  double temp,tempD,tempe,derivM2,derivF,f,m2,fe,m2e;

  for(ibeta=0;ibeta<nbeta;ibeta++)
    {
      sum=sum2=sumderiv=sumderiv2=0;
      f=F[ibeta][nblo]/expo[ibeta][nblo];
      fe=Fener[ibeta][nblo]/expo[ibeta][nblo];
      m2=M2[ibeta][nblo]/expo[ibeta][nblo];
      m2e=M2ener[ibeta][nblo]/expo[ibeta][nblo];
      tempe=ener[ibeta][nblo]/expo[ibeta][nblo];
      derivF=volume*(fe-f*tempe);
      derivM2=volume*(m2e-m2*tempe);

      temp=sqrt(m2/f-1.)/geometric;
      tempD=0.5*temp*((derivM2-derivF)/(m2-f) - derivF/f);
      Obs[ibeta]=temp;
      DObs[ibeta]=tempD;
      for(iblo=0;iblo<nblo;iblo++)
        {

          f=F[ibeta][iblo]/expo[ibeta][iblo];
          fe=Fener[ibeta][iblo]/expo[ibeta][iblo];
          m2=M2[ibeta][iblo]/expo[ibeta][iblo];
          m2e=M2ener[ibeta][iblo]/expo[ibeta][iblo];
          tempe=ener[ibeta][iblo]/expo[ibeta][iblo];
          derivF=volume*(fe-f*tempe);
          derivM2=volume*(m2e-m2*tempe);
          
          temp=sqrt(m2/f-1.)/geometric;
          tempD=0.5*temp*((derivM2-derivF)/(m2-f) - derivF/f);

          sum+=temp;
          sum2+=temp*temp;
          sumderiv+=tempD;
          sumderiv2+=tempD*tempD;
        }
      sum/=nblo;
      sum2/=nblo;
      sumderiv/=nblo;
      sumderiv2/=nblo;
      error[ibeta]=sqrt((nblo-1.)*(sum2-sum*sum));
      errorD[ibeta]=sqrt((nblo-1.)*(sumderiv2-sumderiv*sumderiv));
    }
}

//This function normalizes the histogram and calculates errors
void calculate_Histog(double *Hist,double *err,int nblo,double h)
{
  double *vector;
  double sum,sum2,norma,pippo;
  int iblo,ibin;
  //We first normalize the Histogram (for each block!)
  for(iblo=0;iblo<=nblo;iblo++)
    {
      vector=Hist+iblo*NBIN;//vector points now to the data for block iblo
      norma=0;
      for(ibin=0;ibin<NBIN;ibin++)
        norma+=vector[ibin];
      for(ibin=0;ibin<NBIN;ibin++)
        vector[ibin]/=norma*h;
    }
  //Now calculate errors!
  for(ibin=0;ibin<NBIN;ibin++)
    {
      sum=sum2=0;
      for(iblo=0;iblo<nblo;iblo++)
        {
          pippo=Hist[iblo*NBIN+ibin];
          sum+=pippo;
          sum2+=pippo*pippo;
        }
      sum2/=nblo;
      sum/=nblo;
      err[ibin]=sqrt((nblo-1.)*(sum2-sum*sum));
    }
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

