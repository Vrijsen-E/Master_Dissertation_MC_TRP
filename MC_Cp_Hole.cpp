/** Elien Vrijsen | Vrijsen.E@gmail.com

-----------------------------
Before compiling and running:
(settings -> compiler -> compiler flags)
    - set the -std=c++11
    - optimization flag -O3
**/

#include <iostream> //use code from ...
#include <cmath> //math functions
#include <cstdlib> //rand()
#include <utility> //function return more than one value
#include <stdio.h> //for array as input in function/N before main
#include <random> //random numbers better than rand()
#include <time.h> //for random numbers
#include <stdlib.h>
#include <typeinfo> //to check type of variables
#include <chrono> //time in terminal
#include <fstream>
#include <sstream> //make string of double
#include <string>

using namespace std;
using namespace std::chrono; //time in terminal


int main(){
/// timer in terminal
auto timer_start  = std::chrono::high_resolution_clock::now();
auto timer_finish = std::chrono::high_resolution_clock::now();
timer_start = std::chrono::high_resolution_clock::now();
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
        /*int rand_01; (?)*/


/**------------------------------------
        DECLARATION OF VARIABLES
    consider a separate file for this
------------------------------------**/
double DH; double Tm; double R;
int Dwgl; /// !! Dwgl is an integer
std::vector<int> Dwgl_s; //If you know you're going to need a collection of variables but not precisely how many of them
int length_Dwgl_s;
int N; int Npro; int MCcycles;
int MC0; double wH;
long long int counter_non_acc=0;
long long int counter_k=0;
long long int counter_l=0;
/// counter_DE0: when DE=<0
long long int counter_DE0=0;
long long int counter_acc=0;
int ngl_s[7]={0,0,0,0,0,0,0};
int ab; int be; int le; int rig;
int DE_index; int sij;
double DE; double expDE;
double r; //random number
//int MCsteps=N*N*MCcycles;
int length_Ts;
double tN=0; double tNsq=0; double tc=0;

///quantities used in loops
double DEs[14]; double ExponentDE[14];
double a; double b;
int i; int j; int Nl; int Ngl;
double T; double beta;

int i_up; int j_up; int add_j;

/**------------------------------------
            c++ variables
------------------------------------**/
int nxy[16]={};//all initialized to zero
int ngg; int ngc; int ngl; int ngo;
int nll; int nlc; int nlo;
int ncc; int nco; int noo;
int s_ab; int s_be; int s_le;
int s_ri; int s_ar; int s_bl;
int Nlip;
int nxx_s[4][7][6][7][6][7][6][6][6][7][7];

double allDEs[4][7][7][7][7];


/**------------------------------------
            FIXED PARAMETERS
------------------------------------**/
///gas constant
R=8.3144598; //J/(K*mol)


/**------------------------------------
---------------------------------------
            SET PARAMETERS
---------------------------------------
------------------------------------**/
///lattice size
N=65;
///number of protein spots along one axis (square: Npro x Npro)
Npro=14;


///number of Monte Carlo cycles to average over for heat capacity Cp
MCcycles=60000;
///MC cycles before equilibration
MC0=1000;


DH=36400; //J/mol
wH=0; //J/mol


///melting temperature
Tm=310.3; //DPPC melting temp, 41°C //K 37.15°C=310.3K / 41.3°C=314.45K

///evaluated temperatures
double Ts[]={Tm-11,Tm-9,Tm-7.5,Tm-6,Tm-5,Tm-4,Tm-3,Tm-2,Tm-1.5,Tm-1,Tm-0.75,Tm-0.5,Tm-0.25,Tm-0.15,Tm,Tm+0.15,Tm+0.25,Tm+0.5,Tm+0.75,Tm+1,Tm+1.5,Tm+2,Tm+3,Tm+4,Tm+5,Tm+6,Tm+7.5,Tm+9,Tm+11};



///enter Dwgl through file?
bool Dwgl_input_file=false;

///OPTIONAL: set cooperativity parameters (liquid/gel) manually here
if (!Dwgl_input_file){
    Dwgl_s={1381,1255,1297,1339,1381};//J/mol
}


///OPTIONAL: put calor=true for cal instead of in Joules
bool calor=false;



/**------------------------------------
            IMPORT PARAMETERS
------------------------------------**/
/**read in cooperativity parameters (liquid/gel) from file
    Dwgl_requested can be created through input in the terminal**/

if (Dwgl_input_file){
    /// read in values for Dwgl from file
    ifstream inFile; //input file stream (ifstream) variable
    inFile.open("Dwgl_requested.txt");

    /// error if the file Dwgl_requested cannot be opened
    if (!inFile){
        cerr << "Unable to open file Dwgl_Requested.txt";
        exit(1); /// call system to stop
    }

    double x;
    int numVal=1;
    int read_ind=0;
    while (inFile >> x) {
        Dwgl_s.resize(numVal);
        Dwgl_s[read_ind]=x;
        read_ind+=1; numVal+=1;
    }
    inFile.close();
}


/**------------------------------------
            dependent variables
------------------------------------**/
Nlip=N^2-Npro^2;

///number of different Dwgl and temperature values that will be assessed
length_Dwgl_s=sizeof(Dwgl_s)/sizeof(Dwgl_s[0]);
length_Ts=sizeof(Ts)/sizeof(Ts[0]);

///quantities for storage
int Ngl_s[MCcycles]; int Nl_s[MCcycles];

///variables for Cp calculation
double var_N[length_Ts];
double std_N[length_Ts];
double var_N_N2[length_Ts];
double Cp[length_Ts];
int Nsq=N*N; int s[N][N];

///lower bound and upper bound of protein indices
int l_bound=round((N-Npro)/2); //these are the lowest and highest indices at which the 'ion channel' is
int u_bound=l_bound+Npro-1; //N=65,Npro=14 -> l_bound=25

/// change to calories if that is requested
if(calor){
    DH=DH/4.1868; //cal/mol
    Dwgl=Dwgl/4.1868; //cal/mol
    R=R/4.1868; //cal/(K mol)
    wH=wH/4.1868; //cal/mol
}

/**------------------------------------
---------------------------------------
            LOOP OVER Dwgl
---------------------------------------
------------------------------------**/

for (int coop=0;coop<length_Dwgl_s;coop++){
    Dwgl=Dwgl_s[coop];
    ///don't do anything if Dwgl==0 (?)
    if (!Dwgl){}
    else{
        cout<<"--------------------------------------------"<<endl;
        cout<<"--------------------------------------------"<<endl;
        cout<<"            Dwgl = "<<Dwgl<<endl;
        cout<<"--------------------------------------------"<<endl;
        cout<<"--------------------------------------------"<<endl;



        /**------------------------------------
                    LOOP OVER Ts
        ------------------------------------**/
        for(int y=0;y<length_Ts;y++){
            T=Ts[y];
            beta=1/(R*T);

            ///calculate exponentials
            for(int v=0;v<=6;v++){ //0->6, <7->13
                a=DH*(1-T/Tm);
                b=(6-2*v)*Dwgl;
                DEs[v]=a+b;
                ExponentDE[v]=exp(-beta*DEs[v]);
                DEs[v+7]=-a+b;
                ExponentDE[v+7]=exp(-beta*DEs[v+7]);
            }

            for(int sel=-1;sel<=2;sel++){
                    for(int nng=0;nng<=6;nng++){
                        for(int nnc=0;nnc<=6;nnc++){
                            for(int nnl=0;nnl<=6;nnl++){
                                for(int nno=0;nno<=6;nno++){
                                        if(abs(sel)==1){
                                            double DG=DH*(1-T/Tm);
                                            double DHO = 0;
                                            double DSO = 0;
                                            double DHC = 0;
                                            double DSC = 0;
                                            double DGO=DHO-T*DSO;
                                            double DGC=DHO-T*DSO;
                                            allDEs[sel+1][nng][nnc][nnl][nno]=-sel*(DG+(nng-nnl)*Dwgl+nno*DGO+nnc*DGC);
                                        }
                                }
                            }
                        }
                    }
            }


            ///initialize: all lipids are in the gel phase: -1 (liquid=+1)
            /// all are set to -1
            for(int ir=0;ir<=N-1;ir++){ ///indices range from 0 to N-1
                for(int jr=0;jr<=N-1;jr++){s[ir][jr]=-1;
                }
            }

            ///initialize: all lattice sites inside of 'ion channel' are 'closed'=0 (open=2)
            for(int ir=l_bound;ir<=u_bound;ir++){
                for(int jr=l_bound;jr<=u_bound;jr++){s[ir][jr]=0;
                }
            }

            ///Ngl = total of gel-liquid neighbors in the lattice
            Nl=0; Ngl=0;
            Ngl_s[0]=Ngl; Nl_s[0]=Nl;
            tc=0;tN=0;tNsq=0;

            counter_k=0;

            int nr=1;

            /**------------------------------------
                MONTE CARLO SIMULATION
                    for temperature T
            ------------------------------------**/
            for(int k=1;k<=MCcycles;k++){
                counter_k+=1;
                counter_l=0;

                /**------------------------------------
                    ONE MONTE CARLO CYCLE
                ------------------------------------**/
                for(int ls=0;ls<Nsq;ls++){
                    counter_l+=1;
                    ///randomly select the indices of one lattice site
                    i=uniform_int_distribution<int>(0,N-1)(rng);
                    j=uniform_int_distribution<int>(0,N-1)(rng);//0.35s per rand()
                    //in [0,N-1]


                    ///Do not select points 'inside the ion channel'
                    if(!(((l_bound<=i)&&(i<=u_bound))&&((l_bound<=j)&&(j<=u_bound)))){

                        ///reset all n-- to zero
                        for(int i_s=0;i_s<16;i_s++){
                            nxy[i_s]=0;
                        }
                        int nn[4]={};
                        for(int i_s=0;i_s<4;i_s++){
                            nn[i_s]=0;
                        }

                        ///the selected point
                        sij=s[i][j];

                        ///above
                        if(i==N-1){ab=0;} else{ab=i+1;} //indices start counting from zero
                        s_ab=s[ab][j];
                        nn[s_ab+1]++;
                        nxy[(s_ab+1)+(sij+1)*4]++;


                        if(i==0){be=N-1;} else{be=i-1;}
                        s_be=s[be][j];
                        nn[s_be+1]++;
                        nxy[(s_be+1)+(sij+1)*4]++;

                        if(j==0){le=N-1;} else{le=j-1;}
                        s_le=s[i][le];
                        nn[s_le+1]++;
                        nxy[(s_le+1)+(sij+1)*4]++;

                        if(j==N-1){rig=0;} else{rig=j+1;}
                        s_ri=s[i][rig];
                        nn[s_ri+1]++;
                        nxy[(s_ri+1)+(sij+1)*4]++;

                        s_ar=s[ab][rig];
                        nn[s_ar+1]++;
                        nxy[(s_ar+1)+(sij+1)*4]++;

                        s_bl=s[be][le];
                        nn[s_bl+1]++;
                        nxy[(s_bl+1)+(sij+1)*4]++;

                        int nng=nn[0];
                        int nnc=nn[1];
                        int nnl=nn[2];
                        int nno=nn[3];


                        ngg=nxy[0];
                        ngc=nxy[1];
                        ngl=nxy[2];
                        ngo=nxy[3];

                        ngc+=nxy[4];
                        ncc=nxy[5];
                        nlc=nxy[6];
                        nco=nxy[7];

                        ngl+=nxy[8];
                        nlc+=nxy[9];
                        nll=nxy[10];
                        nlo=nxy[11];

                        ngo+=nxy[12];
                        nco+=nxy[13];
                        nlo+=nxy[14];
                        noo=nxy[15];



                        DE=allDEs[sij+1][nng][nnc][nnl][nno];

                        if(DE<=0){s[i][j]=(-sij);//accept, flip
                            Nl-=sij;
                            Ngl+=(6-2*ngl);
                            counter_DE0+=1;
                        }

                        else{expDE=ExponentDE[DE_index];
                            r=uniform_real_distribution<double>(0,1)(rng);
                            //in [0,1)
                            //cout<<"r="<<r<<endl;
                            if(r<=expDE){
                                s[i][j]=(-sij);
                                Nl-=sij;
                                Ngl+=(6-2*ngl);
                                counter_acc+=1;
                            }//3.4s (6000)
                            else{counter_non_acc+=1;}
                        }

                        ///lipid-ion channel

                    }

                    ///the edge of the 'ion channel'
                    else if(((i==l_bound||i==u_bound)&&(l_bound<=j<=u_bound)) || ((j==l_bound||j==u_bound)&&(l_bound<=i<=u_bound) )) {

                    }

                } ///end ls-loop = end of one cycle

               if(k>=MC0 && nr<=60){
                    tc++; tN=tN+Nl;
                    tNsq=tNsq+Nlip*Nlip;

                    ///Output file stream
                    ofstream saveS;
                    std::ostringstream oss;
                    oss << "S/s_T" << round(T) << "_nr_" << nr << ".txt";
                    std::string var = oss.str();
                    ///std::string var = std::string("S/s_T_") + T + std::string("_nr_") + nr + std::string(".txt");
                    saveS.open(var);
                    for(int i=0;i<=N-1;i++){
                            for(int j=0;j<=N-1;j++){
                                /// all N of the first row, followed by all N of the second row...
                                saveS<<s[i][j]<<" ";
                            }
                    }
                    saveS.close();
                    nr+=1;
                }

            }///end for(k=1->MCcycles)


            var_N[y]=tNsq*1./(tc-1)-(tN*1./(tc-1))*(tN*1./tc);

            cout<<"Temperature: T = "<<T<<endl;
            cout<<"--------------------- "<<endl;
            cout<<"varN = "<<var_N[y]<<endl;

            var_N_N2[y]=var_N[y]/Nlip/Nlip;
            cout<<"varN/N^2 = "<<var_N[y]/(N*N-Npro*Npro)<<endl;

            std_N[y]=sqrt(var_N[y]);
            cout<<"std N = "<<sqrt(var_N[y])<<endl;
            //cout<<"tc="<<tc<<endl;
            //double Cp2;

            Cp[y]=DH*DH*var_N[y]/Nlip/Nlip/T/T/R;
            cout<<"Cp = "<<DH*DH*var_N[y]/(N*N-Npro*Npro)/T/T/R<<endl;
            cout<<"--------------------------------------------"<<endl;


        }///end for(y=1->length(Ts))


        /// timer in terminal
        timer_finish = std::chrono::high_resolution_clock::now();
        chrono::duration<double> timer_elapsed       = timer_finish - timer_start;
        cout << "--- Elapsed time: --- \n" ;
        cout << " total:       " << timer_elapsed.count() << " s\n";
        cout <<" "<<endl;
        cout <<"-------------------------"<<endl;
        cout<<"for the last temperature"<<endl;
        cout <<"-------------------------"<<endl;
        cout<<"in k loop = "<<counter_k<<endl;
        cout<<"in l loop = "<<counter_l<<endl;
        cout <<"--------"<<endl;
        cout<<"in DE<=0 => accepted = "<<counter_DE0<<endl;
        cout<<"in DE!=0, but accepted = "<<counter_acc<<endl;
        cout<<"not accepted = "<<counter_non_acc<<endl;


        cout<<"For the last temperature:"<<counter_non_acc<<endl;
        ofstream plotNl; ///Output file stream
        plotNl.open("Nl_s.txt");
        for(int i=0;i<=MCcycles;i++){
            plotNl<<i<<" "<<Nl_s[i]<<"\n";
        }
        plotNl.close();

        ofstream plotNgl;
        plotNgl.open("Ngl_s.txt");
        for(int i=0;i<=MCcycles;i++){
            plotNgl<<i<<" "<<Ngl_s[i]<<"\n";
        }
        plotNgl.close();


        ofstream sdata;
        sdata.open("s2_data.txt");
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                sdata<<s[i][j]<<"\n";
            }
        }
        sdata.close();


        std::ostringstream strs;
        strs << Dwgl;
        std::string Dwgl_str = strs.str();
        string title="Cp_data_new"+Dwgl_str+".txt";

        ofstream Cpdata;
        Cpdata.open(title);
        for(int i=0;i<length_Ts;i++){
            Cpdata<<Ts[i]<<" "<<Cp[i]<<endl;
        }
        Cpdata.close();

    }///end of else (if Dwgl!=0)

}///end loop over Dwgl_s

return 0;
}///end main

