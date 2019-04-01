#ifndef ISING_H
#define ISING_H

#include <iostream>
#include <cmath> //exp()
#include <random>  // class my_random{};
#include <numeric> // accumulate();
using namespace std;
class my_random{
    private:
    std::random_device rd;
    std::mt19937 mt;
    std::uniform_int_distribution<int> dis_bool;
    std::uniform_real_distribution<double> dis_01;
    public:
    my_random():mt(rd()), dis_bool(0,1),dis_01(0.0,1.0) {}
    my_random(int seed):mt(seed), dis_bool(0,1),dis_01(0.0,1.0){}

    bool rd_bool(){
        return dis_bool(mt);
    }
    double rd_01(){
        return dis_01(mt);
    }

    void rd_thermo(bool *s,int Ntotal){
        for(int i=0;i<Ntotal;i++)
            s[i]=rd_bool();
    }

    int rd_getone(int Ntotal){
        return Ntotal*rd_01();
    }
};
my_random rng(2);//global object rng

class my_random_neighbour{
    public:
    int M; //numbers to neighbhour to be selected,  M>=1  M is integer
    double **Cum; //Cum[M][M+2]
    my_random_neighbour(vector<double> Couplings, int M_input){
        M=M_input;
        Cum=new double*[M];
        for(int i=0;i<M;i++){
            Cum[i]=new double[M+2];
        }
        for(int i=0;i<M;i++){
            Cum[i][i+1]=Couplings[i];
            Cum[i][M+1]=1.0;
            for(int j=0;j<i+1;j++){
                Cum[i][j]=0.0;
            }
        }
        for(int i=0;i<M-1;i++){
            for(int j=i+2;j<M+1;j++){
                Cum[i][j]=Cum[i][j-1]+Couplings[j-1];
            }
        }
        for(int i=0;i<M;i++){
            for(int j=i+1;j<M+1;j++){
                Cum[i][j]=1.0-exp(-2.0*Cum[i][j]);
            }
        }
    }
    ~my_random_neighbour(){
        for(int i=0;i<M;i++){
            delete [] Cum[i];
        }
        delete [] Cum;
    }

    vector<int> rd_neighbour(){
        vector<int> output;
        int index=0;
        for(int i=0;i<M; ){
            /*
            index =insert_rand( Cum[i],rng.rd_01(),M+2);
            if(index<M) {output.push_back(index);}//cout<<"push "<<index<<endl;}
            i=index+1;
           */
            index =insert_rand( Cum[i]+i,rng.rd_01(),M+2-i);
            if(index+i <M) {output.push_back(index+i);}//cout<<"push "<<index+i<<endl;}
            i=index+1+i;
        }
        return output;   // { output= {x} | x in 0,1,...,M-1}
    }

    int insert_rand(double *a,double r,int M){
        int i=0;
        int j=M;
        int mid=0;
        while(i<j-1){
            mid=(i+j)/2;
            if(r< *(a+mid)){j=mid;}
            else{i=mid;}
        }
    return i;
    }
}; //Cumulative probability method
class My_random{
    public:
    my_random_neighbour a;
    My_random(int N0,double K0,double A):a(couples(N0,K0,A),N0-1){
    }
    My_random(int N0,int N1,double K0,double K1,double A):a(couples(N0,N1,K0,K1,A),N0+1){
    }
    My_random(int N0,int N1,int N2,double K0,double K1,double K2,double A):a(couples(N0,N1,N2,K0,K1,K2,A),N0+3){
    }


    std::vector<double> couples(int N0,double K0,double A){
        vector<double> output;
        double pi=atan(1.0)*4.0;//3.14159265
        for(int i=1;i<N0;i++){
            output.push_back(   0.5*A*pow(pi/N0,2)/pow(sin(pi*i/N0),2) );
        }
        output[0]+=K0;
        output[N0-2]+=K0;
        return output;
    }


    vector<double> couples(int N0,int N1,double K0,double K1,double A){
        vector<double> output;
        double pi=atan(1.0)*4.0;//3.14159265
        for(int i=1;i<N0;i++){
            output.push_back(   0.5*A*pow(pi/N0,2)/pow(sin(pi*i/N0),2) );
        }
        output[0]+=K0;
        output[N0-2]+=K0;
        output.push_back(K1);
        output.push_back(K1);
        return output;
    }

    vector<double> couples(int N0,int N1,int N2,double K0,double K1, double K2,double A){
        vector<double> output;
        double pi=atan(1.0)*4.0;//3.14159265
        for(int i=1;i<N0;i++){
            output.push_back(   0.5*A*pow(pi/N0,2)/pow(sin(pi*i/N0),2) );
        }
        output[0]+=K0;
        output[N0-2]+=K0;
        output.push_back(K1);
        output.push_back(K1);
        output.push_back(K2);
        output.push_back(K2);


        return output;
    }

    /*
    My_random(int N0,int N1,double K0,double K1):a(couples(N0,N1,K0,K1),4){
    }
    vector<double> couples(int N0,int N1,double K0,double K1){
        //cout<<"I'm running during initialization"<<endl;
        return {K0,K0,K1,K1};
    }
    */
    vector<int> get_rand(){
        return a.rd_neighbour();
    }
};

//h[N0*N1][N1+1]
void hashing_1d(int **h,int N0){
    int s;
    for(int i=0;i<N0;i++){
            for(int k=0;k<N0-1;k++){
                h[i][k]=(i+k+1)%N0;
            }
    }
}
void hashing_2d(int **h,int N0,int N1){
    int s;
    for(int i=0;i<N0;i++){
        for(int j=0;j<N1;j++){
            s=i+N0*j;
            for(int k=0;k<N0-1;k++){    // N0-1
                h[s][k]=(i+k+1)%N0+N0*j;
            }
            h[s][N0-1]=i+N0*((j-1+N1)%N1);
            h[s][N0]=i+N0*((j+1)%N1);
        }
    }
}
void hashing_3d(int **h,int N0,int N1,int N2){
    int s;
    for(int i=0;i<N0;i++){
        for(int j=0;j<N1;j++){
            for (int n = 0; n < N2; n++) {
                s=i+N0*j+N0*N1*n;

                for(int k=0;k<N0-1;k++){    // N0-1
                    h[s][k]=(i+k+1)%N0+N0*j+N0*N1*n;
                }
                h[s][N0-1]=i+N0*((j-1+N1)%N1)+N0*N1*n;
                h[s][N0]=i+N0*((j+1)%N1)+N0*N1*n;
                h[s][N0+1]=i+N0*j+N0*N1*((n-1+N2)%N2);
                h[s][N0+2]=i+N0*j+N0*N1*((n+1+N2)%N2);
            }
        }
    }


}

double mag_m2(vector<double> m){
    double output=0.0;
    for(auto i:m){
        output+=pow(i,2);
    }
    return output/m.size();
}
double mag_m4(vector<double> m){
    double output=0.0;
    for(auto i:m){
        output+=pow(i,4);
    }
    return output/m.size();
}
void best_and_error(vector<double> v,double &best,double &error){
    best=0.0;
    error=0.0;
    int N=v.size();
    if(N>0){
        for(auto i:v){
            best+=i;
        }
        best=best/N;
        for(auto i:v){
            error+=pow(i-best,2);
        }
        error=pow(error/N,0.5);
    }
}


#include <queue>   // class Ising.updating()  Q
#include <vector> // class Ising.Ns
class Ising{
    public:
    int dimension;
    int Ntotal;
    vector<int> Ns;
    vector<double> Ks;//next neighbour coupling strength
    double A;//longe range coupling strength
    bool *s;//pointer to the Ising field

    int **hashing;//pointer to a 2D array  hashing[i][j]  site i's j-th neighbour
                  //j-th neighbour

    //another possible method
    //bool * **hash;     hash[i][j]=&s[f(i,j)];    *hash[i][j] can be operated directly

    My_random pp;
/*
    Ising(int N0,double K0,double A0):pp(N0,K0,A0){

    }
*/
Ising(int N0,double K0,double A0):pp(N0,K0,A0){
    //neither work for N1=1 nor N1=2,  we should set N1>=3
    //N0 is imaginary time, should be larger
    dimension=1; Ntotal=N0;
    Ns.push_back(N0);
    Ks.push_back(K0);A=A0;
    s=new bool[Ntotal];
    rng.rd_thermo(s,Ntotal);// RANDOM 1
    hashing=new int*[Ntotal];
    for(int i=0;i<Ntotal;i++){
        hashing[i]=new int[(dimension-1)*2+N0-1];
    }
    hashing_1d(hashing,N0);
    //updating(Ntotal*100);
}
    Ising(int N0,int N1,double K0,double K1,double A0):pp(N0,N1,K0,K1,A0){
        //neither work for N1=1 nor N1=2,  we should set N1>=3
        //N0 is imaginary time, should be larger
        dimension=2; Ntotal=N0*N1;
        Ns.push_back(N0);Ns.push_back(N1);
        Ks.push_back(K0);Ks.push_back(K1);A=A0;
        s=new bool[Ntotal];
        rng.rd_thermo(s,Ntotal);// RANDOM 1
        hashing=new int*[Ntotal];
        for(int i=0;i<Ntotal;i++){
            hashing[i]=new int[(dimension-1)*2+N0-1];
        }
        hashing_2d(hashing,N0,N1);
        //updating(Ntotal*100);
    }
    Ising(int N0,int N1,int N2,double K0,double K1,double K2,double A0):pp(N0,N1,N2,K0,K1,K2,A0){
        dimension=3; Ntotal=N0*N1*N2;
        Ns.push_back(N0);Ns.push_back(N1);Ns.push_back(N2);
        Ks.push_back(K0);Ks.push_back(K1);Ks.push_back(K2);A=A0;
        s=new bool[Ntotal];
        rng.rd_thermo(s,Ntotal);// RANDOM 1
        hashing=new int*[Ntotal];
        for(int i=0;i<Ntotal;i++){
            hashing[i]=new int[(dimension-1)*2+N0-1];
        }
        hashing_3d(hashing,N0,N1,N2);
        //updating(Ntotal*100);
    }
    ~Ising(){
        delete [] s;
        for(int i=0;i<Ntotal;i++){
            delete [] hashing[i];
        }
        delete [] hashing;
    }


    void updating(int loops){
      for (size_t i = 0; i < loops; i++) {
        updating();
      }
    }

    void updating(){
          queue<int> Q;
          vector<int> bb;
          int temp=rng.rd_getone(Ntotal);// RANDOM 2

          Q.push(temp);
          s[temp]=!s[temp];
          int site, neighbour;
          while(!Q.empty()){
              site=Q.front();
              Q.pop();
              bb=pp.get_rand(); // RANDOM 3
              for(int i=0;i<bb.size();i++){
                  neighbour=hashing[site][bb[i]];
                  if(s[neighbour]!=s[temp]){         //*hash[site][bb[i]]==color;
                      Q.push(neighbour);
                      s[neighbour]=!s[neighbour];
                  }
              }
          }
    }



    // it still works if <bool> is replaced by <int>

    void printall( ) {
      for (size_t i = 0; i < Ntotal; i++) {
        cout<<s[i];
      }
      cout<<endl;
    }

    double mag_per_site(){
      /*
        std::cout << Ks[0] <<" " << Ks[1] << '\n';
        for (size_t i = 0; i < Ns[0]; i++) {
          for (size_t j = 0; j < Ns[1]; j++) {
            cout<< ((s[i+j*Ns[0]]==true) ? 'X' : '.');
          }
          std::cout << '\n';
        } std::cout << "************" << '\n';
  */
        return accumulate(s,s+Ntotal,0)*2.0/Ntotal-1.0;
    }

    double binder(double threshold){
        //this is a self adaptive calculation
        //if the relative error is larger than the threshold,
        //then increase sec
        int substeps=1000;
        int sec=5;

        vector <vector<double>> magnetization_series(sec);
        double temp_m2;
        double temp_m4;
        vector <double> binder_series(sec);


        double best=1.0;
        double error=1.0;

        for(int i=0;i<substeps;i++){
            updating();
        }
        while( abs(error) >  abs(best*threshold)     ){

            for(int j=0;j<sec;j++){
                for(int i=0;i<substeps;i++){
                    updating(5);
                    magnetization_series[j].push_back(mag_per_site());
                }
                temp_m2=mag_m2(magnetization_series[j]);
                temp_m4=mag_m4(magnetization_series[j]);
                binder_series[j]=( 1.5-0.5*temp_m4/pow(temp_m2,2) );
            }
            best_and_error(binder_series,best,error);
            substeps=substeps*2;
            //std::cout << substeps << '\n';
            //std::cout << best <<  "\t"<< error<< '\n';
        }

        std::cout << best <<  "\t"<< error<< '\n';

        return best;
    }
};



#endif
