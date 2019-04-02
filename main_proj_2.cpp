#include "file_generator.h"
#include "calculate_powerSpectrum.h"
#include "Ising.h"

#include <iostream>
#include <string>
#include <ctime>

int main(int argc, char const *argv[]) {

  //large enough
  int thremoLOOPS=1000*500;
  int LOOPS=1000000;
    
  //set to be optimal by  main_proj_1.cpp
  int sumN=800000; 
  int J=610;
    
  //set seed  from arguments   
  int seed=3;
  seed=std::stoi(argv[1]);  
  rng.mt.seed(seed);

    

  //this is the physics parameters  
  int N0=128;
  int N1=64;
  double K0=0.136;
  double K1=0.2;
  double alpha=0.2;

  Ising e1(N0,N1,  K0, K1 , alpha);
  powerSpectrum pw(N0,N1 );

  /***************************  parameters are above this line  ****************************/
    
  Summation s(e1.Ntotal,sumN);
  writeFiles w(seed);

  e1.updating(thremoLOOPS);

  for (size_t i = 0; i < LOOPS; i++) {
    std::cout << i << '\n';
    s.zero();
    for (size_t j = 0; j < sumN; j++) {
      e1.updating(J);
      pw.calculate(e1.s);
      s.add(pw.m_spectrum);
    }
    s.ave();
    w.save(s.c);
  }

  return 0;
}
