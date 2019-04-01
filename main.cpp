#include "Ising.h"

int main() {


std::vector<int> L0s={50,144}; //{16,36,64};
std::vector<int> L1s={30,12};   //{4,6,8};

std::vector<double> Ks;
for (size_t i = 10; i >0; i--) {
  Ks.push_back(i*0.01+0.40);
}

for(int i=0;i<3;i++){
	for(auto K:Ks){   //10 8 16
	  Ising a(L0s[i],L1s[i],0.136,K,0.2);  // N0*N1*N0  large memory for hashing 0.1361707344557614
	  cout<<K<<'\t';
	  a.binder(0.02);
	}
}



}
