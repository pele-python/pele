
#include "_lbfgs_cpp.h"
#include <vector>
#include <iostream>

using std::vector;
using std::cout;

int main(){
  int N = 3;
  int M = 4;
  vector<double> x0(N);

  x0[0] = .13;
  x0[1] = 1.23;
  x0[2] = -5.03;


  cout << "starting\n";
  LBFGS_ns::LBFGS lbfgs(&x0[0], N, M);
  lbfgs.run();
  //for (int i=0; i<30; ++i){
    //lbfgs.one_iteration();
  //}

}
