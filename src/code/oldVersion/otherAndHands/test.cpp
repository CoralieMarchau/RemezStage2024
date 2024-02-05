#include "ibex.h"
#include <glpk.h>
#include "cmath"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>  
#include<string>

using namespace std;
using namespace ibex;
  
int main(int argc, char** argv) { 
  System sys("step2_model.txt");
  DefaultOptimizer optimizer(sys,1e-01);
  optimizer.optimize(sys.box);
  optimizer.report();
}

