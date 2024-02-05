/*-----------------------------------------------------------------*/
/* DEPENDENCIES */
/*-----------------------------------------------------------------*/
#include "ibex.h"
#include <soplex.h>
#include "cmath"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>  
#include <sciplot/sciplot.hpp>
#include <string>

using namespace std;
using namespace ibex;
using namespace soplex;
using namespace sciplot;

/*-----------------------------------------------------------------*/
/* CUSTOM OBJECTS */
/*-----------------------------------------------------------------*/
//Description of the multivariate polynomial
struct Multivariate_Polynomiale {
  std::vector<double> value; //Value of a
  double (*P)(std::vector<double> x, struct Multivariate_Polynomiale A); //Polynomial function
} MultivariatepolynomialDescription;

//Description of the function to approximate
struct Multivariate_functionDescription {
  string functionString; // With the unknown x1 and x2
  std::vector<string> derivedFunctionStrings;
  double (*f)(std::vector<double> x); // Must correspond to the string version
} MultivariatefunctionDescription;

//Parameters used by both Remez I and Remez I+
struct Multivariate_RemezIParameters {
  std::vector<double>(*phi)(std::vector<double> x, int n);
  struct Multivariate_functionDescription fdesc; //description of the function to approximate
  int approximationDegree; //a_i with i in [0 ; approximationDegree], as many i as nbVar
  int sizeD0; //Advised to be >= to approximationDegree+2 for quickest approximation
  int maxNbTurns; //Failsafe
  int nbVar;
  std::vector<std::vector<double>> bornersVar;
  double approximationResult;
  double approximationPoints; 
} Multivariate_RemezIParameters;


//Results of RemezI multivariate for comparisons and graphs
struct Multivariate_RemezIResult {
  int nbTurns;
  std::vector<double> errorStep1;
  std::vector<double> errorStep2;
  bool remezPlus;
  bool exchange;
  string path;
} Multivariate_RemezIResult;

/*-----------------------------------------------------------------*/
/* FUNCTIONS TO APPROXIMATE */
/*-----------------------------------------------------------------*/

std::vector<double> phi(std::vector<double> x, int n);
void write_introduction(string path, struct Multivariate_RemezIParameters remezdesc);
string createPath(struct Multivariate_RemezIParameters remezdesc);
std::vector<std::vector<double>> step0Multi(struct Multivariate_RemezIParameters remezdesc);
void writeD0(struct Multivariate_RemezIParameters remezdesc, std::vector<std::vector<double>> D, string path);
std::vector<double> floorPhi(std::vector<double> x, int degree);
std::vector<double>  step1(std::vector<std::vector<double>> D, int degree, string path);


double f(std::vector<double> x) {
  return pow(x[0],2) + pow(x[1], 2);
}

double P(struct Multivariate_Polynomiale A , std::vector<double> x){
  double y = 0;
  double temp = 0;
  std::vector<double> p = phi(x, A.value.size());
  for(int i = 0 ; i < A.value.size(); i++){
    temp = A.value[i];
    temp = temp * p[i];
    y += temp;
  }
  return y;
}

double error(struct Multivariate_Polynomiale A , std::vector<double> x){
  return f(x)-P(A, x);
}

std::vector<double> phi(std::vector<double> x, int degree){
  if(degree>0){
    std::vector<double> p0 = phi(x, degree-1);
    std::vector<double> p1 = floorPhi(x, degree);
    p0.reserve(p0.size()+p1.size());
    p0.insert(p0.end(), p1.begin(), p1.end());
    return p0;
  } else {
    return {1};
  }
}

std::vector<double> floorPhi(std::vector<double> x, int degree){
  if(degree == 0){
  return {1};
  } else {
      if(degree == 1){
        return x;
      } else {
        std::vector<double> prev = floorPhi(x, degree-1);
        std::vector<double> actual = {};
        for(int i = 0; i<x.size(); i++){
          for(int j=0; j<prev.size(); j++){
            actual.push_back(x[i]*prev[j]);
          }
        }
        return actual;
      }
  }
}

/*-----------------------------------------------------------------*/
/* PARAMETERS AND INPUTS */
/*-----------------------------------------------------------------*/

struct Multivariate_functionDescription initialize_fdesc_multi(){
  struct Multivariate_functionDescription fdesc;
  fdesc.functionString = "x1^2 + x2^2";
  fdesc.derivedFunctionStrings = {"2*x1 + 2*x2"};
  fdesc.f = (*f);
  return fdesc;
}

struct Multivariate_RemezIParameters initialize_remezdesc_multi(struct Multivariate_functionDescription fdesc){
  struct Multivariate_RemezIParameters remezdesc;
  remezdesc.approximationDegree = 1;
  remezdesc.phi = (*phi);
  remezdesc.sizeD0 = 3 ; //Warning, if exchange, sizeD0 value will be replaced by degree+2
  remezdesc.maxNbTurns = 100;
  remezdesc.nbVar = 2;
  remezdesc.approximationResult = 1.e-10;
  remezdesc.approximationPoints = 1.e-15;
  remezdesc.bornersVar = {{-1 , 1}, {-1, 1}}; //Warning, remember to put in crescent order
  remezdesc.fdesc = fdesc;
  return remezdesc;
}
/*-----------------------------------------------------------------*/
/* MAIN AND COMMON FUNCTIONS */
/*-----------------------------------------------------------------*/

int main(int argc, char** argv) { 
  struct Multivariate_functionDescription fdesc;
  fdesc = initialize_fdesc_multi();
  struct Multivariate_RemezIParameters remezdesc;
  remezdesc = initialize_remezdesc_multi(fdesc);
  string path=createPath(remezdesc);
  write_introduction(path, remezdesc);
  std::vector<std::vector<double>> D = step0Multi(remezdesc);
  writeD0(remezdesc, D, path);
  cout << phi(D[0], remezdesc.approximationDegree).size() << endl;
  step1(D, remezdesc.approximationDegree, path);
}

std::vector<std::vector<double>> step0Multi(struct Multivariate_RemezIParameters remezdesc){
  std::vector<double> steps = {};
  for (int i = 0; i<remezdesc.nbVar; i++){
    steps.push_back((remezdesc.bornersVar[i][1] - remezdesc.bornersVar[i][0])/(remezdesc.sizeD0-1));
  }
  std::vector<std::vector<double>> D0 = {};
  std::vector<double> d = {};
  for(int i = 0; i < remezdesc.sizeD0; i ++){
    for(int j = 0; j < remezdesc.nbVar; j++){
      d.push_back(remezdesc.bornersVar[j][0] + steps[j]*i);
    }
    D0.push_back(d);
    d = {};
  }
  return D0;
}

std::vector<double> step1(std::vector<std::vector<double>> D, int degree, string path){
   std::vector<std::vector<double>> p = {};
   for(int j = 0; j< D.size(); j++){
     p.push_back(phi(D[j], degree));
   }
   SoPlex mysoplex;
   /* parameters */
   mysoplex.setIntParam(SoPlex::READMODE, SoPlex::READMODE_REAL);
   mysoplex.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_REAL);
   mysoplex.setIntParam(SoPlex::CHECKMODE, SoPlex::CHECKMODE_REAL);
   mysoplex.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
   mysoplex.setRealParam(SoPlex::FEASTOL, 0.0);
   mysoplex.setRealParam(SoPlex::OPTTOL, 0.0);
   /* Objective sense */
   mysoplex.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
   /* Variables */
   DSVectorReal dummycol(0);  
   for(int i = 0; i<=p[0].size(); i++){
     mysoplex.addColReal(LPColReal(0, dummycol, infinity, -infinity));
   }
   mysoplex.addColReal(LPColReal(1, dummycol, infinity, -infinity));
   /* Constraints */
   double x;
   Real r;
   for(int i = 0 ; i<D.size(); i++){ 
     DSVectorReal row1(0);
     for(int j = 0 ; j<=p[i].size(); j++){
       x = p[i][j];
       r = x;
       row1.add(j, r);
     }
     r = 1;
     row1.add(p[i].size()+1, r);
     r = f(D[i]);
     mysoplex.addRowReal(LPRowReal(r, row1, infinity));
   }
   for(int i = 0 ; i<D.size(); i++){
     DSVectorReal row1(0);
     for(int j = 0 ; j<=p[i].size(); j++){
       x = p[i][j];
       r = -x;
       row1.add(j, r);
     }
     r = 1;
     row1.add(p[i].size()+1, r);
     r = -f(D[i]);
     mysoplex.addRowReal(LPRowReal(r, row1, infinity));
   }
   string fullPathToModel = path + "models/step1Soplex.lp";
   mysoplex.writeFileReal(fullPathToModel.c_str(), NULL, NULL, NULL);
   SPxSolver::Status stat;
   DVectorRational prim(degree+2);
   DVectorRational dual(D.size()*2);
   stat = mysoplex.optimize();
   std::vector<double> primDouble;
   /*Solution */
   if(stat == SPxSolver::OPTIMAL)
   {
      mysoplex.getPrimalRational(prim);
      mysoplex.getDualRational(dual);
      std::cout << "LP solved to optimality.\n";
      std::cout << "Objective value is " << mysoplex.objValueRational() << ".\n";
      cout << "Primal solution is [ "; 
      for(int j = 0 ;  j<=degree; j++){
        cout << prim[j] << ", ";
        primDouble.push_back(prim[j]);
      }
      cout << " ]" << endl;
      std::cout << "Error solution is [" << prim[degree+1] << "].\n";
      primDouble.push_back(prim[degree+1]);
      return primDouble;
   }
   else
   {
      std::cout << "Error: SoPlex returned with status " << stat << ".\n";
      return primDouble;
   }   
}


/*-----------------------------------------------------------------*/
/* FUNCTIONS GRAPHS CREATIONS*/
/*-----------------------------------------------------------------*/


/*-----------------------------------------------------------------*/
/* FUNCTIONS FILE MANAGEMENT*/
/*-----------------------------------------------------------------*/
void write_introduction(string path, struct Multivariate_RemezIParameters remezdesc){
  std::ofstream summary;
  summary.open(path + "summary/" + remezdesc.fdesc.functionString + ".txt", std::ios_base::app);
  summary << "Welcome to the summary of ";
  summary << "Remez 1 multi ";
  summary << "for multivariate polynomial approximation. \n";
  summary << "With the following parameters : \n";
  summary << "degree of approximation = " << remezdesc.approximationDegree << "\n";
  summary << "size of first discretization =" << remezdesc.sizeD0 << "\n";
  summary << "maximum number of turns = " << remezdesc.maxNbTurns << "\n";
  summary << "approximation of error = " << remezdesc.approximationResult << "\n";
  summary << "approximation of points = " << remezdesc.approximationPoints << "\n";
  summary << "x in = [";
  for(int i = 0; i<remezdesc.bornersVar.size(); i++){
      summary << "[" << remezdesc.bornersVar[i][0] << "," << remezdesc.bornersVar[i][1] << "]";
  }
  summary << "]\n";
  summary << "For the function " << remezdesc.fdesc.functionString << ".\n";
}

string createPath(struct Multivariate_RemezIParameters remezdesc){
  string path = "../resultsOfSteps/multi/";
  return path;
}

void writeD0(struct Multivariate_RemezIParameters remezdesc, std::vector<std::vector<double>> D, string path){
  std::ofstream summary;
  summary.open(path + "summary/" + remezdesc.fdesc.functionString + ".txt", std::ios_base::app);
  summary << "Generated D0 is : ";
  for(int i=0; i<D.size(); i++){
    summary << "[";
    for(int j=0; j<D[i].size(); j++){
      summary << D[i][j] << ", ";
    }
    summary << "] ";
  }
}

























