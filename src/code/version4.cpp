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
#include "gmp.h"

using namespace std;
using namespace ibex;
using namespace soplex;
using namespace sciplot;
/*-----------------------------------------------------------------*/
/* CUSTOM OBJECTS */
/*-----------------------------------------------------------------*/

//Description of the multivariate polynomial
struct Multivariate_Polynomiale {
  std::vector<double> a; //Value of a
  std::vector<double>(*phi)(std::vector<double> x, int degree);
  std::vector<std::vector<int>> (*derivatedPhi)(int nbX, int degree);
  int degree;
} MultivariatepolynomialDescription;

//Description of the function to approximate
struct Multivariate_functionDescription {
  string functionString; // With the unknown x1 and x2
  int nbX;
  std::vector<string> derivedFunctionStrings;
  std::vector<std::vector<double>> bornersVar;
  std::vector<bool> notConstantDerived;
  double (*f)(std::vector<double> x); // Must correspond to the string version
} MultivariatefunctionDescription;

struct Multivariate_RemezIIParameters {
  struct Multivariate_Polynomiale poly;
  struct Multivariate_functionDescription fdesc; //description of the function to approximate
  int approximationDegree; //a_i with i in [0 ; approximationDegree], as many i as nbVar
  int maxNbTurns; //Failsafe
  double approximationResult;
  double approximationPoints; 
} Multivariate_RemezIIParameters;

//Results of RemezII multivariate for comparisons and graphs
struct Multivariate_RemezIIResult {
  int nbTurns;
  int currentTurn;
  std::vector<double> errorStep1;
  std::vector<double> errorStep2;
  string path;
} Multivariate_RemezIIResult;

/*-----------------------------------------------------------------*/
/* FUNCTIONS DECLARATION */
/*-----------------------------------------------------------------*/
struct Multivariate_functionDescription initialize_fdesc(string fString, std::vector<string> df, double (*f)(std::vector<double> x), std::vector<std::vector<double>> bornersVar, std::vector<bool> nCD);
struct Multivariate_RemezIIParameters initialize_remezdesc(struct Multivariate_functionDescription fdesc, int degree, int nbTurns);
struct Multivariate_Polynomiale initialize_poly(std::vector<double>(*phi)(std::vector<double> x, int degree), std::vector<std::vector<int>> (*derivatedPhi)(int nbX, int degree));
std::vector<std::vector<double>> step0(int nbX, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_functionDescription fdesc);
string createPath();
std::vector<double> step1(std::vector<std::vector<double>> D, struct Multivariate_functionDescription fdesc, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_Polynomiale poly);
void write_introduction(string path, struct Multivariate_RemezIIParameters remezdesc);
void writeD(struct Multivariate_RemezIIParameters remezdesc, std::vector<std::vector<double>> D, string path);
std::vector<std::vector<double>> step2(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path, int turn);
void createErrorGraph(struct Multivariate_RemezIIParameters remezdesc, string path);
void writeA(std::vector<double> a, string path, string functionString);
std::vector<std::vector<double>> step0_random(int nbPoints, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_functionDescription fdesc);
void write_data_for_graphs(std::vector<double> A, std::vector<std::vector<double>> bornersVar, int turn);
std::vector<std::vector<double>> step2_center(struct Multivariate_RemezIIParameters remezdesc, string path, int turn);
void createModelStep2(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path);
void createModelStep2_side(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc,double borner, string path, int var);
std::vector<std::vector<double>> step2_side(struct Multivariate_RemezIIParameters remezdesc, string path, int turn);
std::vector<std::vector<double>> noCopyAdd(std::vector<std::vector<double>> v1, std::vector<std::vector<double>> v2);
double step2Error(std::vector<std::vector<double>> newPoints, std::vector<double> a, struct Multivariate_RemezIIParameters remezdesc, string path);
void writePolynom(std::vector<double> A, string path, struct Multivariate_RemezIIParameters remezdesc);
void writeErrors(double error1, double error2, struct Multivariate_RemezIIParameters remezdesc, string path);
std::vector<std::vector<double>> step2_centerDIRECT(struct Multivariate_RemezIIParameters remezdesc, string path, int turn);
void createModelStep2Error(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path);
void grapheConvergenceComparisons(struct Multivariate_RemezIIResult remezResults);
bool cheking_side_not_constant(struct Multivariate_RemezIIParameters remezdesc, int side, int var);
/*-----------------------------------------------------------------*/
/* FUNCTIONS TO APPROXIMATE */
/*-----------------------------------------------------------------*/
//Didactical Example

double f0(std::vector<double> x) {
  return pow(x[0],2) + pow(x[1],2);
}

struct Multivariate_functionDescription initializeF0(){
  return initialize_fdesc("x1^2 + x2^2",{"2*x1","2*x2", "2", "2", "0"},f0, {{0,1},{0,1}}, {true, true, true, true});
}

//Reemtsen Examples
//Constant for x1 = 0 only
double f1(std::vector<double> x) {
  return log(x[0]+x[1])*sin(x[0]);
}

struct Multivariate_functionDescription initializeF1(){
  return initialize_fdesc("ln(x1+x2)*sin(x1)",{"((sin(x1)/(x1+x2))+(ln((x1+x2))*cos(x1)))","(sin(x1)/(x1+x2))"},f1, {{0,1},{1, 2.5}}, {false, true, true, true});
}

//Constant for x1 = 0
double f2(std::vector<double> x) {
  return pow((1+x[0]), x[1]);
}
struct Multivariate_functionDescription initializeF2(){
  return initialize_fdesc("(1+x1)^(x2)",{"x2*(x1+1)^(x2-1)","ln(1+x1)*(1+x1)^(x2)"},f2, {{0,1},{1, 2.5}}, {false, true, true, true});
}

//Constant for x3 = 1
double f3(std::vector<double> x) {
  return cos(x[2])*pow(1+x[0],x[1]);
}
struct Multivariate_functionDescription initializeF3(){
  return initialize_fdesc("cos(x3)*(1+x1)^(x2)",{"x2*cos(x3)*(1+x1)^(x2-1)"," cos(x3)*(1+x1)^(x2)*ln(1+x1)","(-sin(x3)*(1+x1)^(x2))"},f3, {{0,1},{1, 2},{0,1}}, {true, true, true, true, false, true});
}

double f4(std::vector<double> x) {
  return 1/(x[0]+2*x[1]+4);
}
struct Multivariate_functionDescription initializeF4(){
  return initialize_fdesc("1/(x1+2*x2+4)",{"(-(1/((4+x1)+(2*x2))^2))","(-(2/((4+x1)+(2*x2))^2))"},f4, {{-1,1},{-1,1}}, {true, true, true, true});
}

//Constant for x1 = 0 (not a side so not important)
double f5(std::vector<double> x) {
  return exp(pow(x[0],2)+x[0]*x[1]);
}
struct Multivariate_functionDescription initializeF5(){
  return initialize_fdesc("exp(x1^2 + x1*x2)",{"((2*(x1*exp((x1^2+(x1*x2)))))+(x2*exp((x1^2+(x1*x2)))))", "(x1*exp((x1^2+(x1*x2))))"}, f5, {{-1,1},{-1,1}}, {false, true, true, true});
} 

double f6_7(std::vector<double> x) {
  return sqrt(x[0]+2*x[1]+4);
}
struct Multivariate_functionDescription initializeF6_7(){
  return initialize_fdesc("sqrt(x1+2*x2+4)",{"1/(2*sqrt(x1+2*x2+4))","1/(sqrt(x1+2*x2+4))"},f6_7, {{-1,1},{-1,1}}, {true, true, true, true});
}


//Constant for x2 = 0
double f8(std::vector<double> x) {
  return abs(log((x[0]*x[1]+1)/(x[0]+0.5)))*pow(x[1],(x[2]+1)/2);
}
struct Multivariate_functionDescription initializeF8(){
  return initialize_fdesc("",{},f8, {{0,1},{0,1},{0,1}}, {true, true, false, true, true, true});
}

/*-----------------------------------------------------------------*/
/* BASICS */
/*-----------------------------------------------------------------*/

double P(struct Multivariate_Polynomiale A , std::vector<double> x){
  double y = 0;
  std::vector<double> p = A.phi(x, A.a.size());
  for(int i = 0 ; i < A.a.size(); i++){
    y += A.a[i] * p[i];;
  }
  return y;
}

double error(struct Multivariate_Polynomiale A , std::vector<double> x, double (*f)(std::vector<double> x)){
  return abs(f(x)-P(A, x));
}

/*-----------------------------------------------------------------*/
/* POLYNOMIALS */
/*-----------------------------------------------------------------*/

std::vector<double> phiBasique(std::vector<double> x, int degree){
  std::vector<double> p = {};
  for(int i = 0; i<=degree; i=i+1){
    for(int j = 0; j<=degree-i; j=j+1){
      p.push_back(pow(x[0],i)*pow(x[1],j));
    }
  }
  return p;
}

std::vector<std::vector<int>> phiBasiqueDerivated(int nbX, int degree){
  std::vector<std::vector<int>> powOfX = {};
  std::vector<std::vector<int>> temp = {};
  std::vector<std::vector<int>>  prev = {std::vector<int>(nbX, 0.0)};
  for(int i = 0; i<=degree; i=i+1){
    for(int j = 0; j<=degree-i; j=j+1){
      powOfX.push_back({i,j});
    }
  }
  
  return powOfX;
  
}

struct Multivariate_Polynomiale initializePhiBasic(){
  return initialize_poly(phiBasique, phiBasiqueDerivated);
}

/*-----------------------------------------------------------------*/
/* INITIALIZES OBJECTS */
/*-----------------------------------------------------------------*/
struct Multivariate_functionDescription initialize_fdesc(string fString, std::vector<string> df, double (*f)(std::vector<double> x), std::vector<std::vector<double>> bornersVar, std::vector<bool> nCD){
  struct Multivariate_functionDescription fdesc;
  fdesc.functionString = fString;
  fdesc.derivedFunctionStrings = df;
  fdesc.f = f;
  fdesc.notConstantDerived = nCD;
  fdesc.bornersVar = bornersVar;
  fdesc.nbX = bornersVar.size();
  return fdesc;
}

struct Multivariate_RemezIIParameters initialize_remezdesc(struct Multivariate_functionDescription fdesc, int degree, int nbTurns){
  struct Multivariate_RemezIIParameters remezdesc;
  remezdesc.approximationDegree = degree;
  remezdesc.maxNbTurns = nbTurns;
  remezdesc.approximationResult = 1.e-10;
  remezdesc.approximationPoints = 1.e-15;
  remezdesc.fdesc = fdesc;
  return remezdesc;
}

struct Multivariate_Polynomiale initialize_poly(std::vector<double>(*phi)(std::vector<double> x, int degree), std::vector<std::vector<int>> (*derivatedPhi)(int nbX, int degree)){
  struct Multivariate_Polynomiale poly;
  poly.a = {};
  poly.phi = (*phi);
  poly.derivatedPhi = (*derivatedPhi);
  return poly;
}

struct Multivariate_RemezIIResult initialize_results(int nbTurnsMax, string path){
  struct Multivariate_RemezIIResult result;
  result.nbTurns = nbTurnsMax;
  result.currentTurn = 0;
  result.errorStep1 = {};
  result.errorStep2 = {};
  result.path = path;
  return result;
}
/*-----------------------------------------------------------------*/
/* MAIN */
/*-----------------------------------------------------------------*/

int main(int argc, char** argv) { 
  int nbTurns = 15;
  int degree = 2;
  double errorStep2 = 0;
  double errorStep1 = 0; 
  string path = createPath();
  std::vector<std::vector<double>> newPoints;
  struct Multivariate_RemezIIParameters remezdesc;
  struct Multivariate_RemezIIResult results;
  results = initialize_results(nbTurns, path);
  int typeOfD0 = 1;
  remezdesc = initialize_remezdesc(initializeF6_7(), degree, nbTurns);
  remezdesc.poly = initializePhiBasic();
  remezdesc.poly.degree = degree;
  write_introduction(path, remezdesc);
  
  //for(int i = typeOfD0; i<=4; i++){
    //int i = 5;
    std::vector<std::vector<double>> D ;
    int nbPoints = 5;
      
      D = noCopyAdd(step0_random(5-nbPoints, remezdesc, remezdesc.fdesc), step0(remezdesc.fdesc.nbX, remezdesc, remezdesc.fdesc));
    		std::vector<double> phi;
  for(int j = 0; j<nbTurns; j++){
    writeD(remezdesc, D, path);
    remezdesc.poly.a = step1(D, remezdesc.fdesc, remezdesc, remezdesc.poly);
      
  //Taking the error out of the polynomial
    errorStep1 = remezdesc.poly.a[remezdesc.poly.a.size()-1];
    remezdesc.poly.a.pop_back();
  
    writePolynom(remezdesc.poly.a, path, remezdesc);
  //createErrorGraph(remezdesc, path);
    write_data_for_graphs(remezdesc.poly.a, remezdesc.fdesc.bornersVar, j+1);
     
    newPoints = step2(remezdesc.poly.a, remezdesc, path, j+1);
    errorStep2 = step2Error(newPoints, remezdesc.poly.a, remezdesc, path);
    D = noCopyAdd(newPoints, D);
    writeErrors(errorStep1, errorStep2, remezdesc, path);
    results.errorStep1.push_back(errorStep1);
    results.errorStep2.push_back(errorStep2);
   //}
    cout << "\n---------------------------------------------------------------\n";
  }
  grapheConvergenceComparisons(results);
}

/*-----------------------------------------------------------------*/
/* REMEZ II */
/*-----------------------------------------------------------------*/

std::vector<std::vector<double>> step0(int nbX, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_functionDescription fdesc){
  std::vector<std::vector<double>> D0 = {};
  std::vector<double> d = {};
  for(int i = 0; i < nbX; i ++){ //middle
    d.push_back(((fdesc.bornersVar[i][1] - fdesc.bornersVar[i][0])/2)+fdesc.bornersVar[i][0]);
  }
  D0.push_back(d);
  d = {};
  for(int i = 0; i < nbX; i ++){ //all min
    d.push_back(fdesc.bornersVar[i][0]);
  }
  D0.push_back(d);
  d = {};
  
  for(int i = 0; i < nbX; i ++){ //all max
    d.push_back(fdesc.bornersVar[i][1]);
  }
  D0.push_back(d);
  d = {};
  
  for(int j = 0; j < nbX; j ++){ 
    //j only max
    for(int i = 0; i < j; i ++){ 
      d.push_back(fdesc.bornersVar[i][0]);
    }
    d.push_back(fdesc.bornersVar[j][1]);
    for(int i = j+1; i < nbX; i ++){ 
      d.push_back(fdesc.bornersVar[i][0]);
    }
    D0.push_back(d);
    d = {};
    //all until j at max all after at min
    for(int i = 0; i < j; i ++){ 
      d.push_back(fdesc.bornersVar[i][1]);
    }
    d.push_back(fdesc.bornersVar[j][1]);
    for(int i = j+1; i < nbX; i ++){ 
      d.push_back(fdesc.bornersVar[i][0]);
    }
    D0.push_back(d);
    d = {};
    //all except j at max
    for(int i = 0; i < j; i ++){ 
      d.push_back(fdesc.bornersVar[i][1]);
    }
    d.push_back(fdesc.bornersVar[j][0]);
    for(int i = j+1; i < nbX; i ++){ 
      d.push_back(fdesc.bornersVar[i][1]);
    }
    D0.push_back(d);
    d = {};
  }
  
  return D0;
}

std::vector<std::vector<double>> step0_random(int nbPoints, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_functionDescription fdesc){
  srand (42);
  std::vector<std::vector<double>> D0 = {};
  std::vector<double> temp = {};
  for(int i = 0; i<nbPoints; i++){
    for(int j = 0; j < fdesc.bornersVar.size(); j ++){ 
      temp.push_back(((double)rand()/RAND_MAX)*(fdesc.bornersVar[j][1] - fdesc.bornersVar[j][0])+fdesc.bornersVar[j][0]);
    }
    D0.push_back(temp);
    temp = {};
  }
  return D0;
}

std::vector<double> step1(std::vector<std::vector<double>> D, struct Multivariate_functionDescription fdesc, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_Polynomiale poly){
  string path = createPath();
  std::vector<std::vector<double>> p = {};
  for(int j = 0; j< D.size(); j++){
    p.push_back(poly.phi(D[j], poly.degree));
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
  for(int i = 0; i<p[0].size(); i++){
     mysoplex.addColReal(LPColReal(0, dummycol, infinity, -infinity));
   }
   mysoplex.addColReal(LPColReal(1, dummycol, infinity, -infinity));
   /* Constraints */
   Real r;
   for(int i = 0 ; i<D.size(); i++){ 
     DSVectorReal row1(0);
     for(int j = 0 ; j<p[i].size(); j++){
       r = p[i][j];
       row1.add(j, r);
     }
     r = 1;
     row1.add(p[i].size(), r);
     r = fdesc.f(D[i]);
     mysoplex.addRowReal(LPRowReal(r, row1, infinity));
   }
   for(int i = 0 ; i<D.size(); i++){
     DSVectorReal row1(0);
     for(int j = 0 ;  j<p[i].size(); j++){
       r = -p[i][j];
       row1.add(j, r);
     }
     r = 1;
     row1.add(p[i].size(), r);
     r = -fdesc.f(D[i]);
     mysoplex.addRowReal(LPRowReal(r, row1, infinity));
   }
   string fullPathToModel = path + "models/step1Soplex.lp";
   mysoplex.writeFileReal(fullPathToModel.c_str(), NULL, NULL, NULL);
   SPxSolver::Status stat;
   DVectorRational prim(p[0].size()+1);
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
      for(int j = 0 ;  j<p[0].size(); j++){
        cout << prim[j] << ", ";
        primDouble.push_back(prim[j]);
      }
      cout << " ]" << endl;
      std::cout << "Error solution is [" << prim[primDouble.size()] << "].\n";
      primDouble.push_back(prim[primDouble.size()]);
      return primDouble;
   }
   else
   {
      std::cout << "Error: SoPlex returned with status " << stat << ".\n";
      return primDouble;
   }   
}

std::vector<std::vector<double>> step2(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path, int turn){
  std::vector<std::vector<double>> pointsSide = {};
  std::vector<std::vector<double>> pointsCenter = {};
  //Getting the points
  pointsCenter = step2_center(remezdesc, path, turn);
  pointsSide = step2_side(remezdesc, path, turn);
  
  //Getting every points together
  pointsCenter.insert(pointsCenter.end(), pointsSide.begin(), pointsSide.end());
 
  return pointsCenter; 
}

std::vector<std::vector<double>> step2_side(struct Multivariate_RemezIIParameters remezdesc, string path, int turn){
  std::ofstream data;
  data.open("data/"+to_string(turn)+"points.txt",  std::ios_base::app);
  std::vector<std::vector<double>> newPoints = {};
  std::vector<double> point;
  string fullpath;
  for(int i = 0; i  < remezdesc.fdesc.bornersVar.size(); i++){
    if(cheking_side_not_constant(remezdesc, i*2 + 0, i)){
      createModelStep2_side(remezdesc.poly.a, remezdesc, remezdesc.fdesc.bornersVar[i][0], path, i);
    
      fullpath = path + "models/step2_model_side.txt";
      System system(fullpath.c_str());
      DefaultSolver solver(system,1e-09, 1e-15);
      solver.solve(system.box);
      solver.report();
      for(int i = 0; i<solver.get_data().size(); i++){
        for(int j = 0; j<solver.get_data()[i].size(); j++){
          point.push_back(solver.get_data()[i][j].mid());
          data << solver.get_data()[i][j].mid() << " ";
        }
        data << "\n";
        newPoints.push_back(point);
        point = {};
      }
    }
    if(cheking_side_not_constant(remezdesc, i*2 + 1, i)){
      createModelStep2_side(remezdesc.poly.a, remezdesc, remezdesc.fdesc.bornersVar[i][1], path, i);
    
      fullpath = path + "models/step2_model_side.txt";
      System system2(fullpath.c_str());
      DefaultSolver solver2(system2,1e-09, 1e-15);
      solver2.solve(system2.box);
      solver2.report();
    
      for(int i = 0; i<solver2.get_data().size(); i++){
        for(int j = 0; j<solver2.get_data()[i].size(); j++){
          point.push_back(solver2.get_data()[i][j].mid());
          data << solver2.get_data()[i][j].mid() << " ";
        }
        data << "\n";
        newPoints.push_back(point);
        point = {};
      }
    }
  }
  return newPoints;
}

void createModelStep2_side(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, double borner, string path, int var){
  ofstream step2(path + "models/step2_model_side.txt");
  step2 << setprecision(13);
  
  step2 << "Variables\n";
  for(int i = 0; i<remezdesc.fdesc.bornersVar.size(); i++){
    if(i != var){
      step2 << "  x"<< i+1 <<" in [" << remezdesc.fdesc.bornersVar[i][0] <<" , " << remezdesc.fdesc.bornersVar[i][1] <<"];\n";
    } else {
      step2 << "  x"<< i+1 <<" in [" << borner <<" , " << borner <<"];\n";
    }
  }
  step2 << "\nConstraints\n";
  std::vector<std::vector<int>> p = remezdesc.poly.derivatedPhi(remezdesc.fdesc.bornersVar.size(), remezdesc.approximationDegree);
  
  
  int i = 1 - var;
  step2 << remezdesc.fdesc.derivedFunctionStrings[i] << "-("; //derivated of f depending xi
  for(int j=0; j<p.size(); j++){//P'(x) depending of each x
    if(p[j][i]>0){
      step2 << " + " << A[j] ;
      for(int k=0; k<p[j].size(); k++){
        if(k==i){
          if(p[j][k]>1){
            step2 << "*((x" << k+1 << "*" << p[j][k] << ")^" << p[j][k]-1 << ")";
           }
         } else {
           if(k==var){
             if(p[j][k]!=0){
               step2 << "*((" << borner << ")^" << p[j][k] << ")";
              }
            }else{
              if(p[j][k]>=1){
                step2 << "*((x" << k+1 << ")^" << p[j][k] << ")";
              }
            }
          }
      }
    }
  }
    step2 << ")=0;\n"; 
  step2 << "end";
}

std::vector<std::vector<double>> step2_center(struct Multivariate_RemezIIParameters remezdesc, string path, int turn){
  std::ofstream data;
  data.open("data/"+to_string(turn)+"points.txt",std::ofstream::trunc);
  std::vector<std::vector<double>> newPoints = {};
  std::vector<double> point= {};
  createModelStep2(remezdesc.poly.a, remezdesc, path);
  string fullpath = path + "models/step2_model.txt";
  System system(fullpath.c_str());
  DefaultSolver solver(system,1e-09, 1e-15);
  solver.solve(system.box);
  solver.report();
  for(int i = 0; i<solver.get_data().size(); i++){
    for(int j = 0; j<solver.get_data()[i].size(); j++){
      point.push_back(solver.get_data()[i][j].mid());
      data << solver.get_data()[i][j].mid() << " ";
    }
    data << "\n";
    newPoints.push_back(point);
    point = {};
  }
  return newPoints;
}
/*
std::vector<std::vector<double>> step2_centerDIRECT(struct Multivariate_RemezIIParameters remezdesc, string path, int turn){
  std::ofstream data;
  data.open("data/"+to_string(turn)+"points.txt",std::ofstream::trunc);
  std::vector<std::vector<double>> newPoints = {};
  std::vector<double> point= {};
  Function f("x1","x2", remezdesc.fdesc.functionString);
  Function df(f,Function::DIFF);

  Variable x1, x2;
  SystemFactory fac;
  
  fac.add_var(x1);
  fac.add_var(x2);
  
  NumConstraint c(df,ibex::EQ); // the constraint x+1<=0
  System system(fac);
  DefaultSolver solver(system,1e-09, 1e-15);
  solver.solve(system.box);
  solver.report();
  for(int i = 0; i<solver.get_data().size(); i++){
    for(int j = 0; j<solver.get_data()[i].size(); j++){
      point.push_back(solver.get_data()[i][j].mid());
      data << solver.get_data()[i][j].mid() << " ";
    }
    data << "\n";
    newPoints.push_back(point);
    point = {};
  }
  return newPoints;
}*/

void createModelStep2(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path){
  ofstream step2(path + "models/step2_model.txt");
  step2 << setprecision(13);
  step2 << "Variables\n";
  for(int i = 0; i<remezdesc.fdesc.bornersVar.size(); i++){
    step2 << "  x"<< i+1 <<" in [" << remezdesc.fdesc.bornersVar[i][0] <<" , " << remezdesc.fdesc.bornersVar[i][1] <<"];\n";
  }
  step2 << "\nConstraints\n";
  std::vector<std::vector<int>> p = remezdesc.poly.derivatedPhi(remezdesc.fdesc.bornersVar.size(), remezdesc.approximationDegree);
  //1st derivative test
  for(int i = 0; i<remezdesc.fdesc.nbX; i=i+1){
    step2 << remezdesc.fdesc.derivedFunctionStrings[i] << "-("; //derivated of f depending xi
    for(int j=0; j<p.size(); j++){//P'(x) depending of each x
      if(p[j][i]>0){
        step2 << " + " << A[j] ;
        for(int k=0; k<p[j].size(); k++){
          if(k==i){
            if(p[j][k]>1){
              step2 << "*(x" << k+1 << "*" << p[j][k] << ")^" << p[j][k]-1;
            }
          } else {
            if(p[j][k]>=1){
              step2 << "*(x" << k+1 << ")^" << p[j][k];
            }
          }
        }
      }
    }
    step2 << ")=0;\n";
  }
  
  step2 << "end";
  step2.close();

}

double step2Error(std::vector<std::vector<double>> newPoints, std::vector<double> a, struct Multivariate_RemezIIParameters remezdesc, string path){
  double errorStep = 0;
  createModelStep2Error(a, remezdesc, path);
  string fullpath = path + "models/step2_model_error.txt";
  System sys(fullpath.c_str());
  DefaultOptimizer optimizer(sys,1e-09, 1e-10);
  optimizer.optimize(sys.box);
  optimizer.report();
  errorStep = optimizer.get_loup();
  return -errorStep;
}

void createModelStep2Error(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path){
  ofstream step2(path + "models/step2_model_error.txt");
  step2 << setprecision(13);
  step2 << "Variables\n";
  for(int i = 0; i<remezdesc.fdesc.bornersVar.size(); i++){
    step2 << "  x"<< i+1 <<" in [" << remezdesc.fdesc.bornersVar[i][0] <<" , " << remezdesc.fdesc.bornersVar[i][1] <<"];\n";
  }
  step2 << "\nMinimize\n";
  step2 << "-abs(";
  step2 << remezdesc.fdesc.functionString << "-(";
  std::vector<std::vector<int>> p = remezdesc.poly.derivatedPhi(remezdesc.fdesc.bornersVar.size(), remezdesc.approximationDegree);
  for(int i = 0; i < A.size(); i++){
    step2 << "+(" << A[i] << ")" ;
    for(int j = 0; j < p[i].size(); j++){
      if(p[i][j] != 0){
        step2 << "*(x" << j+1 << "^(" << p[i][j] << "))";
      }
    }
  }
  step2 << "));";
  step2.close(); 
}

bool cheking_side_not_constant(struct Multivariate_RemezIIParameters remezdesc, int side, int var){
  bool doSide = remezdesc.fdesc.notConstantDerived[side];
  if(!remezdesc.fdesc.notConstantDerived[side]){
    doSide = false;
    std::vector<std::vector<int>> p = remezdesc.poly.derivatedPhi(remezdesc.fdesc.bornersVar.size(), remezdesc.approximationDegree);
    for(int i = 0; i < p.size(); i++){
      if(p[i][(var+1)%2] > 1 && remezdesc.poly.a[i] != 0){
        doSide = true;
      }
    }
  }
  
  return doSide;
}

/*-----------------------------------------------------------------*/
/* OUTPUT */
/*-----------------------------------------------------------------*/

string createPath(){
  string path = "../resultsOfSteps/multi/";
  return path;
}

void write_introduction(string path, struct Multivariate_RemezIIParameters remezdesc){
  std::ofstream summary;
  summary.open(path + "summary/" + remezdesc.fdesc.functionString + ".txt",std::ofstream::trunc);
  summary << "Welcome to the summary of ";
  summary << "Remez 1 multi ";
  summary << "for multivariate polynomial approximation. \n";
  summary << "With the following parameters : \n";
  summary << "degree of approximation = " << remezdesc.approximationDegree << "\n";
  summary << "maximum number of turns = " << remezdesc.maxNbTurns << "\n";
  summary << "approximation of error = " << remezdesc.approximationResult << "\n";
  summary << "approximation of points = " << remezdesc.approximationPoints << "\n";
  summary << "x in = [";
  for(int i = 0; i<remezdesc.fdesc.bornersVar.size(); i++){
      summary << "[" << remezdesc.fdesc.bornersVar[i][0] << "," << remezdesc.fdesc.bornersVar[i][1] << "]";
  }
  summary << "]\n";
  summary << "For the function " << remezdesc.fdesc.functionString << ".\n";
}

void writeD(struct Multivariate_RemezIIParameters remezdesc, std::vector<std::vector<double>> D, string path){
  std::ofstream summary;
  summary.open(path + "summary/" + remezdesc.fdesc.functionString + ".txt", std::ios_base::app);
  summary << "\nGenerated D0 is : ";
  for(int i=0; i<D.size(); i++){
    summary << "[";
    for(int j=0; j<D[i].size(); j++){
      summary << D[i][j] << ", ";
    }
    summary << "] ";
  }
}

void writePolynom(std::vector<double> A, string path, struct Multivariate_RemezIIParameters remezdesc){
  std::ofstream summary;
  summary.open(path + "summary/" + remezdesc.fdesc.functionString + ".txt", std::ios_base::app);
  summary << "Polynom = [" << A[0] ;
  for(int j = 1; j< A.size(); j++){
    summary << "," << setprecision(13) << A[j] ;
  }
  summary << "] \n" ;
}

void writeErrors(double error1, double error2, struct Multivariate_RemezIIParameters remezdesc, string path){
  std::ofstream summary;
  summary.open(path + "summary/" + remezdesc.fdesc.functionString + ".txt", std::ios_base::app);
  summary << "Error step 1 : " << setprecision(13) << error1;
  summary << "\nError step 2 : " << setprecision(13) << error2 ;
  summary << "\nDistance between errors :" << setprecision(13) << error2-error1 << "\n";
}

void write_data_for_graphs(std::vector<double> A, std::vector<std::vector<double>> bornersVar, int turn){
  std::ofstream data;
  data << setprecision(13) ;
  data.open("data/"+ to_string(turn) +"data.txt",std::ofstream::trunc);
  for(int i = 0; i<A.size(); i++){
    data << A[i] << " ";
  }
  //for ease, we put 0 until having the same number of "case" for the python tab
  for(int i = A.size()-1 ; i<=bornersVar.size(); i++){
    data << "0 ";
  }
  
  data << "\n";
  for(int i = 0; i<bornersVar.size(); i++){
    data << bornersVar[i][0] << " " << bornersVar[i][1] << " ";
  }
  for(int i = bornersVar.size()+2; i<A.size(); i++){
    data << "0" << " ";
  }
}

std::vector<std::vector<double>> noCopyAdd(std::vector<std::vector<double>> v1, std::vector<std::vector<double>> v2){
  bool toNotAdd = true;
  std::vector<std::vector<double>> v = v1;
  
  v.insert(v.end(), v2.begin(), v2.end());
  std::sort(v.begin(), v.end());
  // Remove duplicate values from vector
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return v;
}

/*-----------------------------------------------------------------*/
/* GRAPHS */
/*-----------------------------------------------------------------*/

void createErrorGraph(struct Multivariate_RemezIIParameters remezdesc, string path){
     //TO DEBUG
  Vec x0 = linspace(remezdesc.fdesc.bornersVar[0][0], remezdesc.fdesc.bornersVar[0][1], 100);
  Vec x1 = linspace(remezdesc.fdesc.bornersVar[1][0], remezdesc.fdesc.bornersVar[1][1], 100);
  std::vector<double> e={};
  std::vector<double> X0={};
  std::vector<double> X1={};
  std::vector<double> x = {};
  int i = 0;
  for (int i = 0; i<x0.size(); i++) {
    for (int j=0; j<x1.size(); j++) {
      X0.push_back(x0[i]);
      X1.push_back(x1[j]);
      x.push_back(x0[i]);
      x.push_back(x1[j]);
      e.push_back(error(remezdesc.poly, x, remezdesc.fdesc.f));
      x = {};
    }
  }
  Plot3D plot;
  plot.xlabel("x0");
  plot.ylabel("x1");
  plot.zlabel("error");
  
  plot.border().clear();
  plot.border().bottomLeftFront();
  plot.border().bottomRightFront();
  plot.border().leftVertical();
  
  plot.palette("dark2");
  plot.drawCurve(X0, X1, e).label("error");
  Figure fig = {{plot}};
  Canvas canvas = {{fig}};
  canvas.show();
  canvas.save(path + "errorEndStep1.pdf");
}



void grapheConvergenceComparisons(struct Multivariate_RemezIIResult remezResults){
  string graphName = "../resultsOfSteps/multi/GrapheConvergenceAll.pdf";
  string output = "";
  std::vector<double> distanceTemp = {};
  std::vector<std::vector<double>> distance = {};
  double dist = 0;
  std::vector<double> conv = {};
    for(int i = 0; i<remezResults.errorStep1.size(); i++){
      dist = sqrt(pow((remezResults.errorStep2[i]-remezResults.errorStep1[i]),2));
      distanceTemp.push_back(dist);
      conv.push_back(dist);
    }
    distance.push_back(distanceTemp);
    distanceTemp = {};
  
  int maxNbTurns = remezResults.nbTurns;

  vector<vector<double>> x = {};
  vector<double> xi = {};
    for(int j=1; j<remezResults.nbTurns+1; j++){
      xi.push_back(j);
    }
    x.push_back(xi);

  Plot2D plot;
  plot.xlabel("Turn");
  plot.ylabel("Convergence");
  plot.xrange(1., maxNbTurns-1);
  plot.yrange(*min_element(conv.begin(), conv.end()), *max_element(conv.begin(), conv.end()));

  plot.legend()
    .atOutsideBottom()
    .displayHorizontal()
    .displayExpandWidthBy(2);
  plot.ytics().logscale(2);
  string curveTitle = "convergence logarithmic Remez";
  plot.drawCurve(x[0], distance[0])
        .label(curveTitle);
  Figure fig = {{plot}};
  Canvas canvas = {{fig}};
  canvas.save(graphName);
}












