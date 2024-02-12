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
#include <gmp.h>

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
  std::vector<string> (*derivatedPhi)(std::vector<double> a, int nbX, int degree);
  int degree;
} MultivariatepolynomialDescription;

//Description of the function to approximate
struct Multivariate_functionDescription {
  string functionString; // With the unknown x1 and x2
  int nbX;
  std::vector<string> derivedFunctionStrings;
  std::vector<std::vector<double>> bornersVar;
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
struct Multivariate_functionDescription initialize_fdesc(string fString, std::vector<string> df, double (*f)(std::vector<double> x), std::vector<std::vector<double>> bornersVar);
struct Multivariate_RemezIIParameters initialize_remezdesc(struct Multivariate_functionDescription fdesc, int degree, int nbTurns);
struct Multivariate_Polynomiale initialize_poly(std::vector<double>(*phi)(std::vector<double> x, int degree), std::vector<string> (*derivatedPhi)(std::vector<double> a, int nbX, int degree));
std::vector<std::vector<double>> step0(int nbX, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_functionDescription fdesc);
string createPath();
std::vector<double> step1(std::vector<std::vector<double>> D, struct Multivariate_functionDescription fdesc, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_Polynomiale poly);
void write_introduction(string path, struct Multivariate_RemezIIParameters remezdesc);
void writeD(struct Multivariate_RemezIIParameters remezdesc, std::vector<std::vector<double>> D, string path);
double step2(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path);
void createErrorGraph(struct Multivariate_RemezIIParameters remezdesc, string path);
void writeA(std::vector<double> a, string path, string functionString);
std::vector<std::vector<double>> step0_random(int nbPoints, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_functionDescription fdesc);
void write_data_for_graphs(std::vector<double> A, std::vector<std::vector<double>> bornersVar);
/*-----------------------------------------------------------------*/
/* FUNCTIONS TO APPROXIMATE */
/*-----------------------------------------------------------------*/
//Didactical Example

double f0(std::vector<double> x) {
  return pow(x[0],2) + pow(x[1],2);
}

struct Multivariate_functionDescription initializeF0(){
  return initialize_fdesc("x1^2 + x2^2",{"2*x1","2*x2", "0"},f0, {{0,1},{0,1}});
}

//Reemtsen Examples
double f1(std::vector<double> x) {
  return log(x[0]+x[1])*sin(x[0]);
}

struct Multivariate_functionDescription initializeF1(){
  return initialize_fdesc("log(x1+x2)*sin(x1)",{"sin(x1)+log(x1+x2)*cos(x1)*(x1+x2)","sin(x1)/log(x1+x2)"},f1, {{0,1},{1, 2.5}});
}

double f2(std::vector<double> x) {
  return pow((1+x[0]), x[1]);
}
struct Multivariate_functionDescription initializeF2(){
  return initialize_fdesc("(1+x1)^(x2)",{"x2*(1+x1)^(x2-1)","(1+x1)^x2*ln(1+x1)"},f2, {{0,1},{1, 2.5}});
}

double f3(std::vector<double> x) {
  return cos(x[2])*pow(1+x[0],x[1]);
}
struct Multivariate_functionDescription initializeF3(){
  return initialize_fdesc("cos(x3)*(1+x1)^(x2)",{"x2*cos(x3)*(1+x1)^(x2-1)"," cos(x3)*(1+x1)^(x2)*ln(1+x1)","(-sin(x3)*(1+x1)^(x2))"},f3, {{0,1},{1, 2},{0,1}});
}

double f4(std::vector<double> x) {
  return 1/(x[0]+2*x[1]+4);
}
struct Multivariate_functionDescription initializeF4(){
  return initialize_fdesc("1/(x1+2*x2+4)",{"-1/((x1+2*x2+4)^2)","-1/((x1+2*x2+4)^2)"},f4, {{-1,1},{-1,1}});
}

double f5(std::vector<double> x) {
  return exp(pow(x[0],2)+x[0]*x[1]);
}
struct Multivariate_functionDescription initializeF5(){
  return initialize_fdesc("exp(x1^2 + x1*x2)",{"exp(x1^2 + x1*x2)*(2x+x2)", "exp(x1^2 + x1*x2)*(x)"}, f5, {{-1,1},{-1,1}});
}

double f6_7(std::vector<double> x) {
  return sqrt(x[0]+2*x[1]+4);
}
struct Multivariate_functionDescription initializeF6_7(){
  return initialize_fdesc("sqrt(x1+2*x2+4)",{"1/2*sqrt(x1+2*x2+4)","1/sqrt(x1+2*x2+4)"},f6_7, {{-1,1},{-1,1}});
}

double f8(std::vector<double> x) {
  return abs(log((x[0]*x[1]+1)/(x[0]+0.5)))*pow(x[1],(x[2]+1)/2);
}
struct Multivariate_functionDescription initializeF8(){
  return initialize_fdesc("",{},f8, {{0,1},{0,1},{0,1}});
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
  std::vector<double> p = {1};
  std::vector<double> temp = {};
  std::vector<double> prev = {1};
  for(int i = 1; i<=degree; i=i+1){
    for(int j = 0; j<x.size(); j=j+1){
      for(int k = 0; k<prev.size(); k++){
        temp.push_back(prev[k]*x[j]);
      }
    } 
    for(int j = 0; j<temp.size(); j++){
      p.push_back(temp[j]);
    }
    prev = temp;
    temp = {};
  }
  return p;
}

std::vector<string> phiBasiqueDerivated(std::vector<double> a, int nbX, int degree){
  std::vector<std::vector<double>> powOfX = {{1}};
  std::vector<std::vector<double>> temp = {std::vector<double>(nbX, 0.0)};
  std::vector<std::vector<double>>  prev = {{1}};
  
  for(int i = 1; i<degree; i=i+1){
   for(int j = 0; j<nbX; j=j+1){
      for(int k = 0; k<prev.size(); k++){
	temp.push_back(prev[k]);
	temp[temp.size()-1][j] += 1; 
      }
    }
    for(int j = 0; j<temp.size(); j++){
      powOfX.push_back(temp[j]);
    }
    prev = temp;
    std::fill(temp[0].begin(), temp[0].end(), 0.0);
    temp = {temp[0]};
  }
  string str;
  std::vector<string> result;
  for(int k = 0; k<nbX; k=k+1){
    str = "";
    for(int i = 0; i<powOfX.size(); i++){
      if(powOfX[i][k] >= 1){
        str = "2*" ;
        str += a[i];
        for(int j = 0; j<powOfX[i].size(); j++){
          if(powOfX[i][j]!=0){
            str += "x" + (j+1);
          }
        }
      }
      result.push_back(str);
      str = "";
    }
  }
  
  return result;
  
}

struct Multivariate_Polynomiale initializePhiBasic(){
  return initialize_poly(phiBasique, phiBasiqueDerivated);
}

/*-----------------------------------------------------------------*/
/* INITIALIZES OBJECTS */
/*-----------------------------------------------------------------*/
struct Multivariate_functionDescription initialize_fdesc(string fString, std::vector<string> df, double (*f)(std::vector<double> x), std::vector<std::vector<double>> bornersVar){
  struct Multivariate_functionDescription fdesc;
  fdesc.functionString = fString;
  fdesc.derivedFunctionStrings = df;
  fdesc.f = f;
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

struct Multivariate_Polynomiale initialize_poly(std::vector<double>(*phi)(std::vector<double> x, int degree), std::vector<string> (*derivatedPhi)(std::vector<double> a, int nbX, int degree)){
  struct Multivariate_Polynomiale poly;
  poly.a = {};
  poly.phi = (*phi);
  poly.derivatedPhi = (*derivatedPhi);
  return poly;
}

/*-----------------------------------------------------------------*/
/* MAIN */
/*-----------------------------------------------------------------*/

int main(int argc, char** argv) { 
  int nbTurns = 100;
  int degree = 1;
  double errorStep1 = 0;
  string path = createPath();
  struct Multivariate_RemezIIParameters remezdesc;
  struct Multivariate_functionDescription fdesc;
  //struct Multivariate_Polynomiale poly;
  
  fdesc = initializeF0();
  remezdesc = initialize_remezdesc(fdesc, degree, nbTurns);
  remezdesc.poly = initializePhiBasic();
  remezdesc.poly.degree = degree;
  write_introduction(path, remezdesc);
  
  
  //std::vector<std::vector<double>> D = step0(fdesc.nbX, remezdesc, fdesc);
  int nbPoints = 1000;
  std::vector<std::vector<double>> D = step0_random(nbPoints, remezdesc, fdesc);
  
  
  writeD(remezdesc, D, path);
  remezdesc.poly.a = step1(D, fdesc, remezdesc, remezdesc.poly);
  
  //Taking the error out of the polynomial
  errorStep1 = remezdesc.poly.a[remezdesc.poly.a.size()-1];
  remezdesc.poly.a.pop_back();
  
  //writeA(remezdesc.poly.a, path, fdesc.functionString);
  createErrorGraph(remezdesc, path);
  write_data_for_graphs(remezdesc.poly.a, fdesc.bornersVar);
  //double newPoints = step2(remezdesc.poly.a, remezdesc, path);
}

/*-----------------------------------------------------------------*/
/* REMEZ II */
/*-----------------------------------------------------------------*/

std::vector<std::vector<double>> step0(int nbX, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_functionDescription fdesc){
  std::vector<std::vector<double>> D0 = {};
  std::vector<double> d = {};
  for(int i = 0; i < nbX; i ++){ //middle
    d.push_back((fdesc.bornersVar[i][1] - fdesc.bornersVar[i][0])/2);
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
/*
std::vector<std::vector<double>> step2(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path){
  std::vector<std::vector<double>> pointsSide = {};
  std::vector<std::vector<double>> pointsCenter = {};
  //Getting the points
  pointsSide = step2_side(struct Multivariate_RemezIIParameters remezdesc, string path);
  pointsCenter = step2_center(struct Multivariate_RemezIIParameters remezdesc, string path);
  //Getting every points together
  pointsCenter.insert( pointsCenter.end(), pointsSide.begin(), pointsSide.end() );
 
  return pointsCenter; 
}

std::vector<std::vector<double>> step2_side(struct Multivariate_RemezIIParameters remezdesc, string path){
    
}

std::vector<std::vector<double>> step2_center(struct Multivariate_RemezIIParameters remezdesc, string path){

  Function f("x","y","z","x*y*z");
    Function df(f,Function::DIFF);
  cout << "df=" << df << endl;

}
*/

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
  summary << "Generated D0 is : ";
  for(int i=0; i<D.size(); i++){
    summary << "[";
    for(int j=0; j<D[i].size(); j++){
      summary << D[i][j] << ", ";
    }
    summary << "] ";
  }
}

void write_data_for_graphs(std::vector<double> A, std::vector<std::vector<double>> bornersVar){
  std::ofstream data;
  data.open("data/data.txt",std::ofstream::trunc);
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















