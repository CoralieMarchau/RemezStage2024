#include "ibex.h"
#include <glpk.h>
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
#include<string>

using namespace std;
using namespace ibex;
using namespace soplex;
using namespace sciplot;

/*-----------------------------STRUCT------------------------------------------------*/











/*----------------------FUNCTION DECLARATION-----------------------------------------*/
int RemezUniSoplexExchange(int nbA, double approx, int nbTurns, double up, double down, string function);
std::vector<double> exchange(std::vector<double> D, std::vector<double> A, double newPoint);
void makeGraphUni(std::vector<double> P, string function, int turn, int degree, std::vector<double> borners, string path);
void grapheConvergenceRemez1(std::vector<double> errorStep1, std::vector<double> errorStep2, int nbTurns, string function, string path);
std::vector<int> getIndice(int nbA, int whereWeAre);
void step2PlusSop(std::vector<double> P, int nbA, string derivedFunction, double up, double down);
int RemezUniPlusSoplex(int nbA, double approx, int sizeD, int nbTurns, double up, double down, string function, string derivedFunction);
std::vector<double> step1Soplex(std::vector<double> D, int sizeD, int nbA);
void createMPSFileHorner(std::vector<double> D, int sizeD, int nbA);
std::vector<double> step0Univarie(int sizeD, double up, double low);
void step1(std::vector<double> D, int sizeD, int nbA);
void createMPSFile(std::vector<double> D, int sizeD, int nbA);
void createModelStep2Sop(std::vector<double> P, int nbA, string function, double up, double down);
void createModelStep2PlusSop(std::vector<double> P, int nbA, string derivedFunction, double up, double down);
void step2Sop(std::vector<double> P, int nbA, string function, double up, double down);
void createModelStep2(std::vector<string> P, int nbA, string function, double up, double down);
void step2(std::vector<string> P, int nbA, string function, double up, double down);
std::vector<string> parserResultStep1(int nbA);
std::vector<double> parserResultStep2();
bool tooClose(std::vector<double> D, int sizeD, double newPoint, double approx);
void createModelStep2Plus(std::vector<string> P, int nbA, string derivedFunction, double up, double down);
void step2Plus(std::vector<string> P, int nbA, string derivedFunction, double up, double down);
std::vector<std::vector<double>> parserResultStep2Plus();
int RemezUni(int nbA, double approx, int sizeD, int nbTurns, double up, double down, string function);
int RemezUniPlus();
void convergenceCalculation(std::vector<double> errorStep1, std::vector<double> errorStep2, string path);
std::vector<std::vector<double>> parserResultStep2Plus();
void createMPSFileMultivariate(std::vector<std::vector<double>> D, int sizeD, int nbVar, int nbA);
std::vector<std::vector<double>> step0_multi(int sizeD, int nbVar, std::vector<std::vector<double>> borners);
void step1_multivariate(std::vector<std::vector<double>> D, int sizeD, int nbVar, int nbA);
std::vector<string> parserResultStep1_multivariate(int nbA, int nbVar);
void createModelStep2_multivariate(std::vector<string> P, int nbA, int nbVar, string function, std::vector<std::vector<double>> borners);
void step2_multi(std::vector<string> P, int nbA, int nbVar, string function, std::vector<std::vector<double>> borners);
int RemezMulti();
std::vector<double>  parserResultStep2Multi(int nbVar);
std::vector<bool> tooClosePlus(std::vector<std::vector<double>> newPoints, std::vector<double> D, double approx);
bool tooCloseMulti(std::vector<std::vector<double>> D, int nbVar, int sizeD, std::vector<double> newPoint, double approx);
int RemezUniPlus(int nbA, double approx, int sizeD, int nbTurns, double up, double down, string function, string derivedFunction);
int RemezUniSoplex(int nbA, double approx, int sizeD, int nbTurns, double up, double down, string function);
void makeGraphUniPlus(std::vector<double> P, string function, int turn, std::vector<double> borners, double error, string path);
double Poly(double x, std::vector<double> A);

/*------------------------------MAIN------------------------------------------------*/
//Pour convergence monovariée : exemple cos(x) + sin(5x)
//Pour convergence multivariée ; exemple goutte, un cas est singulier (la tête) et pas l'autre (la queue).

double f(double x) {
  return cos(x) + sin(5*x);
  //return cos(x) + 2*(sin(x)) + x - 1;
}

double Poly(double x, std::vector<double> A) {
  double y = 0;
  for(int i = 0 ; i < A.size(); i++){
    y += A[i] * pow(x, i);
  }
  return y;
}

int main(int argc, char** argv) { 
  int nbA = 10;
  //double approx = 0.0000005;
  double approx = 1.e-10;
  int sizeD = 20;
  int nbTurns = 100;
  double up = 1;
  double down = -1;
  string function = "cos(x)+sin(5*x)";
  string derivedFunction = "-sin(x)+5*cos(5*x)";
  //string function = "cos(x) + 2*(sin(x)) + x - 1";
  //string derivedFunction = "-sin(x) + 1 + 2*cos(x)";

  RemezUniSoplexExchange( nbA,  approx,  nbTurns,  up,  down,  function);
  RemezUniSoplex(nbA, approx, sizeD, nbTurns, up, down, function);
  RemezUniPlusSoplex(nbA, approx, sizeD, nbTurns, up, down, function, derivedFunction);
  //RemezMulti();
  
  //RemezUni(nbA, approx, sizeD, nbTurns, up, down, function);  
  //RemezUniPlus(nbA, approx, sizeD, nbTurns, up, down, function, derivedFunction);
  
}


void convergenceCalculation(std::vector<double> errorStep1, std::vector<double> errorStep2, string path){
  ofstream convText(path+"convergence.txt");
  string output = "";
  std::vector<double> distance = {};
  double dist = 0;
  std::vector<double> conv = {};
  for(int i = 0; i<errorStep1.size(); i++){
    dist = sqrt(pow((errorStep2[i]-errorStep1[i]),2));
    convText << "distance au tour " << i << " : " << dist << "\n";
    distance.push_back(dist);
  }
  conv = distance;
  
  vector<double> x = {};
  for(int i=1; i<conv.size(); i++){
    x.push_back(i);
  }
  Plot2D plot;
  plot.xlabel("Turn");
  plot.ylabel("Convergence");
  plot.xrange(0., conv.size()-1);
  plot.yrange(*min_element(conv.begin(), conv.end()), *max_element(conv.begin(), conv.end()));

  plot.legend()
    .atOutsideBottom()
    .displayHorizontal()
    .displayExpandWidthBy(2);
  plot.ytics().logscale(2);
  plot.drawCurve(x, conv)
        .label("convergence");
  Figure fig = {{plot}};
  Canvas canvas = {{fig}};
  canvas.save(path + "convergenceLog.pdf");
}

void grapheConvergenceRemez1(std::vector<double> errorStep1, std::vector<double> errorStep2, int nbTurns, string function, string path){
  string graphName = function + "GrapheConvergence";
  graphName = path + graphName;
  graphName += ".pdf";
  Vec x = linspace(1, nbTurns, nbTurns);
  Plot2D plot;
  plot.xlabel("x");
  plot.ylabel("y");
  plot.xrange(1, nbTurns);
  plot.yrange(errorStep1[0] , errorStep2[0]);
  plot.legend()
    .atOutsideBottom()
    .displayHorizontal()
    .displayExpandWidthBy(2);
  plot.drawCurve(x, errorStep1).label("error step 1");
  plot.drawCurve(x, errorStep2).label("error step 2");
  Figure fig = {{plot}};
  Canvas canvas = {{fig}};
  canvas.save(graphName);
}


void makeGraphUni(std::vector<double> P, string function, int turn, int degree, std::vector<double> borners, double error, string path){
  string graphName = function + to_string(turn) + "AproximationRemez1";
  graphName = path + graphName;
  graphName += ".pdf";
  Vec x = linspace(borners[0], borners[1], 1000);
  
  std::vector<double> e = {};
  
  for(int i=0; i<x.size(); i++){
    e.push_back(Poly(x[i],P)-f(x[i]));
  }
  
  Plot2D plot;
  plot.xlabel("x");
  plot.ylabel("y");
  plot.xrange(borners[0], borners[1]);
  plot.yrange(-error , error);
  
  plot.legend()
    .atOutsideBottom()
    .displayHorizontal()
    .displayExpandWidthBy(2);
  plot.drawCurve(x, e).label("error");
  Figure fig = {{plot}};
  Canvas canvas = {{fig}};
  canvas.save(graphName);
}


void makeGraphUniPlus(std::vector<double> P, string function, int turn, std::vector<double> borners, double error, string path){
  string graphName = function + to_string(turn) + "AproximationRemez1+";
  graphName = path + graphName;
  graphName += ".pdf";
  Vec x = linspace(borners[0], borners[1], 1000);
  
  std::vector<double> e = {};
  
  for(int i=0; i<x.size(); i++){
    e.push_back(Poly(x[i],P)-f(x[i]));
  }
  
  Plot2D plot;
  plot.xlabel("x");
  plot.ylabel("y");
  plot.xrange(borners[0], borners[1]);
  plot.yrange(-error , error);
  
  plot.legend()
    .atOutsideBottom()
    .displayHorizontal()
    .displayExpandWidthBy(2);
  plot.drawCurve(x, e).label("error");
  Figure fig = {{plot}};
  Canvas canvas = {{fig}};
  canvas.save(graphName);
}

double fMulti(std::vector<double> x) {
  return cos(x[0]) + sin(x[1]);
}

/*---------------------------FUNCTIONS----------------------------------------------*/
/*---------------------------REMEZ 1 UNIVARIE---------------------------------------*/

int RemezUniSoplexExchange(int nbA, double approx, int nbTurns, double up, double down, string function){
  ofstream summary("../resultsOfSteps/summary_of_executions/Summary_Remez_1_Exchange");
  summary << "Welcome to the summary of Remez 1 for univariate polynomial approximation. \n";
  summary << "The parameters of Remez have been initialized as : \n";
  double error, newPoint;
  std::vector<double> errorStep = {};
  int sizeD = nbA + 2;
  double errorStep1 = 0.0;
  int i = 1;
  string path = "../resultsOfSteps/img/Remez1Exchange/";
  bool itIsTheEnd = false;
  std::vector<double> P ;
  std::vector<double> D = step0Univarie(sizeD, up, down);
  double temp;
  temp = D[1];
  D[1] = D[0];
  D.erase(D.begin());
  D.push_back(temp);
  std::vector<double> errorStep2 = {};
  
  summary << "-Degree of the polynomial approximation : " << nbA + 1 << "\n";
  summary << "-Size of the first discretized set : " << sizeD << "\n";
  summary << "-Approximation of x : " << approx << "\n";
  summary << "-Function to approximate : " << function << "\n";
  summary << "-Maximum number of turn allowed : " << nbTurns << "\n";
  summary << "-x borners : [" << up << " , " << down << "]\n";
  summary << "BEGINNING : \n" ;
  while(!itIsTheEnd){
    summary << "\n Turn " << i << "\n";
    summary << "D"<<i<<"=[" << D[0];
    for(int j = 1; j<sizeD; j++){
      summary << ","  << D[j];
    }
    summary << "]\n";
    P = step1Soplex(D, sizeD, nbA);
    summary << "[" << P[0] ;
    for(int j = 1; j<=nbA; j++){
      summary << "," << setprecision(13) << P[j] ;
    }
    summary << "] \n" ;
    for(int j = 0; j<=nbA; j++){
      summary << "-a" << j << "= " << P[j] << "\n";
    }
    errorStep1 = P[nbA+1];
    P.pop_back();
    step2Sop(P, nbA, function, up, down);
    std::vector<double> step2Result = parserResultStep2();
    error = -step2Result[0];
    newPoint = step2Result[1];
    errorStep.push_back(errorStep1);
    errorStep2.push_back(error);
    summary << "\nError : " << error << "\n" ;
    summary << "Error for x in D" << i-1 << ": " << errorStep1 << "\n" ;
    if (tooClose(D, sizeD, newPoint, approx) || i >= nbTurns){
      itIsTheEnd = true ;
      summary << "\n Last error is : " << setprecision(13) << errorStep[i-1] << " \n";
      summary << "\n END \n";
    } else {
        D = exchange(D, P, newPoint);
        summary << "Added x is : " << newPoint << "\n";
    }
     itIsTheEnd = itIsTheEnd || (errorStep[errorStep.size()-1] >= errorStep2[errorStep.size()-1]-approx && errorStep[errorStep.size()-1] <= errorStep2[errorStep.size()-1]+approx);
    std::vector<double> borners = {};
    borners.push_back(down);
    borners.push_back(up);
    makeGraphUni(P, function, i, nbA, borners, error, path);
    i = i + 1;
  }
  convergenceCalculation(errorStep, errorStep2,path);
  grapheConvergenceRemez1(errorStep, errorStep2, i-1, function, path);
  summary << "\n Remez finished after " << i-1 << " turns. \n";
  return i;
}

std::vector<double> exchange(std::vector<double> D, std::vector<double> A, double newPoint){
  std::vector<double> Di = {};
  int i = 0;
  while(!(newPoint >= D[i] && newPoint <= D[i+1])){
    Di.push_back(D[i]);
    i = i + 1;
  }
  double error1, error2, errorNP;
  errorNP = Poly(newPoint,A) - f(newPoint);
  error1 = Poly(D[i],A) - f(D[i]);
  error2 = Poly(D[i+1],A) - f(D[i+1]);
  if(error1<0 && errorNP<0 || error1>0 && errorNP>0 ){
    Di.push_back(newPoint);
    Di.push_back(D[i+1]);
  } else {
  Di.push_back(D[i]);
  Di.push_back(newPoint); }
  i = i + 2;
  for(int j = i; j< D.size(); j++){
    Di.push_back(D[j]);
  }
  return Di;
}

int RemezUni(int nbA, double approx, int sizeD, int nbTurns, double up, double down, string function){
  ofstream summary("../resultsOfSteps/summary_of_executions/Summary_Remez_1");
  summary << "Welcome to the summary of Remez 1 for univariate polynomial approximation. \n";
  summary << "The parameters of Remez have been initialized as : \n";
  double error, newPoint;
  std::vector<double> errorStep = {};
  std::vector<double> errorStep2 = {};
  double errorStep1 = 0.0;
  int i = 1;
  bool itIsTheEnd = false;
  string path = "../resultsOfSteps/img/Remez1/glpk/";
  double a;
  std::vector<string> P ;
  std::vector<double> Pdouble ;
  std::vector<double> D = step0Univarie(sizeD, up, down);
  
  summary << "-Degree of the polynomial approximation : " << nbA + 1 << "\n";
  summary << "-Size of the first discretized set : " << sizeD << "\n";
  summary << "-Approximation of x : " << approx << "\n";
  summary << "-Function to approximate : " << function << "\n";
  summary << "-Maximum number of turn allowed : " << nbTurns << "\n";
  summary << "-x borners : [" << up << " , " << down << "]\n";
  summary << "BEGINNING : \n" ;
  
  while(!itIsTheEnd){
    summary << "\n Turn " << i << "\n";
    summary << "D"<<i<<"=[" << D[0];
    for(int j = 1; j<sizeD; j++){
      summary << ","  << D[j];
    }
    summary << "]\n";
    step1(D, sizeD, nbA);
    P = parserResultStep1(nbA);
    summary << "[" << P[0] ;
    istringstream  istr(P[0]);
    istr >> a;
    Pdouble.push_back(a);
    for(int j = 1; j<=nbA; j++){
      summary << "," << P[j] ;
      istringstream  istr(P[j]);
      istr >> a;
      Pdouble.push_back(a);
    }
    summary << "] \n" ;
    for(int j = 0; j<=nbA; j++){
      summary << "-a" << j << "= " << P[j] << "\n";
    }
    istringstream  istr1(P[nbA+1]);
    istr1 >> errorStep1;
    P.pop_back();
    step2(P, nbA, function, up, down);
    std::vector<double> step2Result = parserResultStep2();
    error = -step2Result[0];
    newPoint = step2Result[1];
    errorStep.push_back(errorStep1);
    errorStep2.push_back(error);
    summary << "\nError : " << error << "\n" ;
    summary << "Error for x in D" << i-1 << ": " << setprecision(13)  << errorStep1 << "\n" ;
    if (tooClose(D, sizeD, newPoint, approx) || i >= nbTurns){
      itIsTheEnd = true ;
      summary << "\n Last error is : " << setprecision(13)  <<  errorStep[i-1] << " \n";
      summary << "\n END \n";
    } else {
        D.push_back(newPoint);
        sizeD = sizeD + 1;
        summary << "Added x is : " << newPoint << "\n";
    }
    itIsTheEnd = itIsTheEnd || (errorStep[errorStep.size()-1] >= errorStep2[errorStep.size()-1]-approx && errorStep[errorStep.size()-1] <= errorStep2[errorStep.size()-1]+approx);
    std::vector<double> borners = {};
    borners.push_back(down);
    borners.push_back(up);
    makeGraphUni(Pdouble, function, i, nbA, borners, error, path);
    i = i + 1;
    Pdouble = {};
  }
  convergenceCalculation(errorStep, errorStep2,path);
  grapheConvergenceRemez1(errorStep, errorStep2, i-1, function, path);

  summary << "\n Remez finished after " << i-1 << " turns. \n";
  return i;
}


int RemezUniSoplex(int nbA, double approx, int sizeD, int nbTurns, double up, double down, string function){
  ofstream summary("../resultsOfSteps/summary_of_executions/Summary_Remez_1_usingSOPLEX");
  summary << "Welcome to the summary of Remez 1 for univariate polynomial approximation. \n";
  summary << "The parameters of Remez have been initialized as : \n";
  double error, newPoint;
  std::vector<double> errorStep = {};
  double errorStep1 = 0.0;
  int i = 1;
  string path = "../resultsOfSteps/img/Remez1/soplex/";
  bool itIsTheEnd = false;
  std::vector<double> P ;
  std::vector<double> D = step0Univarie(sizeD, up, down);
  std::vector<double> errorStep2 = {};
  
  summary << "-Degree of the polynomial approximation : " << nbA + 1 << "\n";
  summary << "-Size of the first discretized set : " << sizeD << "\n";
  summary << "-Approximation of x : " << approx << "\n";
  summary << "-Function to approximate : " << function << "\n";
  summary << "-Maximum number of turn allowed : " << nbTurns << "\n";
  summary << "-x borners : [" << up << " , " << down << "]\n";
  summary << "BEGINNING : \n" ;
  while(!itIsTheEnd){
    summary << "\n Turn " << i << "\n";
    summary << "D"<<i<<"=[" << D[0];
    for(int j = 1; j<sizeD; j++){
      summary << ","  << D[j];
    }
    summary << "]\n";
    P = step1Soplex(D, sizeD, nbA);
    summary << "[" << P[0] ;
    for(int j = 1; j<=nbA; j++){
      summary << "," << setprecision(13) << P[j] ;
    }
    summary << "] \n" ;
    for(int j = 0; j<=nbA; j++){
      summary << "-a" << j << "= " << P[j] << "\n";
    }
    errorStep1 = P[nbA+1];
    P.pop_back();
    step2Sop(P, nbA, function, up, down);
    std::vector<double> step2Result = parserResultStep2();
    error = -step2Result[0];
    newPoint = step2Result[1];
    errorStep.push_back(errorStep1);
    errorStep2.push_back(error);
    summary << "\nError : " << error << "\n" ;
    summary << "Error for x in D" << i-1 << ": " << errorStep1 << "\n" ;
    if (tooClose(D, sizeD, newPoint, approx) || i >= nbTurns){
      itIsTheEnd = true ;
      summary << "\n Last error is : " << setprecision(13) << errorStep[i-1] << " \n";
      summary << "\n END \n";
    } else {
        D.push_back(newPoint);
        sizeD = sizeD + 1;
        summary << "Added x is : " << newPoint << "\n";
    }
     itIsTheEnd = itIsTheEnd || (errorStep[errorStep.size()-1] >= errorStep2[errorStep.size()-1]-approx && errorStep[errorStep.size()-1] <= errorStep2[errorStep.size()-1]+approx);
    std::vector<double> borners = {};
    borners.push_back(down);
    borners.push_back(up);
    makeGraphUni(P, function, i, nbA, borners, error, path);
    i = i + 1;
  }
  convergenceCalculation(errorStep, errorStep2,path);
  grapheConvergenceRemez1(errorStep, errorStep2, i-1, function, path);
  summary << "\n Remez finished after " << i-1 << " turns. \n";
  return i;
}

//Tools : 

bool tooClose(std::vector<double> D, int sizeD, double newPoint, double approx){
  bool isClose = false;
  for(int i=0 ; i<sizeD; i++){
    if(D.at(i)>=newPoint-approx && D.at(i)<=newPoint+approx){ isClose = true;
    }
  }
  return isClose; 
}


std::vector<double> step0Univarie(int sizeD, double up, double low){
   double step = (up-low)/(sizeD-1);
   double temp = low;
   std::vector<double> D0 = {low, up};
   for(int i = 0; i < sizeD-2; i++) {
     temp = temp + step;
     D0.push_back(temp);
   }
   return D0;
}

//STEP 1 :


std::vector<double> step1Soplex(std::vector<double> D, int sizeD, int nbA){
   SoPlex mysoplex;

   /* set parameters for exact solving */
   mysoplex.setIntParam(SoPlex::READMODE, SoPlex::READMODE_REAL);
   mysoplex.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_REAL);
   mysoplex.setIntParam(SoPlex::CHECKMODE, SoPlex::CHECKMODE_REAL);
   mysoplex.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
   mysoplex.setRealParam(SoPlex::FEASTOL, 0.0);
   mysoplex.setRealParam(SoPlex::OPTTOL, 0.0);
   
   /* set the objective sense */
   mysoplex.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
   
   /* we first add variables (the integer data is converted to type Rational) */
   DSVectorReal dummycol(0);  
   for(int i = 0; i<=nbA; i++){
     mysoplex.addColReal(LPColReal(0, dummycol, infinity, -infinity));
   }
   mysoplex.addColReal(LPColReal(1, dummycol, infinity, -infinity));
   
   double x;
   Real r;
   /* then constraints one by one (here we show how Rationals can be used directly) */
   for(int i = 0 ; i<D.size(); i++){
     DSVectorReal row1(0);
     for(int j = 0 ;  j<=nbA; j++){
       x = pow(D[i], j);
       r = x;
       row1.add(j, r);
     }
     r = 1;
     row1.add(nbA+1, r);
     r = f(D[i]);
     mysoplex.addRowReal(LPRowReal(r, row1, infinity));
   }
   
   
   for(int i = 0 ; i<D.size(); i++){
     DSVectorReal row1(0);
     for(int j = 0 ;  j<=nbA; j++){
       x = pow(D[i], j);
       r = -x;
       row1.add(j, r);
     }
     r = 1;
     row1.add(nbA+1, r);
     r = -f(D[i]);
     mysoplex.addRowReal(LPRowReal(r, row1, infinity));
   }
   
   mysoplex.writeFileReal("step1Soplex.lp", NULL, NULL, NULL);
   SPxSolver::Status stat;
   DVectorRational prim(nbA+2);
   DVectorRational dual(D.size()*2);
   stat = mysoplex.optimize();
   std::vector<double> primDouble;
   /* get solution */
   if(stat == SPxSolver::OPTIMAL)
   {
      mysoplex.getPrimalRational(prim);
      mysoplex.getDualRational(dual);
      std::cout << "LP solved to optimality.\n";
      std::cout << "Objective value is " << mysoplex.objValueRational() << ".\n";
      cout << "Primal solution is [ "; 
      for(int j = 0 ;  j<=nbA; j++){
        cout << prim[j] << ", ";
        primDouble.push_back(prim[j]);
      }
      cout << " ]" << endl;
      std::cout << "Error solution is [" << prim[nbA+1] << "].\n";
      primDouble.push_back(prim[nbA+1]);
      return primDouble;
   }
   else
   {
      std::cout << "Error: SoPlex returned with status " << stat << ".\n";
      return primDouble;
   }
   
}

//step 1 using glpk
void step1(std::vector<double> D, int sizeD, int nbA){
  createMPSFile(D, sizeD, nbA);
  glp_prob *mip;
  glp_tran *tran;
  int ret;
  mip = glp_create_prob();
  tran = glp_mpl_alloc_wksp();
  ret = glp_mpl_read_model(tran, "models/step1.mod", 1);
  ret = glp_mpl_generate(tran, NULL);
  glp_mpl_build_prob(tran, mip);
  glp_simplex(mip, NULL);
  glp_intopt(mip, NULL);
  glp_print_sol(mip, "../resultsOfSteps/stepResults/step1Result.txt");
  glp_mpl_free_wksp(tran);
  glp_delete_prob(mip);   
}

//Creates the model to solve
void createMPSFile(std::vector<double> D, int sizeD, int nbA){
  ofstream step1("models/step1.mod");
  double x = 0.0;
  step1 << "# Parameters \n";
  step1 << "\n";
  
  step1 << "# Variable \n";
  for(int i = 0; i<=nbA; i++){
    step1 << "var a" << i << ";\n";
  }
  step1 << "var y; \n";
  
  step1 << "# Objectif \n";  
  step1 << "minimize z : y; \n";
  step1 << "\n";
  
  step1 << "# Constraint \n";  
  for(int i = 0; i<sizeD; i++){
    step1 << "d" << i << " : a0 ";
    for(int j = 1; j<=nbA; j++){
      x = pow(D[i], j);
      step1 << " + (" << x << "*a" << j << ")" ;
    }
    step1 << "-(" << setprecision(13)  << f(D[i]) << ")<= y; \n" ;
  }
    
  for(int i = 0; i<sizeD; i++){
    step1 << "d" << i+sizeD << " : a0 ";
    for(int j = 1; j<=nbA; j++){
      x = pow(D[i], j);
      step1 << " + (" << x << "*a" << j << ")" ;
    }
    step1 << "-(" << setprecision(13)  << f(D[i]) << ")>= -y; \n";
  }
  
    step1 << "\n solve ; \n";  
    for(int j = 0; j<=nbA; j++){
      step1 << "display a" << j << "; \n";
    }
    step1.close();
}


//Parser recuperating the string of the polynomial 
//(no need to change them in double since it is to write them in step 2)
std::vector<string> parserResultStep1(int nbA){
  ifstream resultStep1("../resultsOfSteps/stepResults/step1Result.txt");
  std::vector<string> P = {};
  if(resultStep1) {
    bool varReached = false;
    string str;
    resultStep1 >> str;
    //To reach the first a
    while(str != "a0"){
      resultStep1 >> str;
    }
    //Catching the value of each ai
    for(int i = 0; i <= nbA+1; i++){
      resultStep1 >> str;
      //Depending of the state of the variable, the number of necessary jump 
      //to reach next variable changes.
      if (str == "NF"){ 
        // If the variable isn't used, its value is 0
        resultStep1 >> str;
        P.push_back("0");
        resultStep1 >> str;
        resultStep1 >> str;
        resultStep1 >> str;
        resultStep1 >> str;

      } else {      
        resultStep1 >> str; 
        P.push_back(str);
        resultStep1 >> str;
        resultStep1 >> str;
      }
    }
  }
  else {
    cout << "ERROR FILE NOT OPENING" << endl;
  }
  return P;
}

//STEP 2 :

void createModelStep2Horner(std::vector<string> P, int nbA, string function, double up, double down){
  ofstream step2("models/step2_model.txt");
  step2 << "Variables\n";
  step2 << "  x in [" << down <<" , " << up <<"];\n";
  step2 << "\nMinimize\n";
  step2 << "  -abs(" ;
  step2 << "(" << P[0] << ")";
  for(int i = 1; i<nbA+1; i++){
    step2 << "+x*((" << P[i] << ")";
  }
  for(int i = 1; i<nbA+1; i++){
    step2 << ")";
  }
  step2 << "-(" << function << ")); \n";

  step2 << "\nConstraints\n";
  step2 << "end";
  step2.close();
}

void createModelStep2(std::vector<string> P, int nbA, string function, double up, double down){
  ofstream step2("models/step2_model.txt");
  step2 << "Variables\n";
  step2 << "  x in [" << down <<" , " << up <<"];\n";
  step2 << "\nMinimize\n";
  step2 << "  -abs(" ;
  step2 << "(" << P[0] << ") +";
  for(int i = 1; i<nbA+1; i++){
    step2 << "((x^" << i << ")*(" << P[i] << ")) +";
  }
  step2 << "(-(" << function << "))); \n";

  step2 << "\nConstraints\n";
  step2 << "end";
  step2.close();
}

void createModelStep2Sop(std::vector<double> P, int nbA, string function, double up, double down){
  ofstream step2("models/step2_model.txt");
  step2 << "Variables\n";
  step2 << "  x in [" << down <<" , " << up <<"];\n";
  step2 << "\nMinimize\n";
  step2 << "  -abs(" ;
  step2 << "(" << setprecision(13) << P[0] << ") +";
  for(int i = 1; i<nbA+1; i++){
    step2 << "((x^" << i << ")*(" << setprecision(13) << P[i] << ")) +";
  }
  step2 << "(-(" << function << "))); \n";

  step2 << "\nConstraints\n";
  step2 << "end";
  step2.close();
}

void step2(std::vector<string> P, int nbA, string function, double up, double down){
  createModelStep2(P, nbA, function, up, down);
  ofstream step2R("../resultsOfSteps/stepResults/step2_result.txt");

  System sys("models/step2_model.txt");
  DefaultOptimizer optimizer(sys,1e-09, 1e-10);
  optimizer.optimize(sys.box);
  optimizer.report();
  
  //cout << "Writing step2 solution in : ../resultsOfSteps/step2_result.txt\n";
  step2R << "minimizer: " << optimizer.get_loup_point() << endl;
  step2R << "uplo " << optimizer.get_uplo() << "\n" << endl;
  step2R << "loup " << optimizer.get_loup() << endl;
  step2R << optimizer.get_data();
  step2R.close();
}

void step2Sop(std::vector<double> P, int nbA, string function, double up, double down){
  createModelStep2Sop(P, nbA, function, up, down);
  ofstream step2R("../resultsOfSteps/stepResults/step2_result.txt");

  System sys("models/step2_model.txt");
  DefaultOptimizer optimizer(sys,1e-09, 1e-10);
  optimizer.optimize(sys.box);
  optimizer.report();
  
  //cout << "Writing step2 solution in : ../resultsOfSteps/step2_result.txt\n";
  step2R << "minimizer: " << optimizer.get_loup_point() << endl;
  step2R << "uplo " << optimizer.get_uplo() << "\n" << endl;
  step2R << "loup " << optimizer.get_loup() << endl;
  step2R << optimizer.get_data();
  step2R.close();
}

std::vector<double> parserResultStep2(){
  ifstream resultStep2("../resultsOfSteps/stepResults/step2_result.txt");  //Ouverture d'un fichier en lecture
  std::vector<double> returnValue = {};
  double result = 0.0;
  double error = 0.0;
  if(resultStep2) {
    string str;
    resultStep2 >> str;
    resultStep2 >> str;
    resultStep2 >> str;
    str.erase(std::remove(str.begin(), str.end(), '>'), str.end());
    str.erase(std::remove(str.begin(), str.end(), ')'), str.end());
    istringstream  istr(str);
    istr >> result;
    
    resultStep2 >> str;
    resultStep2 >> str;
    resultStep2 >> str; 

    resultStep2 >> str;
    resultStep2 >> str;
    resultStep2 >> str;  
    resultStep2 >> str;   
    
    string stringDown, stringUp;
    double bigger, smaller;
    
    str.erase(std::remove(str.begin(), str.end(), ']'), str.end());
    str.erase(std::remove(str.begin(), str.end(), '['), str.end());
    stringDown = str.substr(0, str.find(","));
    stringUp = str.substr(str.find(",")+1, str.length());
    istringstream  istr2(stringDown);
    istringstream  istr3(stringUp);
    istr2 >> smaller;
    istr3 >> bigger;
    error = (bigger + smaller)/2;
  }
  else {
    cout << "ERROR FILE NOT OPENING" << endl;
  }
  returnValue.push_back(error);
  returnValue.push_back(result);
  return returnValue;
}


/*---------------------------REMEZ 1+ UNIVARIE--------------------------------------*/


int RemezUniPlusSoplex(int nbA, double approx, int sizeD, int nbTurns, double up, double down, string function, string derivedFunction){
  ofstream summary("../resultsOfSteps/summary_of_executions/Summary_Remez_1_Plus_Soplex");
  summary << "Welcome to the summary of Remez 1+ for univariate polynomial approximation. \n";
  summary << "The parameters of Remez have been initialized as : \n";
  int i = 1;
  bool temp = false;
  bool itIsTheEnd = false;
  std::vector<bool> pointsToKeep;
  double errorStep1;
  string path = "../resultsOfSteps/img/Remez1+/soplex/";
  double error2;
  std::vector<double> P;
  std::vector<double> errorStep={};
  std::vector<double> errorStep2={};
  std::vector<std::vector<double>> newPoints;
  std::vector<double> D = step0Univarie(sizeD, up, down);
  summary << "-Degree of the polynomial approximation : " << nbA + 1 << "\n";
  summary << "-Size of the first discretized set : " << sizeD << "\n";
  summary << "-Approximation of x : " << approx << "\n";
  summary << "-Function to approximate : " << function << "\n";
  summary << "-Maximum number of turn allowed : " << nbTurns << "\n";
  summary << "-x borners : [" << up << " , " << down << "]\n";
  summary << "BEGINNING : \n" ;
  while(!itIsTheEnd){
    summary << "\nTurn " << i << "\n";
    summary << "[" << D[0];
    for(int j = 1; j<sizeD; j++){
      summary << ","  << D[j];
    }
    summary << "]\n";
    P=step1Soplex(D, sizeD, nbA);
    summary << "[" << P[0] ;
    for(int j = 1; j<=nbA; j++){
      summary << "," << P[j] ;
    }
    summary << "] \n" ;
    for(int j = 0; j<=nbA; j++){
      summary << "-a" << j << "= " << P[j] << "\n";
    }
    errorStep1 = P[nbA+1];
    errorStep.push_back(errorStep1);
    P.pop_back();
    step2PlusSop(P, nbA, derivedFunction, up, down);
    newPoints = parserResultStep2Plus();
    pointsToKeep = tooClosePlus(newPoints, D, approx);
    summary << "new points : [";
    std::vector<double> np = {};
    error2 = 0;
    for(int j=0; j<newPoints.size(); j++){
      if (pointsToKeep[j]==1){
        D.push_back(newPoints[j][0]);
	summary << newPoints[j][0] << ", ";
        sizeD = sizeD + 1;
      }
      error2 = std::max(error2, abs(Poly(newPoints[j][0], P)-f(newPoints[j][0])));
    }
    summary << "]\n";
    for(int j=0; j<pointsToKeep.size(); j++){
      temp = temp || pointsToKeep[j];
    }
    itIsTheEnd = (!temp) || (i>nbTurns);
    temp = false;
    summary << "Error for x in D" << i << ": " << setprecision(13)  << errorStep1 << "\n" ;
    std::vector<double> borners = {};
    borners.push_back(down);
    borners.push_back(up);
    makeGraphUniPlus(P, function, i, borners, error2, path);
    errorStep2.push_back(error2);
    itIsTheEnd = itIsTheEnd || (errorStep[errorStep.size()-1] >= errorStep2[errorStep.size()-1]-approx && errorStep[errorStep.size()-1] <= errorStep2[errorStep.size()-1]+approx);
    i = i + 1; 
  }
  convergenceCalculation(errorStep, errorStep2, path);
  grapheConvergenceRemez1(errorStep, errorStep2, i-1, function, path);
  summary << "\n Remez finished after " << i-1 << " turns. \n";
  return i;
}


int RemezUniPlus(int nbA, double approx, int sizeD, int nbTurns, double up, double down, string function, string derivedFunction){
  ofstream summary("../resultsOfSteps/summary_of_executions/Summary_Remez_1_Plus");
  summary << "Welcome to the summary of Remez 1+ for univariate polynomial approximation. \n";
  summary << "The parameters of Remez have been initialized as : \n";
  int i = 1;
  bool temp = false;
  bool itIsTheEnd = false;
  std::vector<bool> pointsToKeep;
  double errorStep1;
  std::vector<double> errorStep = {};
  std::vector<double> errorStep2 = {};
  std::vector<string> P;
  string path = "../resultsOfSteps/img/Remez1+/glpk/";
  double a;
  double error2;
  std::vector<double> Pdouble = {};
  std::vector<std::vector<double>> newPoints;
  std::vector<double> D = step0Univarie(sizeD, up, down);
  summary << "-Degree of the polynomial approximation : " << nbA + 1 << "\n";
  summary << "-Size of the first discretized set : " << sizeD << "\n";
  summary << "-Approximation of x : " << approx << "\n";
  summary << "-Function to approximate : " << function << "\n";
  summary << "-Maximum number of turn allowed : " << nbTurns << "\n";
  summary << "-x borners : [" << up << " , " << down << "]\n";
  summary << "BEGINNING : \n" ;
  while(!itIsTheEnd){
    summary << "\nTurn " << i << "\n";
    summary << "[" << D[0];
    for(int j = 1; j<sizeD; j++){
      summary << ","  << D[j];
    }
    summary << "]\n";
    step1(D, sizeD, nbA);
    P = parserResultStep1(nbA);
    summary << "[" << P[0] ;
    istringstream  istr(P[0]);
    istr >> a;
    Pdouble.push_back(a);
    for(int j = 1; j<=nbA; j++){
      summary << "," << P[j] ;
      istringstream  istr(P[j]);
      istr >> a;
      Pdouble.push_back(a);
    }
    summary << "] \n" ;
    for(int j = 0; j<=nbA; j++){
      summary << "-a" << j << "= " << P[j] << "\n";
    }
    istringstream  istr1(P[nbA+1]);
    istr1 >> errorStep1;
    errorStep.push_back(errorStep1);
    P.pop_back();
    step2Plus(P, nbA, derivedFunction, up, down);
    newPoints = parserResultStep2Plus();
    pointsToKeep = tooClosePlus(newPoints, D, approx);
    summary << "new points : [";
    error2 = 0;
    for(int j=0; j<newPoints.size(); j++){
      if (pointsToKeep[j]==1){
        D.push_back(newPoints[j][0]);
	summary << newPoints[j][0] << ", ";
        sizeD = sizeD + 1;
      }
      error2 = std::max(error2, abs(Poly(newPoints[j][0], Pdouble)-f(newPoints[j][0])));
    }
    errorStep2.push_back(error2);
    summary << "]\n";
    for(int j=0; j<pointsToKeep.size(); j++){
      temp = temp || pointsToKeep[j];
    }
    itIsTheEnd = (!temp) || (i>nbTurns);
    itIsTheEnd = itIsTheEnd || (errorStep[errorStep.size()-1] >= errorStep2[errorStep.size()-1]-approx && errorStep[errorStep.size()-1] <= errorStep2[errorStep.size()-1]+approx);
    temp = false;
    summary << "Error for x in D" << i << ": " << setprecision(13)  << errorStep1 << "\n" ;
    std::vector<double> borners = {};
    borners.push_back(down);
    borners.push_back(up);
    makeGraphUniPlus(Pdouble, function, i, borners, error2, path);
    Pdouble = {};
    i = i + 1; 
  }
  convergenceCalculation(errorStep, errorStep2, path);
  grapheConvergenceRemez1(errorStep, errorStep2, i-1, function, path);
  summary << "\n Remez finished after " << i-1 << " turns. \n";
  return i;
}

std::vector<bool> tooClosePlus(std::vector<std::vector<double>> newPoints, std::vector<double> D, double approx) {
  std::vector<bool> whichAreNTClose;
  bool temp = false;
  for(int i = 0; i<newPoints.size(); i++){
    for(int j = 0; j<D.size(); j++){
      if(newPoints[i][0]<=D[j]+approx && newPoints[i][0]>=D[j]-approx){
        temp = true;
      }
      if(newPoints[i][1]<=D[j]+approx && newPoints[i][1]>=D[j]-approx){
        temp = temp && true;
      }
    }
    whichAreNTClose.push_back(!temp);
    temp = false;
  }
  return whichAreNTClose;
}


std::vector<std::vector<double>> parserResultStep2Plus(){
  ifstream resultStep2("../resultsOfSteps/stepResults/step2Plus_result.txt");  //Ouverture d'un fichier en lecture
  std::vector<double> R = {};
  std::vector<std::vector<double>> S = {};
  string smaller, bigger;
  double small, big;
  double result = 0.0;
  if(resultStep2) {
    string str;
    resultStep2 >> str;
    
    while(str == "solution"){
      resultStep2 >> str;
      resultStep2 >> str;
      resultStep2 >> str;
      smaller = str;
      smaller.erase(std::remove(smaller.begin(), smaller.end(), '['), smaller.end());
      smaller.erase(std::remove(smaller.begin(), smaller.end(), '('), smaller.end());
      smaller.erase(std::remove(smaller.begin(), smaller.end(), ','), smaller.end());
      resultStep2 >> str;
      bigger = str;
      bigger.erase(std::remove(bigger.begin(), bigger.end(), ']'), bigger.end());
      bigger.erase(std::remove(bigger.begin(), bigger.end(), ')'), bigger.end());
      istringstream  istr2(bigger);
      istringstream  istr3(smaller);
      istr2 >> big;
      istr3 >> small;
      R.push_back(small);
      R.push_back(big);
      resultStep2 >> str;
      S.push_back(R);
      R = {};
    }
  }
  else {
    cout << "ERROR FILE NOT OPENING" << endl;
  }
  return S;
}

void step2PlusSop(std::vector<double> P, int nbA, string derivedFunction, double up, double down){
  ofstream step2RPlus("../resultsOfSteps/stepResults/step2Plus_result.txt");
  createModelStep2PlusSop(P, nbA, derivedFunction, up, down);
  System system("models/step2Plus_model.txt");

  DefaultSolver solver(system,1e-07, 1e-15);
  solver.solve(system.box);
  solver.report();
  step2RPlus << solver.get_data() << endl;
}

void step2Plus(std::vector<string> P, int nbA, string derivedFunction, double up, double down){
  ofstream step2RPlus("../resultsOfSteps/stepResults/step2Plus_result.txt");
  createModelStep2Plus(P, nbA, derivedFunction, up, down);
  System system("models/step2Plus_model.txt");

  DefaultSolver solver(system,1e-07, 1e-15);
  solver.solve(system.box);
  solver.report();
  step2RPlus << solver.get_data() << endl;
}

void createModelStep2PlusSop(std::vector<double> P, int nbA, string derivedFunction, double up, double down){
  ofstream step2("models/step2Plus_model.txt");
  step2 << "Variables\n";
  step2 << "  x in [" << down <<" , " << up <<"];\n";
  step2 << "\nConstraints\n";
  
  step2 << P[1] << " + ";
  for(int i = 2 ; i<= nbA; i++){
    step2 << "((" << P[i] <<"*" << i <<")*x" << "^" << i-1 << ")+";
  }
  step2 << "-(" << derivedFunction << ") = 0 ;\n";

  step2 << "end";
  step2.close();
}


void createModelStep2Plus(std::vector<string> P, int nbA, string derivedFunction, double up, double down){
  ofstream step2("models/step2Plus_model.txt");
  step2 << "Variables\n";
  step2 << "  x in [" << down <<" , " << up <<"];\n";
  step2 << "\nConstraints\n";
  
  step2 << P[1] << " + ";
  for(int i = 2 ; i<= nbA; i++){
    step2 << "((" << P[i] <<"*" << i <<")*x" << "^" << i-1 << ")+";
  }
  step2 << "-(" << derivedFunction << ") = 0 ;\n";

  step2 << "end";
  step2.close();
}


/*---------------------------REMEZ 1 MULTIVARIE------------------------------------*/

int RemezMulti(){
  ofstream summary("../resultsOfSteps/summary_of_executions/Summary_Remez_1_Multivariate");
  summary << "Welcome to the summary of Remez 1 for multivariate polynomial approximation. \n";
  summary << "The parameters of Remez have been initialized as : \n";
  int sizeD = 5;
  int i = 1;
  int maxTurn = 100;
  int nbVar = 2;
  double approx = 0.00000005;
  std::vector<double> S;
  bool itIsTheEnd = false;
  string functionMulti = "cos(x0) + sin(x1)";
  std::vector<std::vector<double>> borners = {{-1, 1},{-1, 1}};
  std::vector<std::vector<double>> D = step0_multi(sizeD, nbVar, borners);
  int nbA = 5;
 while(!itIsTheEnd){
    step1_multivariate(D, sizeD, nbVar, nbA);
    std::vector<string>P = parserResultStep1_multivariate(nbA, nbVar);
    
    summary << "\nTurn " << i << ": \n";
    summary << "Error Step1 : " << P[P.size()-1] << endl;
    P.pop_back(); 
       
    summary << "["<< P[0] << ", ";
    for(int j=1; j<P.size()-1; j++){
      summary << P[j] << ", ";
    }
    summary << P[P.size()-1] << "]\n";
    
    step2_multi(P, nbA, nbVar, functionMulti, borners);
    S = parserResultStep2Multi(nbVar);
    
    summary << "Error Step2 = " << S[S.size()-1];
    S.pop_back();
    itIsTheEnd = (i>maxTurn);
    if(!tooCloseMulti(D, nbVar, sizeD, S, approx)){
      D.push_back(S);
    } else {
      itIsTheEnd = true;
    }
    i = i + 1;
  }
  return i;
}


bool tooCloseMulti(std::vector<std::vector<double>> D, int nbVar, int sizeD, std::vector<double> newPoint, double approx){
  bool isClose;
  for(int i=0 ; i<D.size(); i++){
    isClose = true;
    for(int j = 0; j<nbVar; j++){
      if(!(D[i][j]<=newPoint[j]+approx && D[i][j]>=newPoint[j]-approx)){ isClose = false; }
    }
    if(isClose == true){
      return isClose;
    }
  }
  return isClose; 
}

std::vector<std::vector<double>> step0_multi(int sizeD, int nbVar, std::vector<std::vector<double>> borners){

  std::vector<double> steps = {};
  for (int i = 0; i<nbVar; i++){
    steps.push_back((borners[i][1] - borners[i][0])/(sizeD-1));
  }
  std::vector<std::vector<double>> S = {};
  
  for(int k = 0; k<pow(sizeD, nbVar-1); k++){
    for (int i = 0; i<sizeD; i++){
      S.push_back({borners[0][0]+steps[0]*(i)});
      
    }
  }
  
  for (int j = 1; j<nbVar; j++){
    for (int i = 0; i<pow(sizeD, nbVar); i++){
      S[i].push_back(borners[j][0]+steps[j]*((i+(j*i/sizeD))%sizeD));
    }
  }
  
  /* Debug D0
  for(int i = 0; i< S.size(); i++){
    cout << "[";
    for(int j=0; j<S[i].size(); j++){
      cout << S[i][j] << ", ";
    }
    cout << "] \n";
  }*/
  return S;
}


std::vector<double> parserResultStep2Multi(int nbVar){
  ifstream resultStep2("../resultsOfSteps/stepResults/step2_mult_result.txt");  //Ouverture d'un fichier en lecture
  std::vector<string> P = {};
  std::vector<double> result = {};
  double temp;
  //double error = 0.0;
  if(resultStep2) {
  
    string str;
    resultStep2 >> str;
    resultStep2 >> str;
    resultStep2 >> str;
    for(int i = 0; i<nbVar; i++){
      istringstream  istr(str);
      istr >> temp;
      result.push_back(temp);
      resultStep2 >> str;
      resultStep2 >> str;
      resultStep2 >> str; 
    }
    resultStep2 >> str;
    istringstream  istr(str);
    istr >> temp;
    result.push_back(temp);
  }
  else {
    cout << "ERROR FILE NOT OPENING" << endl;
  }
  return result;
}

void step2_multi(std::vector<string> P, int nbA, int nbVar, string function, std::vector<std::vector<double>> borners){
  createModelStep2_multivariate(P, nbA, nbVar, function, borners);
  ofstream step2R("../resultsOfSteps/stepResults/step2_mult_result.txt");

  System sys("models/step2_multi_model.txt");
  DefaultOptimizer optimizer(sys,1e-10);
  optimizer.optimize(sys.box);
  optimizer.report();
  
  //step2R <<  "f* in " << Interval(optimizer.get_uplo(),optimizer.get_loup()) << endl;
  //cout << "Writing step2 solution in : ../resultsOfSteps/step2_result.txt\n";
  step2R << "minimizer: " << optimizer.get_loup_point() << endl;
  step2R << "uplo " << optimizer.get_uplo() << "\n" << endl;
  step2R << "loup " << optimizer.get_loup() << endl;
  step2R << optimizer.get_data();
  step2R.close();

}


void createModelStep2_multivariate(std::vector<string> P, int nbA, int nbVar, string function, std::vector<std::vector<double>> borners){
  ofstream step2("models/step2_multi_model.txt");
  step2 << "Variables\n";
    std::vector<int> index = {};
  for(int i=0; i<nbVar; i++){
    step2 << "  x" << i <<" in [" << borners[i][0] <<" , " << borners[i][1] <<"];\n";
  }
  
  step2 << "\nMinimize\n";
  step2 << "  -abs(" ;
  
  for(int i = 0; i<P.size(); i++){
    step2 << "(" << P[i] << ")";
    index = getIndice(nbA, i);
    for(int k = index.size(); k < nbVar; k++){
      index.push_back(0);
    }
    for(int j = 0; j<nbVar; j++){
      if (index[j]>0){
        step2 << "*(x" << j << "^" << index[j] << ")";
      }
    }
    step2 << "+" ;
  }
  
  step2 << "(-" << function << ")); \n";
  step2 << "\nConstraints\n";
  step2 << "end";
  step2.close();
}

std::vector<string> parserResultStep1_multivariate(int nbA, int nbVar){
  ifstream resultStep1("../resultsOfSteps/stepResults/step1_multivariate_Result.txt");  //Ouverture d'un fichier en lecture
  std::vector<string> P = {};
  string firstVar = "a";
  for(int i = 0; i<nbVar; i++){
    firstVar += "0";
  }
  if(resultStep1) {
    bool varReached = false;
    string str;
    resultStep1 >> str;
    while(str != firstVar){
      resultStep1 >> str;
    }
    for(int i = 0; i <= pow(nbA+1,nbVar); i++){
      resultStep1 >> str;
      
      if (str == "NF"){
        resultStep1 >> str;
        P.push_back("0");
        resultStep1 >> str;
        resultStep1 >> str;
        resultStep1 >> str;
        resultStep1 >> str;

      } else {
      
        resultStep1 >> str;
      
        P.push_back(str);
        resultStep1 >> str;
        resultStep1 >> str;
      }
    }
  }
  else {
    cout << "ERROR FILE NOT OPENING" << endl;
  }
  return P;
}


void step1_multivariate(std::vector<std::vector<double>> D, int sizeD, int nbVar, int nbA){
  createMPSFileMultivariate(D, sizeD, nbVar ,nbA);
  glp_prob *mip;
  glp_tran *tran;
  int ret;
  mip = glp_create_prob();
  tran = glp_mpl_alloc_wksp();
  ret = glp_mpl_read_model(tran, "models/step1_multi.mod", 1);
  ret = glp_mpl_generate(tran, NULL);
  glp_mpl_build_prob(tran, mip);
  glp_simplex(mip, NULL);
  glp_intopt(mip, NULL);
  glp_print_sol(mip, "../resultsOfSteps/stepResults/step1_multivariate_Result.txt");
  glp_mpl_free_wksp(tran);
  glp_delete_prob(mip);  
}

void createMPSFileMultivariate(std::vector<std::vector<double>> D, int sizeD, int nbVar, int nbA){
  ofstream step1("models/step1_multi.mod");
  double x;
  step1 << "# Parameters \n";
  step1 << "\n";
  std::vector<int> index;
  step1 << "# Variable";
  
  for(int m = 0; m<pow(nbA+1,nbVar); m++){
    index = getIndice(nbA, m);
    for(int k = index.size(); k < nbVar; k++){
      index.push_back(0);
    }
    step1 << "\n var a";
    for(int j = 0; j<nbVar; j++){
      step1 << index[j] ;
    } 
    step1 << ";";
  }
  step1 << "\nvar y; \n";
  
  step1 << "# Objectif \n";  
  step1 << "minimize z : y; \n";
  step1 << "\n";
  
  step1 << "# Constraint \n";  
  
  for(int i = 0; i<D.size(); i++){
    step1 << "d" << i << " : ";
    for(int m = 0; m<pow(nbA+1,nbVar); m++){
      index = getIndice(nbA, m);
      
      for(int k = index.size(); k < nbVar; k++){
        index.push_back(0);
      }
      
      step1 << "a";
      x = 1;
      for(int j = 0; j<nbVar; j++){
        step1 << index[j] ;
        x *= pow(D[i][j], index[j]);
      } 
      step1 << "*("<<x<<")+"; 
    }
    step1 << "(" << -fMulti(D[i]) << ") <= y; \n" ;
  }
  
    for(int i = 0; i<D.size(); i++){
    step1 << "d" << i+D.size() << " : ";
    for(int m = 0; m<pow(nbA+1,nbVar); m++){
      index = getIndice(nbA, m);
      
      for(int k = index.size(); k < nbVar; k++){
        index.push_back(0);
      }
      
      step1 << "a" ;
      x = 1;
      for(int j = 0; j<nbVar; j++){
        step1 << index[j] ;
        x *= pow(D[i][j], index[j]); 
      } 
      step1 << "*("<<x<<")+"; 
    }
    step1 << "(" << -fMulti(D[i]) << ") >= -y; \n" ;
  }

  step1 << "\n solve ; \n";  
    step1.close();
}

std::vector<int> getIndice(int nbA, int whereWeAre){
  std::vector<int> I = {};
  std::vector<int> L = {};
  if (whereWeAre <= nbA) {
    I.push_back(whereWeAre);
    return I;
  } else {
    int i = whereWeAre/(nbA+1);
    int reste = whereWeAre % (nbA+1);
    L = getIndice(nbA, i);
    I.push_back(reste);
    I.insert(I.end(), L.begin(), L.end());
    return I;
  }
  return I;
}



/*---------------------------REMEZ 1+ MULTIVARIE------------------------------------*/























//--------------------Trash-----------------//


/* Create a ibex file for step 2 : unused, untested
void createIbexFile(std::vector<double> P, int nbA, string function, double up, double down){
  ofstream step2("../step2/step2.cpp");
  step2 << "#include <ibex.h>" << "\n";
  step2 << "using namespace ibex;" << "\n";
  step2 << "using namespace std;" << "\n" << "\n";
  step2 << "int main(int argc, char** argv) { " << "\n";
  step2 << "  ofstream step2R(\"../resultsOfSteps/step2_result.txt\");";
  step2 << "  Variable x;" << "\n";
  step2 << "  SystemFactory fac;" << "\n";
  step2 << "  fac.add_var(x);" << "\n";
  step2 << "fac.add_goal(-abs(";
  for(int i = 0; i<nbA; i++){
    step2 << "(pow(x," << i << ")*" << P[i] << ") +";
  }
  step2 << "(-" << function << "))); \n";
  step2 << "  fac.add_ctr(x >=" << down << ");" << "\n";
  step2 << "  fac.add_ctr(x <=" << up << ");" << "\n";
  step2 << "  System sys(fac);" << "\n";
  step2 << "  DefaultOptimizer optimizer(sys,1e-01);";
  step2 << "  step2R <<  optimizer.get_loup_point() << endl;";
  step2 << "}";
  step2.close();
} */


/*Version 1 : glpk direct. A un probleme dans la matrix
int main(int argc, char** argv) { 
  bool debug = true;
  if(debug){
    cout << "Welcome to Remez_1 DEBUG_MODE !"<< "\n";
    cout << "The very first step will be making the first Discretized domain D" << "\n";
  }
  
  std::vector<double> D = {-1., -0.5, 0., 0.5, 1.} ;
  
  glp_prob *lp;
  int ia[1+1000], ja[1+1000];
  double ar[1+1000], z, y, a0, a1, a2, f;
  lp = glp_create_prob();
  glp_set_prob_name(lp, "step 1");
  glp_set_obj_dir(lp, GLP_MIN);
  
  glp_add_rows(lp, 10);
  glp_set_row_name(lp, 1, "d1+");
  glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 0.0);
  glp_set_row_name(lp, 2, "d2+");
  glp_set_row_bnds(lp, 2, GLP_UP, 0.0, 0.0);
  glp_set_row_name(lp, 3, "d3+");
  glp_set_row_bnds(lp, 3, GLP_UP, 0.0, 0.0);
  glp_set_row_name(lp, 4, "d4+");
  glp_set_row_bnds(lp, 4, GLP_UP, 0.0, 0.0);
  glp_set_row_name(lp, 5, "d5+");
  glp_set_row_bnds(lp, 5, GLP_UP, 0.0, 0.0);
  
  glp_set_row_name(lp, 6, "d1-");
  glp_set_row_bnds(lp, 6, GLP_LO, 0.0, 0.0);
  glp_set_row_name(lp, 7, "d2-");
  glp_set_row_bnds(lp, 7, GLP_LO, 0.0, 0.0);
  glp_set_row_name(lp, 8, "d3-");
  glp_set_row_bnds(lp, 8, GLP_LO, 0.0, 0.0);
  glp_set_row_name(lp, 9, "d4-");
  glp_set_row_bnds(lp, 9, GLP_LO, 0.0, 0.0);
  glp_set_row_name(lp, 10, "d5-");
  glp_set_row_bnds(lp, 10, GLP_LO, 0.0, 0.0);
  
  glp_add_cols(lp, 5);
  glp_set_col_name(lp, 1, "a0");
  glp_set_col_bnds(lp, 1, GLP_FR, 0.0, 0.0);
  glp_set_obj_coef(lp, 1, 0.0);
  glp_set_col_name(lp, 2, "a1");
  glp_set_col_bnds(lp, 2, GLP_FR, 0.0, 0.0);
  glp_set_obj_coef(lp, 2, 0.0);
  glp_set_col_name(lp, 3, "a2");
  glp_set_col_bnds(lp, 3, GLP_FR, 0.0, 0.0);
  glp_set_obj_coef(lp, 3, 0.0);
  glp_set_col_name(lp, 4, "y");
  glp_set_col_bnds(lp, 4, GLP_FR, 0.0, 0.0);
  glp_set_obj_coef(lp, 4, 1.0);
  glp_set_col_name(lp, 5, "f");
  glp_set_col_bnds(lp, 5, GLP_FR, 0.0, 0.0);
  glp_set_obj_coef(lp, 5, 0.0);
  
  ia[1] = 1, ja[1] = 1, ar[1] = -1.0; //a0
  ia[2] = 1, ja[2] = 2, ar[2] = 1.0;  //a1
  ia[3] = 1, ja[3] = 3, ar[3] = -1.0;  //a2
  ia[4] = 1, ja[4] = 4, ar[4] = 1.0;  //y
  ia[5] = 1, ja[5] = 5, ar[5] = cos(-1.0); //f
  ia[6] = 2, ja[6] = 1, ar[6] = -1.0; //a0
  ia[7] = 2, ja[7] = 2, ar[7] = 0.5;
  ia[8] = 2, ja[8] = 3, ar[8] = -0.25;
  ia[9] = 2, ja[9] = 4, ar[9] = 1.0;
  ia[10] = 2, ja[10] = 5, ar[10] = cos(-0.5); //f
  ia[11] = 3, ja[11] = 1, ar[11] = -1.0; //a0
  ia[12] = 3, ja[12] = 2, ar[12] = 0.0;
  
  ia[13] = 3, ja[13] = 3, ar[13] = 0.0;
  ia[14] = 3, ja[14] = 4, ar[14] = 1.0;
  ia[15] = 3, ja[15] = 5, ar[15] = cos(0.0); //f
  ia[16] = 4, ja[16] = 1, ar[16] = -1.0; //a0
  ia[17] = 4, ja[17] = 2, ar[17] = -0.5;
  ia[18] = 4, ja[18] = 3, ar[18] = -0.25;
  ia[19] = 4, ja[19] = 4, ar[19] = 1.0;
  ia[20] = 4, ja[20] = 5, ar[20] = cos(0.5); //f
  ia[21] = 5, ja[21] = 1, ar[21] = -1.0; //a0
  ia[22] = 5, ja[22] = 2, ar[22] = -1.0;
  ia[23] = 5, ja[23] = 3, ar[23] = -1.0;
  ia[24] = 5, ja[24] = 4, ar[24] = 1.0;
  ia[25] = 5, ja[25] = 5, ar[25] = cos(1.0); //f
  
  ia[26] = 6, ja[26] = 1, ar[26] = 1.0; //a0
  ia[27] = 6, ja[27] = 2, ar[27] = -1.0;
  ia[28] = 6, ja[28] = 3, ar[28] = 1.0;
  ia[29] = 6, ja[29] = 4, ar[29] = 1.0;
  ia[30] = 6, ja[30] = 5, ar[30] = -cos(-1.0); //f
  ia[31] = 7, ja[31] = 1, ar[31] = 1.0; //a0
  ia[32] = 7, ja[32] = 2, ar[32] = -0.5;
  ia[33] = 7, ja[33] = 3, ar[33] = 0.25;
  ia[34] = 7, ja[34] = 4, ar[34] = 1.0;
  ia[35] = 7, ja[35] = 5, ar[35] = -cos(-0.5); //f
  ia[36] = 8, ja[36] = 1, ar[36] = 1.0; //a0
  ia[37] = 8, ja[37] = 2, ar[37] = 0.0;
  
  ia[38] = 8, ja[38] = 3, ar[38] = 0.0;
  ia[39] = 8, ja[39] = 4, ar[39] = 1.0;
  ia[40] = 8, ja[40] = 5, ar[40] = -cos(0.0); //f
  ia[41] = 9, ja[41] = 1, ar[41] = 1.0; //a0
  ia[42] = 9, ja[42] = 2, ar[42] = 0.5;
  ia[43] = 9, ja[43] = 3, ar[43] = 0.25;
  ia[44] = 9, ja[44] = 4, ar[44] = 1.0;
  ia[45] = 9, ja[45] = 5, ar[45] = -cos(0.5); //f
  ia[46] = 10, ja[46] = 1, ar[46] = 1.0; //a0
  ia[47] = 10, ja[47] = 2, ar[47] = 1.0;
  ia[48] = 10, ja[48] = 3, ar[48] = 1.0;
  ia[49] = 10, ja[49] = 4, ar[49] = 1.0;
  ia[50] = 10, ja[50] = 5, ar[50] = -cos(1.0); //f
  
  glp_load_matrix(lp, 50, ia, ja, ar);
  glp_simplex(lp, NULL);
  z = glp_get_obj_val(lp);
  a0 = glp_get_col_prim(lp, 1);
  a1 = glp_get_col_prim(lp, 2);
  a2 = glp_get_col_prim(lp, 3);
  y = glp_get_col_prim(lp, 4);
  f = glp_get_col_prim(lp, 5);
  
  cout << z << " " << a0 << " " << a1 << " " << a2 << " " << y << "\n";
  glp_delete_prob(lp);
  return 0;
}*/

/*
//Creates the model to solve
void createMPSFileHorner(std::vector<double> D, int sizeD, int nbA){
  ofstream step1("models/step1.mod");
  double x = 0.0;
  step1 << "# Parameters \n";
  step1 << "\n";
  
  step1 << "# Variable \n";
  for(int i = 0; i<=nbA; i++){
    step1 << "var a" << i << ";\n";
  }
  step1 << "var y; \n";
  
  step1 << "# Objectif \n";  
  step1 << "minimize z : y; \n";
  step1 << "\n";
  
  step1 << "# Constraint \n";  
  for(int i = 0; i<sizeD; i++){
    step1 << "d" << i << " : a0 ";
    
   for(int j = 1; j<=nbA; j++){
      step1 << "+(" << D[i] << ")*(a" << j ;
    }
    for(int j = 1; j<=nbA; j++){
      step1 << ")" ;
    }
    step1 << "-(" << f(D[i]) << ") <= y; \n" ;
  }
    
 for(int i = 0; i<sizeD; i++){
    step1 << "d" << i+sizeD << " : a0 ";
    
   for(int j = 1; j<=nbA; j++){
      step1 << "+(" << D[i] << ")*(a" << j ;
    }
    for(int j = 1; j<=nbA; j++){
      step1 << ")" ;
    }
    step1 << "-(" << f(D[i]) << ") >= -y; \n" ;
  }

  
  step1 << "\n solve ; \n";  
  for(int j = 0; j<=nbA; j++){
    step1 << "display a" << j << "; \n";
  }
  step1.close();
}
*/


  /*
  double X[1000] = {-1. , -0.997998  , -0.995996  , -0.99399399, -0.99199199, -0.98998999,
 -0.98798799 ,-0.98598599 ,-0.98398398 ,-0.98198198 ,-0.97997998 ,-0.97797798,
 -0.97597598 ,-0.97397397 ,-0.97197197 ,-0.96996997 ,-0.96796797 ,-0.96596597,
 -0.96396396 ,-0.96196196 ,-0.95995996 ,-0.95795796 ,-0.95595596 ,-0.95395395,
 -0.95195195 ,-0.94994995 ,-0.94794795 ,-0.94594595 ,-0.94394394 ,-0.94194194,
 -0.93993994 ,-0.93793794 ,-0.93593594 ,-0.93393393 ,-0.93193193 ,-0.92992993,
 -0.92792793 ,-0.92592593 ,-0.92392392 ,-0.92192192 ,-0.91991992 ,-0.91791792,
 -0.91591592 ,-0.91391391 ,-0.91191191 ,-0.90990991 ,-0.90790791 ,-0.90590591,
 -0.9039039  ,-0.9019019  ,-0.8998999  ,-0.8978979  ,-0.8958959  ,-0.89389389,
 -0.89189189 ,-0.88988989 ,-0.88788789 ,-0.88588589 ,-0.88388388 ,-0.88188188,
 -0.87987988 ,-0.87787788 ,-0.87587588 ,-0.87387387 ,-0.87187187 ,-0.86986987,
 -0.86786787 ,-0.86586587 ,-0.86386386 ,-0.86186186 ,-0.85985986 ,-0.85785786,
 -0.85585586 ,-0.85385385 ,-0.85185185 ,-0.84984985 ,-0.84784785 ,-0.84584585,
 -0.84384384 ,-0.84184184 ,-0.83983984 ,-0.83783784 ,-0.83583584 ,-0.83383383,
 -0.83183183 ,-0.82982983 ,-0.82782783 ,-0.82582583 ,-0.82382382 ,-0.82182182,
 -0.81981982 ,-0.81781782 ,-0.81581582 ,-0.81381381 ,-0.81181181 ,-0.80980981,
 -0.80780781 ,-0.80580581 ,-0.8038038  ,-0.8018018  ,-0.7997998  ,-0.7977978,
 -0.7957958  ,-0.79379379 ,-0.79179179 ,-0.78978979 ,-0.78778779 ,-0.78578579,
 -0.78378378 ,-0.78178178 ,-0.77977978 ,-0.77777778 ,-0.77577578 ,-0.77377377,
 -0.77177177 ,-0.76976977 ,-0.76776777 ,-0.76576577 ,-0.76376376 ,-0.76176176,
 -0.75975976 ,-0.75775776 ,-0.75575576 ,-0.75375375 ,-0.75175175 ,-0.74974975,
 -0.74774775 ,-0.74574575 ,-0.74374374 ,-0.74174174 ,-0.73973974 ,-0.73773774,
 -0.73573574 ,-0.73373373 ,-0.73173173 ,-0.72972973 ,-0.72772773 ,-0.72572573,
 -0.72372372 ,-0.72172172 ,-0.71971972 ,-0.71771772 ,-0.71571572 ,-0.71371371,
 -0.71171171 ,-0.70970971 ,-0.70770771 ,-0.70570571 ,-0.7037037  ,-0.7017017,
 -0.6996997  ,-0.6976977  ,-0.6956957  ,-0.69369369 ,-0.69169169 ,-0.68968969,
 -0.68768769 ,-0.68568569 ,-0.68368368 ,-0.68168168 ,-0.67967968 ,-0.67767768,
 -0.67567568 ,-0.67367367 ,-0.67167167 ,-0.66966967 ,-0.66766767 ,-0.66566567,
 -0.66366366 ,-0.66166166 ,-0.65965966 ,-0.65765766 ,-0.65565566 ,-0.65365365,
 -0.65165165 ,-0.64964965 ,-0.64764765 ,-0.64564565 ,-0.64364364 ,-0.64164164,
 -0.63963964 ,-0.63763764 ,-0.63563564 ,-0.63363363 ,-0.63163163 ,-0.62962963,
 -0.62762763 ,-0.62562563 ,-0.62362362 ,-0.62162162 ,-0.61961962 ,-0.61761762,
 -0.61561562 ,-0.61361361 ,-0.61161161 ,-0.60960961 ,-0.60760761 ,-0.60560561,
 -0.6036036  ,-0.6016016  ,-0.5995996  ,-0.5975976  ,-0.5955956  ,-0.59359359,
 -0.59159159 ,-0.58958959 ,-0.58758759 ,-0.58558559 ,-0.58358358 ,-0.58158158,
 -0.57957958 ,-0.57757758 ,-0.57557558 ,-0.57357357 ,-0.57157157 ,-0.56956957,
 -0.56756757 ,-0.56556557 ,-0.56356356 ,-0.56156156 ,-0.55955956 ,-0.55755756,
 -0.55555556 ,-0.55355355 ,-0.55155155 ,-0.54954955 ,-0.54754755 ,-0.54554555,
 -0.54354354 ,-0.54154154 ,-0.53953954 ,-0.53753754 ,-0.53553554 ,-0.53353353,
 -0.53153153 ,-0.52952953 ,-0.52752753 ,-0.52552553 ,-0.52352352 ,-0.52152152,
 -0.51951952 ,-0.51751752 ,-0.51551552 ,-0.51351351 ,-0.51151151 ,-0.50950951,
 -0.50750751 ,-0.50550551 ,-0.5035035  ,-0.5015015  ,-0.4994995  ,-0.4974975,
 -0.4954955  ,-0.49349349 ,-0.49149149 ,-0.48948949 ,-0.48748749 ,-0.48548549,
 -0.48348348 ,-0.48148148 ,-0.47947948 ,-0.47747748 ,-0.47547548 ,-0.47347347,
 -0.47147147 ,-0.46946947 ,-0.46746747 ,-0.46546547 ,-0.46346346 ,-0.46146146,
 -0.45945946 ,-0.45745746 ,-0.45545546 ,-0.45345345 ,-0.45145145 ,-0.44944945,
 -0.44744745 ,-0.44544545 ,-0.44344344 ,-0.44144144 ,-0.43943944 ,-0.43743744,
 -0.43543544 ,-0.43343343 ,-0.43143143 ,-0.42942943 ,-0.42742743 ,-0.42542543,
 -0.42342342 ,-0.42142142 ,-0.41941942 ,-0.41741742 ,-0.41541542 ,-0.41341341,
 -0.41141141 ,-0.40940941 ,-0.40740741 ,-0.40540541 ,-0.4034034  ,-0.4014014,
 -0.3993994  ,-0.3973974  ,-0.3953954  ,-0.39339339 ,-0.39139139 ,-0.38938939,
 -0.38738739 ,-0.38538539 ,-0.38338338 ,-0.38138138 ,-0.37937938 ,-0.37737738,
 -0.37537538 ,-0.37337337 ,-0.37137137 ,-0.36936937 ,-0.36736737 ,-0.36536537,
 -0.36336336 ,-0.36136136 ,-0.35935936 ,-0.35735736 ,-0.35535536 ,-0.35335335,
 -0.35135135 ,-0.34934935 ,-0.34734735 ,-0.34534535 ,-0.34334334 ,-0.34134134,
 -0.33933934 ,-0.33733734 ,-0.33533534 ,-0.33333333 ,-0.33133133 ,-0.32932933,
 -0.32732733 ,-0.32532533 ,-0.32332332 ,-0.32132132 ,-0.31931932 ,-0.31731732,
 -0.31531532 ,-0.31331331 ,-0.31131131 ,-0.30930931 ,-0.30730731 ,-0.30530531,
 -0.3033033  ,-0.3013013  ,-0.2992993  ,-0.2972973  ,-0.2952953  ,-0.29329329,
 -0.29129129 ,-0.28928929 ,-0.28728729 ,-0.28528529 ,-0.28328328 ,-0.28128128,
 -0.27927928 ,-0.27727728 ,-0.27527528 ,-0.27327327 ,-0.27127127 ,-0.26926927,
 -0.26726727 ,-0.26526527 ,-0.26326326 ,-0.26126126 ,-0.25925926 ,-0.25725726,
 -0.25525526 ,-0.25325325 ,-0.25125125 ,-0.24924925 ,-0.24724725 ,-0.24524525,
 -0.24324324 ,-0.24124124 ,-0.23923924 ,-0.23723724 ,-0.23523524 ,-0.23323323,
 -0.23123123 ,-0.22922923 ,-0.22722723 ,-0.22522523 ,-0.22322322 ,-0.22122122,
 -0.21921922 ,-0.21721722 ,-0.21521522 ,-0.21321321 ,-0.21121121 ,-0.20920921,
 -0.20720721 ,-0.20520521 ,-0.2032032  ,-0.2012012  ,-0.1991992  ,-0.1971972,
 -0.1951952  ,-0.19319319 ,-0.19119119 ,-0.18918919 ,-0.18718719 ,-0.18518519,
 -0.18318318 ,-0.18118118 ,-0.17917918 ,-0.17717718 ,-0.17517518 ,-0.17317317,
 -0.17117117 ,-0.16916917 ,-0.16716717 ,-0.16516517 ,-0.16316316 ,-0.16116116,
 -0.15915916 ,-0.15715716 ,-0.15515516 ,-0.15315315 ,-0.15115115 ,-0.14914915,
 -0.14714715 ,-0.14514515 ,-0.14314314 ,-0.14114114 ,-0.13913914 ,-0.13713714,
 -0.13513514 ,-0.13313313 ,-0.13113113 ,-0.12912913 ,-0.12712713 ,-0.12512513,
 -0.12312312 ,-0.12112112 ,-0.11911912 ,-0.11711712 ,-0.11511512 ,-0.11311311,
 -0.11111111 ,-0.10910911 ,-0.10710711 ,-0.10510511 ,-0.1031031  ,-0.1011011,
 -0.0990991  ,-0.0970971  ,-0.0950951  ,-0.09309309 ,-0.09109109 ,-0.08908909,
 -0.08708709 ,-0.08508509 ,-0.08308308 ,-0.08108108 ,-0.07907908 ,-0.07707708,
 -0.07507508 ,-0.07307307 ,-0.07107107 ,-0.06906907 ,-0.06706707 ,-0.06506507,
 -0.06306306 ,-0.06106106 ,-0.05905906 ,-0.05705706 ,-0.05505506 ,-0.05305305,
 -0.05105105 ,-0.04904905 ,-0.04704705 ,-0.04504505 ,-0.04304304 ,-0.04104104,
 -0.03903904 ,-0.03703704 ,-0.03503504 ,-0.03303303 ,-0.03103103 ,-0.02902903,
 -0.02702703 ,-0.02502503 ,-0.02302302 ,-0.02102102 ,-0.01901902 ,-0.01701702,
 -0.01501502 ,-0.01301301 ,-0.01101101 ,-0.00900901 ,-0.00700701 ,-0.00500501,
 -0.003003   ,-0.001001  ,  0.001001   , 0.003003   , 0.00500501 , 0.00700701,
  0.00900901 , 0.01101101 , 0.01301301 , 0.01501502 , 0.01701702 , 0.01901902,
  0.02102102 , 0.02302302 , 0.02502503 , 0.02702703 , 0.02902903 , 0.03103103,
  0.03303303 , 0.03503504 , 0.03703704 , 0.03903904 , 0.04104104 , 0.04304304,
  0.04504505 , 0.04704705 , 0.04904905 , 0.05105105 , 0.05305305 , 0.05505506,
  0.05705706 , 0.05905906 , 0.06106106 , 0.06306306 , 0.06506507 , 0.06706707,
  0.06906907 , 0.07107107 , 0.07307307 , 0.07507508 , 0.07707708 , 0.07907908,
  0.08108108 , 0.08308308 , 0.08508509 , 0.08708709 , 0.08908909 , 0.09109109,
  0.09309309 , 0.0950951  , 0.0970971  , 0.0990991  , 0.1011011  , 0.1031031,
  0.10510511 , 0.10710711 , 0.10910911 , 0.11111111 , 0.11311311 , 0.11511512,
  0.11711712 , 0.11911912 , 0.12112112 , 0.12312312 , 0.12512513 , 0.12712713,
  0.12912913 , 0.13113113 , 0.13313313 , 0.13513514 , 0.13713714 , 0.13913914,
  0.14114114 , 0.14314314 , 0.14514515 , 0.14714715 , 0.14914915 , 0.15115115,
  0.15315315 , 0.15515516 , 0.15715716 , 0.15915916 , 0.16116116 , 0.16316316,
  0.16516517 , 0.16716717 , 0.16916917 , 0.17117117 , 0.17317317 , 0.17517518,
  0.17717718 , 0.17917918 , 0.18118118 , 0.18318318 , 0.18518519 , 0.18718719,
  0.18918919 , 0.19119119 , 0.19319319 , 0.1951952  , 0.1971972  , 0.1991992,
  0.2012012  , 0.2032032  , 0.20520521 , 0.20720721 , 0.20920921 , 0.21121121,
  0.21321321 , 0.21521522 , 0.21721722 , 0.21921922 , 0.22122122 , 0.22322322,
  0.22522523 , 0.22722723 , 0.22922923 , 0.23123123 , 0.23323323 , 0.23523524,
  0.23723724 , 0.23923924 , 0.24124124 , 0.24324324 , 0.24524525 , 0.24724725,
  0.24924925 , 0.25125125 , 0.25325325 , 0.25525526 , 0.25725726 , 0.25925926,
  0.26126126 , 0.26326326 , 0.26526527 , 0.26726727 , 0.26926927 , 0.27127127,
  0.27327327 , 0.27527528 , 0.27727728 , 0.27927928 , 0.28128128 , 0.28328328,
  0.28528529 , 0.28728729 , 0.28928929 , 0.29129129 , 0.29329329 , 0.2952953,
  0.2972973  , 0.2992993  , 0.3013013  , 0.3033033  , 0.30530531 , 0.30730731,
  0.30930931 , 0.31131131 , 0.31331331 , 0.31531532 , 0.31731732 , 0.31931932,
  0.32132132 , 0.32332332 , 0.32532533 , 0.32732733 , 0.32932933,  0.33133133,
  0.33333333 , 0.33533534 , 0.33733734 , 0.33933934 , 0.34134134,  0.34334334,
  0.34534535 , 0.34734735 , 0.34934935 , 0.35135135 , 0.35335335,  0.35535536,
  0.35735736 , 0.35935936 , 0.36136136 , 0.36336336 , 0.36536537,  0.36736737,
  0.36936937 , 0.37137137 , 0.37337337 , 0.37537538 , 0.37737738,  0.37937938,
  0.38138138 , 0.38338338 , 0.38538539 , 0.38738739  ,0.38938939,  0.39139139,
  0.39339339 , 0.3953954  , 0.3973974  , 0.3993994  , 0.4014014 ,  0.4034034,
  0.40540541 , 0.40740741 , 0.40940941 , 0.41141141 , 0.41341341,  0.41541542,
  0.41741742 , 0.41941942 , 0.42142142 , 0.42342342 , 0.42542543,  0.42742743,
  0.42942943 , 0.43143143 , 0.43343343 , 0.43543544 , 0.43743744,  0.43943944,
  0.44144144 , 0.44344344 , 0.44544545 , 0.44744745 , 0.44944945,  0.45145145,
  0.45345345 , 0.45545546 , 0.45745746 , 0.45945946 , 0.46146146,  0.46346346,
  0.46546547 , 0.46746747 , 0.46946947 , 0.47147147 , 0.47347347,  0.47547548,
  0.47747748 , 0.47947948 , 0.48148148 , 0.48348348 , 0.48548549,  0.48748749,
  0.48948949 , 0.49149149 , 0.49349349 , 0.4954955  , 0.4974975 ,  0.4994995,
  0.5015015  , 0.5035035  , 0.50550551 , 0.50750751 , 0.50950951,  0.51151151,
  0.51351351 , 0.51551552 , 0.51751752 , 0.51951952 , 0.52152152,  0.52352352,
  0.52552553 , 0.52752753 , 0.52952953 , 0.53153153 , 0.53353353,  0.53553554,
  0.53753754 , 0.53953954 , 0.54154154 , 0.54354354 , 0.54554555,  0.54754755,
  0.54954955 , 0.55155155 , 0.55355355 , 0.55555556 , 0.55755756,  0.55955956,
  0.56156156 , 0.56356356 , 0.56556557 , 0.56756757 , 0.56956957,  0.57157157,
  0.57357357 , 0.57557558,  0.57757758 , 0.57957958 , 0.58158158,  0.58358358,
  0.58558559 , 0.58758759 , 0.58958959 , 0.59159159 , 0.59359359,  0.5955956,
  0.5975976  , 0.5995996,   0.6016016  , 0.6036036  , 0.60560561,  0.60760761,
  0.60960961 , 0.61161161,  0.61361361 , 0.61561562 , 0.61761762,  0.61961962,
  0.62162162 , 0.62362362,  0.62562563 , 0.62762763 , 0.62962963,  0.63163163,
  0.63363363 , 0.63563564 , 0.63763764 , 0.63963964 , 0.64164164,  0.64364364,
  0.64564565 , 0.64764765 , 0.64964965 , 0.65165165 , 0.65365365,  0.65565566,
  0.65765766 , 0.65965966 , 0.66166166 , 0.66366366 , 0.66566567,  0.66766767,
  0.66966967 , 0.67167167 , 0.67367367 , 0.67567568 , 0.67767768,  0.67967968,
  0.68168168 , 0.68368368 , 0.68568569 , 0.68768769 , 0.68968969,  0.69169169,
  0.69369369 , 0.6956957  , 0.6976977  , 0.6996997  , 0.7017017 ,  0.7037037,
  0.70570571 , 0.70770771 , 0.70970971 , 0.71171171 , 0.71371371,  0.71571572,
  0.71771772 , 0.71971972 , 0.72172172 , 0.72372372 , 0.72572573,  0.72772773,
  0.72972973 , 0.73173173 , 0.73373373 , 0.73573574 , 0.73773774,  0.73973974,
  0.74174174 , 0.74374374 , 0.74574575 , 0.74774775 , 0.74974975,  0.75175175,
  0.75375375 , 0.75575576 , 0.75775776 , 0.75975976 , 0.76176176,  0.76376376,
  0.76576577 , 0.76776777 , 0.76976977 , 0.77177177 , 0.77377377,  0.77577578,
  0.77777778 , 0.77977978 , 0.78178178 , 0.78378378 , 0.78578579,  0.78778779,
  0.78978979 , 0.79179179 , 0.79379379 , 0.7957958  , 0.7977978 ,  0.7997998,
  0.8018018  , 0.8038038  , 0.80580581 , 0.80780781 , 0.80980981,  0.81181181,
  0.81381381 , 0.81581582 , 0.81781782 , 0.81981982 , 0.82182182,  0.82382382,
  0.82582583  ,0.82782783 , 0.82982983 , 0.83183183 , 0.83383383,  0.83583584,
  0.83783784  ,0.83983984 , 0.84184184 , 0.84384384 , 0.84584585,  0.84784785,
  0.84984985  ,0.85185185 , 0.85385385 , 0.85585586 , 0.85785786,  0.85985986,
  0.86186186 , 0.86386386 , 0.86586587 , 0.86786787 , 0.86986987,  0.87187187,
  0.87387387 , 0.87587588 , 0.87787788 , 0.87987988 , 0.88188188,  0.88388388,
  0.88588589 , 0.88788789 , 0.88988989 , 0.89189189 , 0.89389389,  0.8958959,
  0.8978979  , 0.8998999  , 0.9019019  , 0.9039039  , 0.90590591,  0.90790791,
  0.90990991 , 0.91191191 , 0.91391391 , 0.91591592 , 0.91791792,  0.91991992,
  0.92192192 , 0.92392392 , 0.92592593 , 0.92792793 , 0.92992993,  0.93193193,
  0.93393393 , 0.93593594 , 0.93793794 , 0.93993994 , 0.94194194,  0.94394394,
  0.94594595,  0.94794795 , 0.94994995 , 0.95195195 , 0.95395395,  0.95595596,
  0.95795796 , 0.95995996 , 0.96196196 , 0.96396396 , 0.96596597,  0.96796797,
  0.96996997 , 0.97197197 , 0.97397397 , 0.97597598 , 0.97797798,  0.97997998,
  0.98198198,  0.98398398 , 0.98598599,  0.98798799 , 0.98998999,  0.99199199,
  0.99399399,  0.995996,    0.997998,    1.};
  
  summary << "\n [";
  for(int j = 0; j<1000; j++){
    summary << setprecision(13) << f(X[j]) << ", ";
  }
  summary << "]";*/





