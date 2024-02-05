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

//Description of the function to approximate
struct Univariate_functionDescription {
  string functionString; // With the unknown x
  string derivedFunctionString; // Used for Remez 1+, derived of the function
  double (*f)(double x); // Must correspond to the string version
} functionDescription;

//Parameters used by both Remez I and Remez I+
struct Univariate_RemezIParameters {
  struct Univariate_functionDescription fdesc; //description of the function to approximate
  int approximationDegree; //a_i with i in [0 ; approximationDegree]
  int sizeD0; //Advised to be >= to approximationDegree+2 for quickest approximation
  int maxNbTurns; //Failsafe
  std::vector<double> bornersVar;
  bool remezPlus; // Warning, if remezPlus, exchange is ignored
  bool exchange;  // Warning, if exchange, sizeD0 = approximationDegree+2
  double approximationResult;
  double approximationPoints; 
} Univariate_RemezIParameters;

//Results of RemezI univariate for comparisons and graphs
struct Univariate_RemezIResult {
  int nbTurns;
  std::vector<double> errorStep1;
  std::vector<double> errorStep2;
  bool remezPlus;
  bool exchange;
  string path;
} Univariate_RemezIResult;

/*-----------------------------------------------------------------*/
/* FUNCTIONS DESCRIPTIONS */
/*-----------------------------------------------------------------*/

//Creates the object Univariate_functionDescription.
struct Univariate_functionDescription initialize_fdesc();
//Creates the object Univariate_RemezIParameters.
struct Univariate_RemezIParameters initialize_remezdesc(struct Univariate_functionDescription fdesc);

//Remez 1 depending of remez description for function descripted.
struct Univariate_RemezIResult univariateRemez(struct Univariate_RemezIParameters remezdesc);

//Makes the path in which results will be written.
string createPath(struct Univariate_RemezIParameters remezdesc);
//Writes the introduction of the summary file.
void write_introduction(string path, struct Univariate_RemezIParameters remezdesc);
//Writes D in summary file
void writeD(std::vector<double> D, string path, struct Univariate_RemezIParameters remezdesc, int turn);
//Writes the polynom A
void writePolynom(std::vector<double> A, string path, struct Univariate_RemezIParameters remezdesc);
//Writes the error of Step2
void writeErrorStep2(double errorStep2, string path, struct Univariate_RemezIParameters remezdesc);

//Makes the model of Step 2 
void createModelStep2(std::vector<double> A, struct Univariate_RemezIParameters remezdesc, string path);
//Gets Step2 result
double parserResultStep2(string path);
//Checks if the point is already in D
bool isItANewPoint(double point, std::vector<double> D, double approx);
//Makes the model of Step 2+
void createModelStep2Plus(std::vector<double> A, struct Univariate_RemezIParameters remezdesc, string path);

//Creates the first discretized set
std::vector<double> step0(struct Univariate_RemezIParameters remezdesc);
//Gets the best polynomial for x in D
std::vector<double> step1(std::vector<double> D, int degree, string path);
//Gets point where the error is at its peak
double step2(std::vector<double> A, struct Univariate_RemezIParameters remezdesc, string path);
//Creates the graph of convergence of the two errors
void grapheConvergenceRemez1(std::vector<double> errorStep1, std::vector<double> errorStep2, int nbTurns, string path);
//Creates the graph of error between f(x) and P(a,x)
void grapheError(std::vector<double> A, string function, int turn, int degree, std::vector<double> borners, double errorMax, string path);
//Creates the graph with the logarithmic convergence of the distance between the two errors
void logGrapheConvergence(std::vector<double> errorStep1, std::vector<double> errorStep2, string path);
//Function of exchanges
std::vector<double> exchange(std::vector<double> D, std::vector<double> A, double newPoint);
//Gets all extremum of the error function
std::vector<double> step2Plus(std::vector<double> A, struct Univariate_RemezIParameters remezdesc, string path);
std::vector<double> parserResultStep2Plus(string path);
//Checks if there is any new points in points
bool isThereNewPoints(std::vector<double> D, std::vector<double> points, double approx);
//Adds every new points in D and return the set
std::vector<double> addNewPoints(std::vector<double> D, std::vector<double> points, double approx);
//Return the maximal error in the extremum found by step2Plus
double errorMax(std::vector<double> newPoints, std::vector<double> A);
//Makes a graph with all the distance convergence
void grapheConvergenceComparisons(std::vector<struct Univariate_RemezIResult>remezResults);

/*-----------------------------------------------------------------*/
/* FUNCTIONS TO APPROXIMATE */
/*-----------------------------------------------------------------*/
double f(double x) {
  return cos(x)+sin(5*x)+pow(x,2);
}

double P(std::vector<double> A, double x){
  double y = 0;
  for(int i = 0 ; i < A.size(); i++){
    y += A[i] * pow(x, i);
  }
  return y;
}

double error(std::vector<double> A, double x){
  return abs(P(A, x) - f(x));
}


/*-----------------------------------------------------------------*/
/* PARAMETERS AND INPUTS */
/*-----------------------------------------------------------------*/

struct Univariate_functionDescription initialize_fdesc(){
  struct Univariate_functionDescription fdesc;
  fdesc.functionString = "cos(x)+sin(5*x)+x^2";
  fdesc.derivedFunctionString = "-sin(x)+5*cos(5*x)+2*x";
  fdesc.f = (*f);
  return fdesc;
}

struct Univariate_RemezIParameters initialize_remezdesc(struct Univariate_functionDescription fdesc){
  struct Univariate_RemezIParameters remezdesc;
  remezdesc.approximationDegree = 10;
  remezdesc.sizeD0 = 12; //Warning, if exchange, sizeD0 value will be replaced by degree+2
  remezdesc.maxNbTurns = 100;
  remezdesc.approximationResult = 1.e-10;
  remezdesc.approximationPoints = 1.e-15;
  remezdesc.bornersVar = {-1 , 1}; //Warning, remember to put in crescent order
  remezdesc.fdesc = fdesc;
  remezdesc.remezPlus = false ;
  remezdesc.exchange = false ; // Warning, if remezPlus, exchange must be ignored.
  return remezdesc;
}

/*-----------------------------------------------------------------*/
/* MAIN AND COMMON FUNCTIONS */
/*-----------------------------------------------------------------*/
int main(int argc, char** argv) { 
  //initialize fdesc and remezdesc
  struct Univariate_functionDescription fdesc;
  fdesc = initialize_fdesc();
  struct Univariate_RemezIParameters remezdesc;
  remezdesc = initialize_remezdesc(fdesc);
  
  //execute remez following remezdesc for fdesc in order classical, exchange and plus
  struct Univariate_RemezIResult resultBasic = univariateRemez(remezdesc);
  remezdesc.exchange = true;
  struct Univariate_RemezIResult resultExchange = univariateRemez(remezdesc);
  remezdesc.remezPlus = true;
  remezdesc.exchange = false;
  struct Univariate_RemezIResult resultPlus = univariateRemez(remezdesc);
  
  //create the graph for convergence comparisons
  std::vector<struct Univariate_RemezIResult>remezResults = {};
  remezResults.push_back(resultBasic);
  remezResults.push_back(resultExchange);
  remezResults.push_back(resultPlus);
  grapheConvergenceComparisons(remezResults);
}

struct Univariate_RemezIResult univariateRemez(struct Univariate_RemezIParameters remezdesc){
  //set the summary file
  string path = createPath(remezdesc);
  ofstream summary(path + "summary/" + remezdesc.fdesc.functionString + ".txt");
  write_introduction(path, remezdesc);
  //initializing variables
  struct Univariate_RemezIResult result;
  std::vector<double> errorStep1 = {};
  std::vector<double> errorStep2 = {};
  int turn = 1;
  std::vector<double> A = {}; //Approximation coefficient
  std::vector<double> D = {}; //Discretized set of points
  bool itIsNotOver = true ; //Stopping condition
  
  double newPoint; //Used only if !(remezPlus)
  std::vector<double> newPoints = {}; //Used only if remezPlus
  
  //With exchange, |D| = n+2 at all turns
  if (remezdesc.exchange && !(remezdesc.remezPlus)){
    remezdesc.sizeD0 = remezdesc.approximationDegree + 2;
  }
  
  //begins
  D = step0(remezdesc); //Creates first discretization
  writeD(D, path, remezdesc, turn);
  while(itIsNotOver){
    A = step1(D, remezdesc.approximationDegree, path); //Search a polynomial for discretized set
    writePolynom(A, path, remezdesc);
    errorStep1.push_back(A[remezdesc.approximationDegree+1]); //Will rise toward e*
    A.pop_back();
    if(remezdesc.remezPlus){
      newPoints = step2Plus(A, remezdesc, path); //finds all extremum of error function
      itIsNotOver = isThereNewPoints(newPoints, D, remezdesc.approximationPoints);
      D = addNewPoints(newPoints, D, remezdesc.approximationPoints); //Adds all new points in D
      errorStep2.push_back(errorMax(newPoints, A));
    } else {
      newPoint = step2(A, remezdesc, path); //Search where biggest error for polynomial
      itIsNotOver = isItANewPoint(newPoint, D, remezdesc.approximationPoints);
      if (remezdesc.exchange){
        D = exchange(D, A, newPoint); //With exchange, always n+2 point, so exchange a point of D
      }
      else{
        D.push_back(newPoint); //Otherwhise, add the new point to discretized set
      }
      errorStep2.push_back(error(A, newPoint)); //Will lower toward e*
    }
    writeErrorStep2(errorStep2[errorStep2.size()-1], path, remezdesc);
    turn = turn + 1;
    grapheError(A, remezdesc.fdesc.functionString, turn, remezdesc.approximationDegree, remezdesc.bornersVar, errorStep2[errorStep2.size()-1], path);
    //Checking if a new turn is needed
    itIsNotOver = itIsNotOver && turn < remezdesc.maxNbTurns;
    itIsNotOver = itIsNotOver && ((errorStep2[errorStep2.size()-1]-errorStep1[errorStep1.size()-1]) > remezdesc.approximationResult);
    writeD(D, path, remezdesc, turn);
  }
  grapheConvergenceRemez1(errorStep1, errorStep2, turn, path); //Prints the convergence of errors
  logGrapheConvergence(errorStep1, errorStep2, path); //Prints the convergence of errors distance
  result.errorStep1 = errorStep1;
  result.errorStep2 = errorStep2;
  result.nbTurns = turn - 1;
  result.exchange = remezdesc.exchange;
  result.remezPlus = remezdesc.remezPlus;
  result.path = path; 
  return result;
}

std::vector<double> step0(struct Univariate_RemezIParameters remezdesc){
  double up = remezdesc.bornersVar[1];
  double low = remezdesc.bornersVar[0];
  double step = (up-low)/(remezdesc.sizeD0-1);
  double temp = low;
  std::vector<double> D0 = {low};
  for(int i = 0; i < remezdesc.sizeD0-2; i++) {
    temp = temp + step;
    D0.push_back(temp);
  }
  D0.push_back(up);
  return D0;
}

std::vector<double> step1(std::vector<double> D, int degree, string path){
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
   for(int i = 0; i<=degree; i++){
     mysoplex.addColReal(LPColReal(0, dummycol, infinity, -infinity));
   }
   mysoplex.addColReal(LPColReal(1, dummycol, infinity, -infinity));
   /* Constraints */
   double x;
   Real r;
   for(int i = 0 ; i<D.size(); i++){
     DSVectorReal row1(0);
     for(int j = 0 ;  j<=degree; j++){
       x = pow(D[i], j);
       r = x;
       row1.add(j, r);
     }
     r = 1;
     row1.add(degree+1, r);
     r = f(D[i]);
     mysoplex.addRowReal(LPRowReal(r, row1, infinity));
   }
   
   for(int i = 0 ; i<D.size(); i++){
     DSVectorReal row1(0);
     for(int j = 0 ;  j<=degree; j++){
       x = pow(D[i], j);
       r = -x;
       row1.add(j, r);
     }
     r = 1;
     row1.add(degree+1, r);
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
/* FUNCTIONS REMEZ 1 - UNIVARIE */
/*-----------------------------------------------------------------*/

bool isItANewPoint(double point, std::vector<double> D, double approx){
  bool tooClose = false;
  for(int i=0; i<D.size(); i++){
    if(D.at(i)>=point-approx && D.at(i)<=point+approx){ tooClose = true; }
  }
  return !(tooClose);
}

void createModelStep2(std::vector<double> A, struct Univariate_RemezIParameters remezdesc, string path){
  ofstream step2(path + "models/step2_model.txt");
  step2 << "Variables\n";
  step2 << "  x in [" << remezdesc.bornersVar[0] <<" , " << remezdesc.bornersVar[1] <<"];\n";
  step2 << "\nMinimize\n";
  step2 << "  -abs(" ;
  step2 << "(" << setprecision(13) << A[0] << ") +";
  for(int i = 1; i<=remezdesc.approximationDegree; i++){
    step2 << "((x^" << i << ")*(" << setprecision(13) << A[i] << ")) +";
  }
  step2 << "(-(" << remezdesc.fdesc.functionString << "))); \n";

  step2 << "\nConstraints\n";
  step2 << "end";
  step2.close();
}

double parserResultStep2(string path){
  ifstream resultStep2(path + "stepResult/step2_result.txt");
  double result = 0.0;
  if(resultStep2) {
    string str;
    resultStep2 >> str;
    resultStep2 >> str;
    resultStep2 >> str;
    str.erase(std::remove(str.begin(), str.end(), '>'), str.end());
    str.erase(std::remove(str.begin(), str.end(), ')'), str.end());
    istringstream  istr(str);
    istr >> result;
  }
  else {
    cout << "ERROR FILE NOT OPENING" << endl;
  }
  return result;
}


double step2(std::vector<double> A, struct Univariate_RemezIParameters remezdesc, string path){
  createModelStep2(A, remezdesc, path);
  string fullPathToModel = path + "models/step2_model.txt";
  ofstream step2R(path + "stepResult/step2_result.txt");

  System sys(fullPathToModel.c_str());
  DefaultOptimizer optimizer(sys,1e-09, 1e-10);
  optimizer.optimize(sys.box);
  optimizer.report();
  
  step2R << "minimizer: " << optimizer.get_loup_point() << endl;
  step2R << "uplo " << optimizer.get_uplo() << "\n" << endl;
  step2R << "loup " << optimizer.get_loup() << endl;
  step2R << optimizer.get_data();
  step2R.close();
  return parserResultStep2(path);
  
}

std::vector<double> exchange(std::vector<double> D, std::vector<double> A, double newPoint){
  std::vector<double> Di = {};
  int i = 0;
  while(!(newPoint >= D[i] && newPoint <= D[i+1])){
    Di.push_back(D[i]);
    i = i + 1;
  }
  double error1, error2, errorNP;
  errorNP = P(A,newPoint) - f(newPoint);
  error1 = P(A,D[i]) - f(D[i]);
  error2 = P(A,D[i+1]) - f(D[i+1]);
  if(error1<0 && errorNP<0 || error1>0 && errorNP>0 ){
    Di.push_back(newPoint);
    Di.push_back(D[i+1]);
  } else {
    Di.push_back(D[i]);
    Di.push_back(newPoint);
  }
  i = i + 2;
  for(int j = i; j< D.size(); j++){
    Di.push_back(D[j]);
  }
  return Di;
}


/*-----------------------------------------------------------------*/
/* FUNCTIONS REMEZ 1+ - UNIVARIE */
/*-----------------------------------------------------------------*/

std::vector<double> step2Plus(std::vector<double> A, struct Univariate_RemezIParameters remezdesc, string path){
  ofstream step2RPlus(path + "/stepResult/step2Plus_result.txt");
  createModelStep2Plus(A, remezdesc, path);
  string fullpath = path + "models/step2Plus_model.txt";
  System system(fullpath.c_str());

  DefaultSolver solver(system,1e-09, 1e-15);
  solver.solve(system.box);
  solver.report();
  step2RPlus << solver.get_data() << endl;
  return parserResultStep2Plus(path);
}

void createModelStep2Plus(std::vector<double> A, struct Univariate_RemezIParameters remezdesc, string path){
  ofstream step2(path + "models/step2Plus_model.txt");
  step2 << "Variables\n";
  step2 << "  x in [" << remezdesc.bornersVar[0] <<" , " << remezdesc.bornersVar[1] <<"];\n";
  step2 << "\nConstraints\n";
  step2 << A[1] << " + ";
  for(int i = 2 ; i<= remezdesc.approximationDegree; i++){
    step2 << "((" << A[i] <<"*" << i <<")*x" << "^" << i-1 << ")+";
  }
  step2 << "-(" << remezdesc.fdesc.derivedFunctionString << ") = 0 ;\n";
  step2 << "end";
  step2.close();
}

std::vector<double> parserResultStep2Plus(string path){
  ifstream resultStep2(path + "stepResult/step2Plus_result.txt");
  std::vector<double> R = {};
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
      istringstream  istr3(smaller);
      istr3 >> small;
      R.push_back(small);
      resultStep2 >> str;
    }
  }
  else {
    cout << "ERROR FILE NOT OPENING" << endl;
  }
  return R;
}

bool isThereNewPoints(std::vector<double> D, std::vector<double> points, double approx){
  bool thereIs = false;
  for(int i = 0; i<points.size(); i++){
    thereIs = isItANewPoint(points[i], D, approx) || thereIs;
  }
  return thereIs;
}

std::vector<double> addNewPoints(std::vector<double> D, std::vector<double> points, double approx){
  std::vector<double> Di = {};
  for(int i = 0; i<points.size(); i++){
    if(isItANewPoint(points[i], D, approx)){
      Di.push_back(points[i]);
    }
  }
  
  for(int i = 0; i<D.size(); i++){
    Di.push_back(D[i]);
  }
  return Di;
}

double errorMax(std::vector<double> newPoints, std::vector<double> A){
  double errorMax = 0;
  for(int i=0; i<newPoints.size(); i++){
    errorMax = std::max(errorMax, error(A, newPoints[i]));
  }
  return errorMax;
}

/*-----------------------------------------------------------------*/
/* FUNCTIONS GRAPHS CREATIONS*/
/*-----------------------------------------------------------------*/

void grapheError(std::vector<double> A, string function, int turn, int degree, std::vector<double> borners, double errorMax, string path){
  string graphName = function + to_string(turn) + "AproximationRemez1";
  graphName = path + graphName;
  graphName += ".pdf";
  Vec x = linspace(borners[0], borners[1], 1000);
  
  std::vector<double> e = {};
  
  for(int i=0; i<x.size(); i++){
    e.push_back(P(A,x[i])-f(x[i]));
  }
  
  Plot2D plot;
  plot.xlabel("x");
  plot.ylabel("y");
  plot.xrange(borners[0], borners[1]);
  plot.yrange(-errorMax , errorMax);
  
  plot.legend()
    .atOutsideBottom()
    .displayHorizontal()
    .displayExpandWidthBy(2);
  plot.drawCurve(x, e).label("error");
  Figure fig = {{plot}};
  Canvas canvas = {{fig}};
  canvas.save(graphName);
}

void grapheConvergenceRemez1(std::vector<double> errorStep1, std::vector<double> errorStep2, int nbTurns, string path){
  string graphName = "graph/GrapheConvergence";
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

void logGrapheConvergence(std::vector<double> errorStep1, std::vector<double> errorStep2, string path){
  ofstream convText(path+"summary/convergence.txt");
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
  canvas.save(path + "graph/GrapheConvergenceLogarithmic.pdf");
}


void grapheConvergenceComparisons(std::vector<struct Univariate_RemezIResult>remezResults){
  string graphName = "../resultsOfSteps/GrapheConvergenceAll.pdf";
  string output = "";
  std::vector<double> distanceTemp = {};
  std::vector<std::vector<double>> distance = {};
  double dist = 0;
  std::vector<double> conv = {};
  for(int j = 0; j<remezResults.size(); j++){
    for(int i = 0; i<remezResults[j].errorStep1.size(); i++){
      dist = sqrt(pow((remezResults[j].errorStep2[i]-remezResults[j].errorStep1[i]),2));
      distanceTemp.push_back(dist);
      conv.push_back(dist);
    }
    distance.push_back(distanceTemp);
    distanceTemp = {};
  }
  int maxNbTurns = 0;
  for(int j = 0; j<remezResults.size(); j++){
    maxNbTurns = std::max(maxNbTurns, remezResults[j].nbTurns);
  }
  vector<vector<double>> x = {};
  vector<double> xi = {};
  for(int i=0; i<remezResults.size(); i++){
    for(int j=1; j<remezResults[i].nbTurns+1; j++){
      xi.push_back(j);
    }
    x.push_back(xi);
    xi = {};
  }
  Plot2D plot;
  plot.xlabel("Turn");
  plot.ylabel("Convergence");
  plot.xrange(0., maxNbTurns-1);
  plot.yrange(*min_element(conv.begin(), conv.end()), *max_element(conv.begin(), conv.end()));

  plot.legend()
    .atOutsideBottom()
    .displayHorizontal()
    .displayExpandWidthBy(2);
  plot.ytics().logscale(2);
  for(int j = 0; j<remezResults.size(); j++){
    string curveTitle = "convergence logarithmic Remez1";
    if(remezResults[j].remezPlus){
     curveTitle += "Plus";
    }
    if(remezResults[j].exchange){
     curveTitle += " version Exchange";
    }
    plot.drawCurve(x[j], distance[j])
        .label(curveTitle);

  }
  Figure fig = {{plot}};
  Canvas canvas = {{fig}};
  canvas.save(graphName);
}

/*-----------------------------------------------------------------*/
/* FUNCTIONS FILE MANAGEMENT*/
/*-----------------------------------------------------------------*/

void write_introduction(string path, struct Univariate_RemezIParameters remezdesc){
  std::ofstream summary;
  summary.open(path + "summary/" + remezdesc.fdesc.functionString + ".txt", std::ios_base::app);
  summary << "Welcome to the summary of ";
  if (remezdesc.remezPlus){
    summary << "Remez1+ ";
  } else {
    if(remezdesc.exchange){
      summary << "Exchange ";
    } else {
      summary << "Basic";
    }
    summary << "Remez 1 ";
  }
  summary << "for univariate polynomial approximation. \n";
  summary << "With the following parameters : \n";
  summary << "degree of approximation = " << remezdesc.approximationDegree << "\n";
  summary << "size of first discretization =" << remezdesc.sizeD0 << "\n";
  summary << "maximum number of turns = " << remezdesc.maxNbTurns << "\n";
  summary << "approximation of error = " << remezdesc.approximationResult << "\n";
  summary << "approximation of points = " << remezdesc.approximationPoints << "\n";
  summary<<"x in ["<<remezdesc.bornersVar[0]<<","<<remezdesc.bornersVar[1]<<"]"<<"\n";
  summary << "For the function " << remezdesc.fdesc.functionString << ".\n";
}

string createPath(struct Univariate_RemezIParameters remezdesc){
  string path = "../resultsOfSteps/";
  if (remezdesc.remezPlus){
    path += "remez1+/";
  } else {
    path += "remez1/";
    if(remezdesc.exchange){
      path += "exchange/";
    } else {
      path += "basic/";
    }
  }
  return path;
}

void writeD(std::vector<double> D, string path, struct Univariate_RemezIParameters remezdesc, int turn){
std::ofstream summary;
  summary.open(path + "summary/" + remezdesc.fdesc.functionString + ".txt", std::ios_base::app);
  summary << "Step " << turn << ": \n";
  summary << "D"<<"=[" << D[0];
  for(int j = 1; j<remezdesc.sizeD0; j++){
    summary << ","  << D[j];
  }
  summary << "]\n";
}

void writePolynom(std::vector<double> A, string path, struct Univariate_RemezIParameters remezdesc){
  std::ofstream summary;
  summary.open(path + "summary/" + remezdesc.fdesc.functionString + ".txt", std::ios_base::app);
  summary << "Polynom = [" << A[0] ;
  for(int j = 1; j<=remezdesc.approximationDegree; j++){
    summary << "," << setprecision(13) << A[j] ;
  }
  summary << "] \n" ;
  summary << "Error step1 ="<<A[remezdesc.approximationDegree+1] << "\n";

}

void writeErrorStep2(double errorStep2, string path, struct Univariate_RemezIParameters remezdesc){
  std::ofstream summary;
  summary.open(path + "summary/" + remezdesc.fdesc.functionString + ".txt", std::ios_base::app);
  summary << "Error step2 ="<< setprecision(13) << errorStep2 << "\n";
}












