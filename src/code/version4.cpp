//Coralie Marchau - 2024
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

//Description of a multivariate polynomial
struct Multivariate_Polynomiale {
  std::vector<double> a; //Value of a
  std::vector<double>(*phi)(std::vector<double> x, int degree); //Each x corresponding to the multiplication with coefficient
  std::vector<std::vector<int>> (*derivatedPhi)(int nbX, int degree); //Power of each X for polynomial in same order as coefficient
  int degree; //Degree of the Polynomial
} MultivariatepolynomialDescription;

//Description of the function to approximate
struct Multivariate_functionDescription {
  string functionString; // With the unknown x1, x2, ..., xi, ...
  int nbX; //number of variables
  std::vector<string> derivedFunctionStrings; //derivated function by xi
  std::vector<std::vector<double>> domainX; //Domain of X
  std::vector<bool> notConstantDerived; //Is fixing a side cause f' to be constant
  double (*f)(std::vector<double> x); // Must correspond to the string version
} MultivariatefunctionDescription;

struct Multivariate_RemezIIParameters {
  struct Multivariate_Polynomiale poly;
  struct Multivariate_functionDescription fdesc; //description of the function to approximate
  int approximationDegree; //a_i with i in [0 ; approximationDegree], as many i as nbVar
  int maxNbTurns; //Failsafe
  double approximationResult; //parameter
  double approximationPoints; //parameter
  double approxError = 1.e-13; //If we reach this closeness between phase2 and phase3, we consider we are at optimality.
} Multivariate_RemezIIParameters;

//Results of RemezII multivariate for comparisons and graphs
struct Multivariate_RemezIIResult {
  int nbTurns; 
  int currentTurn;
  std::vector<double> errorphase2;
  std::vector<double> errorphase3;
  string path;
} Multivariate_RemezIIResult;

/*-----------------------------------------------------------------*/
/* FUNCTIONS DECLARATION */
/*-----------------------------------------------------------------*/
struct Multivariate_functionDescription initialize_fdesc(string fString, std::vector<string> df, double (*f)(std::vector<double> x), std::vector<std::vector<double>> domainX, std::vector<bool> nCD);
struct Multivariate_RemezIIParameters initialize_remezdesc(struct Multivariate_functionDescription fdesc, int degree, int nbTurns);
struct Multivariate_Polynomiale initialize_poly(std::vector<double>(*phi)(std::vector<double> x, int degree), std::vector<std::vector<int>> (*derivatedPhi)(int nbX, int degree));
std::vector<std::vector<double>> phase1(int nbX, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_functionDescription fdesc);
string createPath();
std::vector<double> phase2(std::vector<std::vector<double>> D, struct Multivariate_functionDescription fdesc, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_Polynomiale poly);
void write_introduction(string path, struct Multivariate_RemezIIParameters remezdesc);
void writeD(struct Multivariate_RemezIIParameters remezdesc, std::vector<std::vector<double>> D, string path);
std::vector<std::vector<double>> phase3(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path, int turn);
void createErrorGraph(struct Multivariate_RemezIIParameters remezdesc, string path);
void writeA(std::vector<double> a, string path, string functionString);
std::vector<std::vector<double>> phase1_random(int nbPoints, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_functionDescription fdesc);
void write_data_for_graphs(std::vector<double> A, std::vector<std::vector<double>> domainX, int turn);
std::vector<std::vector<double>> phase3_center(struct Multivariate_RemezIIParameters remezdesc, string path, int turn);
void createModelphase3(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path);
void createModelphase3_side(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc,double borner, string path, int var);
std::vector<std::vector<double>> phase3_side(struct Multivariate_RemezIIParameters remezdesc, string path, int turn);
std::vector<std::vector<double>> noCopyAdd(std::vector<std::vector<double>> v1, std::vector<std::vector<double>> v2);
double phase3Error(std::vector<std::vector<double>> newPoints, std::vector<double> a, struct Multivariate_RemezIIParameters remezdesc, string path);
void writePolynom(std::vector<double> A, string path, struct Multivariate_RemezIIParameters remezdesc);
void writeErrors(double error1, double error2, struct Multivariate_RemezIIParameters remezdesc, string path);
std::vector<std::vector<double>> phase3_centerDIRECT(struct Multivariate_RemezIIParameters remezdesc, string path, int turn);
void createModelphase3Error(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path);
void grapheConvergenceComparisons(struct Multivariate_RemezIIResult remezResults);
bool cheking_side_not_constant(struct Multivariate_RemezIIParameters remezdesc, int side, int var);
void replaceAll(std::string& str, const std::string& from, const std::string& to);
bool checkEquioscillationOptimal(struct Multivariate_RemezIIParameters remezdesc, string path);
bool isItOverTest(struct Multivariate_RemezIIParameters remezdesc, string path, double errorphase2, double errorphase3);
void maximaCheck(struct Multivariate_RemezIIParameters remezdesc, string path, double errorphase2, double errorphase3);
std::vector<std::vector<double>> phase3_borders(struct Multivariate_RemezIIParameters remezdesc, string path, int turn);
void createModelphase3_border(struct Multivariate_RemezIIParameters remezdesc, int borner1, int borner2, string path, int var);
void write_nbMaxima(struct Multivariate_RemezIIParameters remezdesc, string path, std::vector<std::vector<double>> newPoints, double errorphase2, double errorphase3);
bool borderNotConstant(struct Multivariate_RemezIIParameters remezdesc, int var, int borner1, int borner2);
void newton(std::vector<double>  point, struct Multivariate_RemezIIParameters remezdesc , string path, double e);
std::vector<std::vector<double>> correctNewPoints(std::vector<std::vector<double>>  points, struct Multivariate_RemezIIParameters remezdesc , string path, double e);
/*-----------------------------------------------------------------*/
/* FUNCTIONS TO APPROXIMATE */
/*-----------------------------------------------------------------*/
//Didactical Examples
double f_testpow4(std::vector<double> x) {
  return pow(x[0],4) + pow(x[0],3)*pow(x[1],3);
}

struct Multivariate_functionDescription initializeFPow4(){
  return initialize_fdesc("x1^4 + (x2^3)*(x1^3)",{"4*(x1^3) + 3*(x1^2)*(x2^3)","3*(x2^2)*(x1^3)"},f_testpow4, {{-1,1},{-1,1}}, {true, true, true, true});
}


double f0(std::vector<double> x) {
  return pow(x[0],2) + pow(x[1],2);
}

struct Multivariate_functionDescription initializeF0(){
  return initialize_fdesc("x1^2 + x2^2",{"2*x1","2*x2", "2", "2", "0"},f0, {{0,1},{0,1}}, {true, true, true, true});
}


double fExp(std::vector<double> x) {
  return exp(x[0])+x[1];
}
struct Multivariate_functionDescription initializeFExp(){
  return initialize_fdesc("exp(x1)+x2",{"exp(x1)", "1"},fExp, {{0,1},{0,1}}, {false, false, true, true});
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
  return initialize_fdesc("cos(x3)*(1+x1)^(x2)",{"x2*cos(x3)*(1+x1)^(x2-1)"," cos(x3)*ln(1+x1)*(1+x1)^(x2)","(-sin(x3)*(1+x1)^(x2))"},f3, {{0,1},{1, 2},{0,1}}, {false, true, true, true, false, true});
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
  return initialize_fdesc("exp(x1^2 + x1*x2)",{"((2*(x1*exp((x1^2+(x1*x2)))))+(x2*exp((x1^2+(x1*x2)))))", "(x1*exp((x1^2+(x1*x2))))"}, f5, {{-1,1},{-1,1}}, {true, true, true, true});
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
  return initialize_fdesc("abs(ln(x1*x2+1)/(x1+0.5))*x2^((x3+1)/2)",{"((x2*(((sign((ln((1+(x1*x2)))/(0.5+x1)))*exp(((0.5+(0.5*x3))*ln(x2))))/(0.5+x1))/(1+(x1*x2))))-(((sign((ln((1+(x1*x2)))/(0.5+x1)))*exp(((0.5+(0.5*x3))*ln(x2))))*ln((1+(x1*x2))))/(0.5+x1)^2))","((x1*(((sign((ln((1+(x1*x2)))/(0.5+x1)))*exp(((0.5+(0.5*x3))*ln(x2))))/(0.5+x1))/(1+(x1*x2))))+(((0.5+(0.5*x3))*(abs((ln((1+(x1*x2)))/(0.5+x1)))*exp(((0.5+(0.5*x3))*ln(x2)))))/x2))","(0.5*((abs((ln((1+(x1*x2)))/(0.5+x1)))*exp(((0.5+(0.5*x3))*ln(x2))))*ln(x2)))"},f8, {{0,1},{0,1},{0,1}}, {false, true, false, false, true, false});
}


/*-----------------------------------------------------------------*/
/* BASICS */
/*-----------------------------------------------------------------*/

double P(struct Multivariate_Polynomiale A , std::vector<double> x){
  double y = 0;
  std::vector<double> p = A.phi(x, A.degree);
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
  if(x.size() == 2){
    for(int i = 0; i<=degree; i=i+1){
      for(int j = 0; j<=degree-i; j=j+1){
        p.push_back(pow(x[0],i)*pow(x[1],j));
      }
    }
  } else {
    if(x.size() == 3){
      for(int i = 0; i<=degree; i=i+1){
        for(int j = 0; j<=degree-i; j=j+1){
          for(int k = 0; k<=degree-i-j; k=k+1){
            p.push_back(pow(x[0],i)*pow(x[1],j)*pow(x[2],k));
          }
        }
      }
    }
  }
  return p;
}

std::vector<std::vector<int>> phiBasiqueDerivated(int nbX, int degree){
  std::vector<std::vector<int>> powOfX = {};
  std::vector<std::vector<int>> temp = {};
  std::vector<std::vector<int>>  prev = {std::vector<int>(nbX, 0.0)};
  if(nbX == 2){
    for(int i = 0; i<=degree; i=i+1){
      for(int j = 0; j<=degree-i; j=j+1){
        powOfX.push_back({i,j});
      }
    }
  }
  if(nbX == 3){
    for(int i = 0; i<=degree; i=i+1){
      for(int j = 0; j<=degree-i; j=j+1){
        for(int k = 0; k<=degree-i-j; k=k+1){
          powOfX.push_back({i,j,k});
        }
      }
    }
  }
  return powOfX;
}

std::vector<double> phiSecond(std::vector<double> x, int degree){
  std::vector<double> p = {};
  
  if(x.size() == 2){
    for(int i = 0; i <= degree ; i++){
      for(int j = 0; j <= degree ; j++){
        p.push_back(pow(x[0],i)*pow(x[1],j));
      }
    }
  } else {
    for(int i = 0; i <= degree ; i++){
      for(int j = 0; j <= degree ; j++){
        for(int k = 0; k <= degree ; k++ ){
          p.push_back(pow(x[0],i)*pow(x[1],j)*pow(x[2],k));
        }
      }
    }
  }
  
  return p;
}

std::vector<std::vector<int>> phiSecondDerivated(int nbX, int degree){
  std::vector<std::vector<int>> powOfX = {};
  if(nbX == 2){
    for(int i = 0; i <= degree ; i++){
      for(int j = 0; j <= degree ; j++){
        powOfX.push_back({i,j});
      }
    }
  } else {
    for(int i = 0; i <= degree ; i++){
      for(int j = 0; j <= degree ; j++){
        for(int k = 0; k <= degree ; k++ ){
          powOfX.push_back({i,j,k});
        }
      }
    }
  }
  
  return powOfX;
}

struct Multivariate_Polynomiale initializePhiSecond(){
  return initialize_poly(phiSecond, phiSecondDerivated);
}



struct Multivariate_Polynomiale initializePhiBasic(){
  return initialize_poly(phiBasique, phiBasiqueDerivated);
}

/*-----------------------------------------------------------------*/
/* INITIALIZES OBJECTS */
/*-----------------------------------------------------------------*/
struct Multivariate_functionDescription initialize_fdesc(string fString, std::vector<string> df, double (*f)(std::vector<double> x), std::vector<std::vector<double>> domainX, std::vector<bool> nCD){
  struct Multivariate_functionDescription fdesc;
  fdesc.functionString = fString;
  fdesc.derivedFunctionStrings = df;
  fdesc.f = f;
  fdesc.notConstantDerived = nCD;
  fdesc.domainX = domainX;
  fdesc.nbX = domainX.size();
  return fdesc;
}

struct Multivariate_RemezIIParameters initialize_remezdesc(struct Multivariate_functionDescription fdesc, int degree, int nbTurns){
  struct Multivariate_RemezIIParameters remezdesc;
  remezdesc.approximationDegree = degree;
  remezdesc.maxNbTurns = nbTurns;
  remezdesc.approximationResult = 1.e-9;
  remezdesc.approximationPoints = 1.e-12;
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
  result.errorphase2 = {};
  result.errorphase3 = {};
  result.path = path;
  return result;
}
/*-----------------------------------------------------------------*/
/* MAIN */
/*-----------------------------------------------------------------*/

int main(int argc, char** argv) { 
  //initialization
  int nbTurns = 23;
  int degree = 2;
  double errorphase3 = 0;
  double errorphase2 = 0; 
  string path = createPath();
  std::vector<std::vector<double>> newPoints;
  struct Multivariate_RemezIIParameters remezdesc;
  struct Multivariate_RemezIIResult results;
  results = initialize_results(nbTurns, path);
  remezdesc = initialize_remezdesc(initializeF4(), degree, nbTurns);
  remezdesc.poly = initializePhiSecond();
  remezdesc.poly.degree = degree;
  write_introduction(path, remezdesc);
  std::vector<double> phi;
  std::vector<std::vector<double>> D ;
  int nbPoints = 200; //minimum advised = a.size + 1
  int turn = 0; 
  bool notFinished = true;
        
  //phase 1 : first discretization
  D = noCopyAdd(phase1_random(5-nbPoints, remezdesc, remezdesc.fdesc), phase1(remezdesc.fdesc.nbX, remezdesc, remezdesc.fdesc));

  while(turn < nbTurns && notFinished){
    writeD(remezdesc, D, path);
    
    //phase 2 : search polynomial optimal for discretized set
    remezdesc.poly.a = phase2(D, remezdesc.fdesc, remezdesc, remezdesc.poly);
       
      //Taking the error out of the polynomial
    errorphase2 = remezdesc.poly.a[remezdesc.poly.a.size()-1];
    remezdesc.poly.a.pop_back();
    writePolynom(remezdesc.poly.a, path, remezdesc);
 
    //createErrorGraph(remezdesc, path); //c++ version of graph no longer used as less readable than python version
    write_data_for_graphs(remezdesc.poly.a, remezdesc.fdesc.domainX, turn+1); //write information for graphs python
    
    //phase 3 : searching points to growth of discretized set
    newPoints = phase3(remezdesc.poly.a, remezdesc, path, turn+1);
    errorphase3 = phase3Error(newPoints, remezdesc.poly.a, remezdesc, path);
    
    //phase 4 : adding the points in discretized set
    D = noCopyAdd(newPoints, D); //to avoid points 'coming back'
    writeErrors(errorphase2, errorphase3, remezdesc, path);
    results.errorphase2.push_back(errorphase2);
    results.errorphase3.push_back(errorphase3);
    
    cout << "\n---------------------------------------------------------------\n";
    notFinished = !(isItOverTest(remezdesc, path, errorphase2, errorphase3)); //needs amelioration, for now, check closeness of errors phase 2 and phase 3
    turn = turn + 1;
  }
  results.nbTurns = turn;
  //grapheConvergenceComparisons(results); //C++ graphe of convergence of distance errors phase 2 and 3
  //maximaCheck(remezdesc, path, errorphase2, errorphase3); //Search all points with an error between error phase 2 - delta and error phase 3
}

void maximaCheck(struct Multivariate_RemezIIParameters remezdesc, string path, double errorphase2, double errorphase3){
  std::vector<double> A = remezdesc.poly.a;
  double delta = 0.0; //delta must stay small
  ofstream phase3(path + "models/maximaCheck.txt");
  phase3 << setprecision(13);
  phase3 << "Variables\n";
  for(int i = 0; i<remezdesc.fdesc.domainX.size(); i++){
    phase3 << "  x"<< i+1 <<" in [" << remezdesc.fdesc.domainX[i][0] <<" , " << remezdesc.fdesc.domainX[i][1] <<"];\n";
  }
  
  phase3 << "\nConstraints\n";
  
  //error(a,x) >= error phase 2
  phase3 << "abs(";
  phase3 << remezdesc.fdesc.functionString << "-(";
  std::vector<std::vector<int>> p = remezdesc.poly.derivatedPhi(remezdesc.fdesc.domainX.size(), remezdesc.approximationDegree);
  for(int i = 0; i < A.size(); i++){
    phase3 << "+(" << A[i] << ")" ;
    for(int j = 0; j < p[i].size(); j++){
      if(p[i][j] != 0){
        phase3 << "*(x" << j+1 << "^(" << p[i][j] << "))";
      }
    }
  }
  phase3 << ")) >=" << errorphase2 - delta << ";\n";
  
  //error(a,x) <= error phase 3. This constraint is normally useless as error phase 3 should be the higher error of the domain.
 /* phase3 << "abs(";
  phase3 << remezdesc.fdesc.functionString << "-(";
  for(int i = 0; i < A.size(); i++){
    phase3 << "+(" << A[i] << ")" ;
    for(int j = 0; j < p[i].size(); j++){
      if(p[i][j] != 0){
        phase3 << "*(x" << j+1 << "^(" << p[i][j] << "))";
      }
    }
  }
  phase3 << ")) <=" << errorphase3 << ";\n";*/
  phase3 << "end";
  phase3.close(); 
  string fullpath = path + "models/maximaCheck.txt";
  System system(fullpath.c_str());
  DefaultSolver solver(system,1e-3);
  solver.solve(system.box);
  solver.report();
  std::ofstream data;
  data.open("data/extremas.txt",  std::ofstream::trunc);
  
  std::vector<double> pt = {};
  std::vector<std::vector<double>> newPoints = {};
  std::vector<std::vector<double>> points = {};
  std::vector<double> point= {};
  cout << "\n";
  for(int i = 0; i<solver.get_data().size(); i++){
    cout << "[";
    for(int j = 0; j<solver.get_data()[i].size(); j++){
      pt.push_back(solver.get_data()[i][j].lb());
      pt.push_back(solver.get_data()[i][j].ub());
      data << solver.get_data()[i][j].lb() << " " << solver.get_data()[i][j].ub() << " ";
      points.push_back(pt);
      p = {};
      cout << solver.get_data()[i][j].mid() << ",";
    }
    data << " \n";
    cout << "] ,";
    newPoints.push_back(point);
    point = {};
  }
  
  cout << points.size() << endl;
}

bool isItOverTest(struct Multivariate_RemezIIParameters remezdesc, string path, double errorphase2, double errorphase3){
  //Possibilities : error distance is too weak or optimality condition reached
  bool over = false; 
  //Check equioscillation optimal = not implemented, should give true if there is equioscillation
  //checkEquioscillationOptimal(remezdesc, path);
  if(false){
    over = true;
    cout << "end by optimality";
  } else {
    over = (errorphase3 - errorphase2 <= remezdesc.approxError);
    if(over){
      cout << "\nend by closeness of errors\n";
    }
  }
  return over; 
}


void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    if(from.empty())
        return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

/*-----------------------------------------------------------------*/
/* REMEZ II */
/*-----------------------------------------------------------------*/

//gets the corners of the domain
std::vector<std::vector<double>> phase1(int nbX, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_functionDescription fdesc){
  std::vector<std::vector<double>> D0 = {};
  std::vector<double> d = {};
  for(int i = 0; i < nbX; i ++){ //middle
    d.push_back(((fdesc.domainX[i][1] - fdesc.domainX[i][0])/2)+fdesc.domainX[i][0]);
  }
  D0.push_back(d);
  d = {};
  for(int i = 0; i < nbX; i ++){ //all min
    d.push_back(fdesc.domainX[i][0]);
  }
  D0.push_back(d);
  d = {};
  
  for(int i = 0; i < nbX; i ++){ //all max
    d.push_back(fdesc.domainX[i][1]);
  }
  D0.push_back(d);
  d = {};
  
  for(int j = 0; j < nbX; j ++){ 
    //j only max
    for(int i = 0; i < j; i ++){ 
      d.push_back(fdesc.domainX[i][0]);
    }
    d.push_back(fdesc.domainX[j][1]);
    for(int i = j+1; i < nbX; i ++){ 
      d.push_back(fdesc.domainX[i][0]);
    }
    D0.push_back(d);
    d = {};
    //all until j at max all after at min
    for(int i = 0; i < j; i ++){ 
      d.push_back(fdesc.domainX[i][1]);
    }
    d.push_back(fdesc.domainX[j][1]);
    for(int i = j+1; i < nbX; i ++){ 
      d.push_back(fdesc.domainX[i][0]);
    }
    D0.push_back(d);
    d = {};
    //all except j at max
    for(int i = 0; i < j; i ++){ 
      d.push_back(fdesc.domainX[i][1]);
    }
    d.push_back(fdesc.domainX[j][0]);
    for(int i = j+1; i < nbX; i ++){ 
      d.push_back(fdesc.domainX[i][1]);
    }
    D0.push_back(d);
    d = {};
  }
  
  return D0;
}

//Gets random points in domain
std::vector<std::vector<double>> phase1_random(int nbPoints, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_functionDescription fdesc){
  srand (42);
  std::vector<std::vector<double>> D0 = {};
  std::vector<double> temp = {};
  for(int i = 0; i<nbPoints; i++){
    for(int j = 0; j < fdesc.domainX.size(); j ++){ 
      temp.push_back(((double)rand()/RAND_MAX)*(fdesc.domainX[j][1] - fdesc.domainX[j][0])+fdesc.domainX[j][0]);
    }
    D0.push_back(temp);
    temp = {};
  }
  return D0;
}

//Search optimal polynomial for discretized set
std::vector<double> phase2(std::vector<std::vector<double>> D, struct Multivariate_functionDescription fdesc, struct Multivariate_RemezIIParameters remezdesc, struct Multivariate_Polynomiale poly){
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
   string fullPathToModel = path + "models/phase2Soplex.lp";
   mysoplex.writeFileReal(fullPathToModel.c_str(), NULL, NULL, NULL);
   SPxSolver::Status stat;
   DVectorRational prim(p[0].size()+1);
   DVectorRational dual(D.size()*2);
   stat = mysoplex.optimize();
   std::vector<double> primDouble;
   mysoplex.saveSettingsFile("soplex_param.set", true);
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

//Searches stationnary points to grow discretized set
std::vector<std::vector<double>> phase3(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path, int turn){
  std::vector<std::vector<double>> pointsSide = {};
  std::vector<std::vector<double>> pointsCenter = {};
  //Getting the points
  pointsCenter = phase3_center(remezdesc, path, turn);
  pointsSide = phase3_side(remezdesc, path, turn);
  if(remezdesc.fdesc.domainX.size() == 3){
    std::vector<std::vector<double>> pointsBorders = {};
    pointsBorders = phase3_borders(remezdesc, path, turn);
    pointsCenter.insert(pointsCenter.end(), pointsBorders.begin(), pointsBorders.end());
  }
  //Getting every points together
  pointsCenter.insert(pointsCenter.end(), pointsSide.begin(), pointsSide.end());
  return pointsCenter; 
}

std::vector<std::vector<double>> phase3_borders(struct Multivariate_RemezIIParameters remezdesc, string path, int turn){
  std::ofstream data;
  data.open("data/"+to_string(turn)+"points.txt",  std::ios_base::app);
  std::vector<std::vector<double>> newPoints = {};
  std::vector<double> point;
  string fullpath = path + "models/phase3_model_border.txt";
  //Seeing which border is checked
  for(int i = 0; i  < remezdesc.fdesc.domainX.size(); i++){
    if(borderNotConstant(remezdesc, i, 0, 0)){
      createModelphase3_border(remezdesc, 0, 0, path, i);
      System system1(fullpath.c_str());
      DefaultSolver solver1(system1,1e-10, 1e-15);
      solver1.solve(system1.box);
      solver1.report();
      for(int i = 0; i<solver1.get_data().size(); i++){
        for(int j = 0; j<solver1.get_data()[i].size(); j++){
          point.push_back(solver1.get_data()[i][j].mid());
          data << solver1.get_data()[i][j].mid() << " ";
        }
        data << "\n";
        newPoints.push_back(point);
        point = {};
      }
    }
    if(borderNotConstant(remezdesc, i, 1, 0)){
      createModelphase3_border(remezdesc, 1, 0, path, i);
      System system2(fullpath.c_str());
      DefaultSolver solver2(system2,1e-10, 1e-15);
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
    if(borderNotConstant(remezdesc, i, 0, 1)){
      createModelphase3_border(remezdesc, 0, 1, path, i);
      System system3(fullpath.c_str());
      DefaultSolver solver3(system3,1e-10, 1e-15);
      solver3.solve(system3.box);
      solver3.report();
      for(int i = 0; i<solver3.get_data().size(); i++){
        for(int j = 0; j<solver3.get_data()[i].size(); j++){
          point.push_back(solver3.get_data()[i][j].mid());
          data << solver3.get_data()[i][j].mid() << " ";
        }
        data << "\n";
        newPoints.push_back(point);
        point = {};
      }
    }
    if(borderNotConstant(remezdesc, i, 1, 1)){
      createModelphase3_border(remezdesc, 1, 1, path, i);
      System system(fullpath.c_str());
      DefaultSolver solver(system,1e-10, 1e-15);
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
  }
  return newPoints;
}

void createModelphase3_border(struct Multivariate_RemezIIParameters remezdesc, int borner1, int borner2, string path, int var){
  ofstream phase3(path + "models/phase3_model_border.txt");
  std::vector<double> A = remezdesc.poly.a;
  phase3 << setprecision(13);
  bool firstX = false;
  phase3 << "Variables\n";
  for(int i = 0; i<remezdesc.fdesc.domainX.size(); i++){
    if(i == var){
      phase3 << "  x"<< i+1 <<" in [" << remezdesc.fdesc.domainX[i][0] <<" , " << remezdesc.fdesc.domainX[i][1] <<"];\n";
    } else {
      if (firstX){
        phase3 << "  x"<< i+1 <<" in [" << remezdesc.fdesc.domainX[i][borner2] <<" , " << remezdesc.fdesc.domainX[i][borner2] <<"];\n";
      } else {
        phase3 << "  x"<< i+1 <<" in [" << remezdesc.fdesc.domainX[i][borner1] <<" , " << remezdesc.fdesc.domainX[i][borner1] <<"];\n";
        firstX = true;
      }
    }
  }
  phase3 << "\nConstraints\n";
  std::vector<std::vector<int>> p = remezdesc.poly.derivatedPhi(remezdesc.fdesc.domainX.size(), remezdesc.approximationDegree);
  
  for(int i = 0; i<remezdesc.fdesc.domainX.size(); i++){
    if(i==var){
      phase3 << remezdesc.fdesc.derivedFunctionStrings[i] << "-("; //derivated of f depending xi
      for(int j=0; j<p.size(); j++){//P'(x) depending of each x
        if(p[j][i]>0){
          phase3 << " + " << A[j] ;
          for(int k=0; k<p[j].size(); k++){
            if(k==i){
              if(p[j][k]>1){
                phase3 << "*((x" << k+1 << "*" << p[j][k] << ")^" << p[j][k]-1 << ")";
               }
             } else {
               if(k!=var){
                 if(p[j][k]!=0){
                   phase3 << "*((x" << k+1 << ")^" << p[j][k] << ")";
                  }
                }else{
                  if(p[j][k]>=1){
                    phase3 << "*((x" << k+1 << ")^" << p[j][k] << ")";
                  }
                }
              }
          }
        }
      }
      phase3 << ")=0;\n"; 
    }
  }
  phase3 << "end";
}

std::vector<std::vector<double>> phase3_side(struct Multivariate_RemezIIParameters remezdesc, string path, int turn){
  std::ofstream data;
  data.open("data/"+to_string(turn)+"points.txt",  std::ios_base::app);
  std::vector<std::vector<double>> newPoints = {};
  std::vector<double> point;
  string fullpath;
  for(int i = 0; i  < remezdesc.fdesc.domainX.size(); i++){
    if(cheking_side_not_constant(remezdesc, i*2 + 0, i)){
      createModelphase3_side(remezdesc.poly.a, remezdesc, remezdesc.fdesc.domainX[i][0], path, i);
    
      fullpath = path + "models/phase3_model_side.txt";
      System system(fullpath.c_str());
      DefaultSolver solver(system,1e-10, 1e-15);
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
      createModelphase3_side(remezdesc.poly.a, remezdesc, remezdesc.fdesc.domainX[i][1], path, i);
    
      fullpath = path + "models/phase3_model_side.txt";
      System system2(fullpath.c_str());
      DefaultSolver solver2(system2,1e-010, 1e-15);
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

void createModelphase3_side(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, double borner, string path, int var){
  ofstream phase3(path + "models/phase3_model_side.txt");
  phase3 << setprecision(13);
  
  phase3 << "Variables\n";
  for(int i = 0; i<remezdesc.fdesc.domainX.size(); i++){
    if(i != var){
      phase3 << "  x"<< i+1 <<" in [" << remezdesc.fdesc.domainX[i][0] <<" , " << remezdesc.fdesc.domainX[i][1] <<"];\n";
    } else {
      phase3 << "  x"<< i+1 <<" in [" << borner <<" , " << borner <<"];\n";
    }
  }
  phase3 << "\nConstraints\n";
  std::vector<std::vector<int>> p = remezdesc.poly.derivatedPhi(remezdesc.fdesc.domainX.size(), remezdesc.approximationDegree);
  
  
  //int i = 1 - var;
  for(int i = 0; i<remezdesc.fdesc.domainX.size(); i++){
    if(i!=var){
      phase3 << remezdesc.fdesc.derivedFunctionStrings[i] << "-("; //derivated of f depending xi
      for(int j=0; j<p.size(); j++){//P'(x) depending of each x
        if(p[j][i]>0){
          phase3 << " + " << A[j] ;
          for(int k=0; k<p[j].size(); k++){
            if(k==i){
              if(p[j][k]>1){
                phase3 << "*((x" << k+1 << "*" << p[j][k] << ")^" << p[j][k]-1 << ")";
               }
             } else {
               if(k==var){
                 if(p[j][k]!=0){
                   phase3 << "*((" << borner << ")^" << p[j][k] << ")";
                  }
                }else{
                  if(p[j][k]>=1){
                    phase3 << "*((x" << k+1 << ")^" << p[j][k] << ")";
                  }
                }
              }
          }
        }
      }
      phase3 << ")=0;\n"; 
    }
  }
  phase3 << "end";
}

std::vector<std::vector<double>> phase3_center(struct Multivariate_RemezIIParameters remezdesc, string path, int turn){
  std::ofstream data;
  data.open("data/"+to_string(turn)+"points.txt",std::ofstream::trunc);
  std::vector<std::vector<double>> newPoints = {};
  std::vector<double> point= {};
  createModelphase3(remezdesc.poly.a, remezdesc, path);
  string fullpath = path + "models/phase3_model.txt";
  System system(fullpath.c_str());
  DefaultSolver solver(system,1e-10, 1e-15);
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

void createModelphase3(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path){
  ofstream phase3(path + "models/phase3_model.txt");
  phase3 << setprecision(13);
  phase3 << "Variables\n";
  for(int i = 0; i<remezdesc.fdesc.domainX.size(); i++){
    phase3 << "  x"<< i+1 <<" in [" << remezdesc.fdesc.domainX[i][0] <<" , " << remezdesc.fdesc.domainX[i][1] <<"];\n";
  }
  phase3 << "\nConstraints\n";
  std::vector<std::vector<int>> p = remezdesc.poly.derivatedPhi(remezdesc.fdesc.domainX.size(), remezdesc.approximationDegree);
  //1st derivative test
  for(int i = 0; i<remezdesc.fdesc.nbX; i=i+1){
    phase3 << remezdesc.fdesc.derivedFunctionStrings[i] << "-("; //derivated of f depending xi
    for(int j=0; j<p.size(); j++){//P'(x) depending of each x
      if(p[j][i]>0){
        phase3 << " + " << A[j] ;
        for(int k=0; k<p[j].size(); k++){
          if(k==i){
            if(p[j][k]>1){
              phase3 << "*(x" << k+1 << "*" << p[j][k] << ")^" << p[j][k]-1;
            }
          } else {
            if(p[j][k]>=1){
              phase3 << "*(x" << k+1 << ")^" << p[j][k];
            }
          }
        }
      }
    }
    phase3 << ")=0;\n";
  }
  
  phase3 << "end";
  phase3.close();

}

//for the best estimation of error phase 3, we search it separetly from the stationnary points
double phase3Error(std::vector<std::vector<double>> newPoints, std::vector<double> a, struct Multivariate_RemezIIParameters remezdesc, string path){
  double errorStep = 0;
  createModelphase3Error(a, remezdesc, path);
  string fullpath = path + "models/phase3_model_error.txt";
  System sys(fullpath.c_str());
  DefaultOptimizer optimizer(sys,1e-09, 1e-10);
  optimizer.optimize(sys.box);
  optimizer.report();
  errorStep = optimizer.get_loup();
  return -errorStep;
}

void createModelphase3Error(std::vector<double> A, struct Multivariate_RemezIIParameters remezdesc, string path){
  ofstream phase3(path + "models/phase3_model_error.txt");
  phase3 << setprecision(13);
  phase3 << "Variables\n";
  for(int i = 0; i<remezdesc.fdesc.domainX.size(); i++){
    phase3 << "  x"<< i+1 <<" in [" << remezdesc.fdesc.domainX[i][0] <<" , " << remezdesc.fdesc.domainX[i][1] <<"];\n";
  }
  phase3 << "\nMinimize\n";
  phase3 << "-abs(";
  phase3 << remezdesc.fdesc.functionString << "-(";
  std::vector<std::vector<int>> p = remezdesc.poly.derivatedPhi(remezdesc.fdesc.domainX.size(), remezdesc.approximationDegree);
  for(int i = 0; i < A.size(); i++){
    phase3 << "+(" << A[i] << ")" ;
    for(int j = 0; j < p[i].size(); j++){
      if(p[i][j] != 0){
        phase3 << "*(x" << j+1 << "^(" << p[i][j] << "))";
      }
    }
  }
  phase3 << "));";
  phase3.close(); 
}

//checks a border isn't constant
bool borderNotConstant(struct Multivariate_RemezIIParameters remezdesc, int var, int borner1, int borner2){
  bool doSide = false;
  int firstToDerive = -1;
  int secondToDerive = -1;
  for(int i = 0; i < remezdesc.fdesc.nbX; i++){
    if(i!=var){
      if(firstToDerive != -1){
        secondToDerive = i;
      } else {
        firstToDerive = i;
      }
    }
  }
  bool needToCheck = !(remezdesc.fdesc.notConstantDerived[firstToDerive*2 + borner1]);
  needToCheck = needToCheck || !(remezdesc.fdesc.notConstantDerived[secondToDerive*2 + borner2]);
  if(needToCheck){
    std::vector<std::vector<int>> p = remezdesc.poly.derivatedPhi(remezdesc.fdesc.domainX.size(), remezdesc.approximationDegree);
    for(int i = 0; i < p.size(); i++){
      if(p[i][var] > 1 && remezdesc.poly.a[i] != 0){
        if(!((p[i][firstToDerive] > 0 && remezdesc.fdesc.domainX[firstToDerive][borner1] == 0) || (p[i][secondToDerive] > 0 && remezdesc.fdesc.domainX[secondToDerive][borner2] == 0))){
          doSide = true; 
        }
      }
    }
    return doSide;
  } else {
    return true;
  }
}

//checks a side isn't constant
bool cheking_side_not_constant(struct Multivariate_RemezIIParameters remezdesc, int side, int var){
  bool doSide = remezdesc.fdesc.notConstantDerived[side];
  
  if(!remezdesc.fdesc.notConstantDerived[side]){
    doSide = false;
    std::vector<std::vector<int>> p = remezdesc.poly.derivatedPhi(remezdesc.fdesc.domainX.size(), remezdesc.approximationDegree);
    if(remezdesc.fdesc.nbX == 2){
      for(int i = 0; i < p.size(); i++){
        if(p[i][(var+1)%2] > 1 && remezdesc.poly.a[i] != 0){
          if(!(p[i][var] > 0 && remezdesc.fdesc.domainX[var][side%remezdesc.fdesc.nbX] == 0)){
            doSide = true; 
          }
        }
      }
    } else {
      bool doOneSide = false;
      bool isVarOneUsed = false;
      int whoIsVarOne = -1;
      bool isVarTwoUsed = false;
      doSide = true;
      int nbX = remezdesc.fdesc.nbX;
      for(int i = 0; i < p.size(); i++){
        for(int j = 0; j < p[i].size(); j++){
          if(j != var){
            if(p[i][j] > 1 && remezdesc.poly.a[i] != 0){
            if(!(p[i][var] > 0 && remezdesc.fdesc.domainX[var][(var*2)%nbX + side%nbX] == 0)){
              doOneSide = true; 
              if(j != whoIsVarOne){
                if(whoIsVarOne != -1){
                  isVarTwoUsed = true;
                }else{
                  whoIsVarOne = j;
                  isVarOneUsed = true;
                }
              }
            }
            }
          }
        }
        doSide = doSide || doOneSide;
      }
      doSide = doSide && isVarOneUsed && isVarTwoUsed;
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
  for(int i = 0; i<remezdesc.fdesc.domainX.size(); i++){
      summary << "[" << remezdesc.fdesc.domainX[i][0] << "," << remezdesc.fdesc.domainX[i][1] << "]";
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
    summary << "," << setprecision(19) << A[j] ;
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

void write_data_for_graphs(std::vector<double> A, std::vector<std::vector<double>> domainX, int turn){
  std::ofstream data;
  data << setprecision(13) ;
  data.open("data/"+ to_string(turn) +"data.txt",std::ofstream::trunc);
  for(int i = 0; i<A.size(); i++){
    data << A[i] << " ";
  }
  //for ease, we put 0 until having the same number of "case" for the python tab
  for(int i = A.size()-1 ; i<=domainX.size(); i++){
    data << "0 ";
  }
  
  data << "\n";
  for(int i = 0; i<domainX.size(); i++){
    data << domainX[i][0] << " " << domainX[i][1] << " ";
  }
  for(int i = domainX.size()+2; i<A.size(); i++){
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
  Vec x0 = linspace(remezdesc.fdesc.domainX[0][0], remezdesc.fdesc.domainX[0][1], 100);
  Vec x1 = linspace(remezdesc.fdesc.domainX[1][0], remezdesc.fdesc.domainX[1][1], 100);
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
  canvas.save(path + "errorEndphase2.pdf");
}


void grapheConvergenceComparisons(struct Multivariate_RemezIIResult remezResults){
  string graphName = "../resultsOfSteps/multi/GrapheConvergenceAll.pdf";
  string output = "";
  std::vector<double> distanceTemp = {};
  std::vector<std::vector<double>> distance = {};
  double dist = 0;
  std::vector<double> conv = {};
    for(int i = 0; i<remezResults.errorphase2.size(); i++){
      dist = sqrt(pow((remezResults.errorphase3[i]-remezResults.errorphase2[i]),2));
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
  cout << "\n" << *min_element(conv.begin(), conv.end());
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
