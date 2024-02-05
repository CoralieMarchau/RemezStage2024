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

int RemezUni(double approx, int sizeD, int nbTurns, double up, double down);
bool tooClose(std::vector<double> D, int sizeD, double newPoint, double approx);
std::vector<double> step0Univarie(int sizeD, double up, double low);
void step1(std::vector<double> D, int sizeD);
void createMPSFile(std::vector<double> D, int sizeD);
std::vector<string> parserResultStep1();
void createModelStep2(std::vector<string> P, double up, double down);
void step2(std::vector<string> P, double up, double down);
std::vector<double> parserResultStep2();




int main(int argc, char** argv) { 
  double approx = 0.000000005;
  int sizeD = 5;
  int nbTurns = 100;
  double up = 1.0;
  double down = -1.0;
  RemezUni(approx, sizeD, nbTurns, up, down);
  //RemezUniPlus(approx, sizeD, nbTurns, up, down, function, derivedFunction);
}


int RemezUni(double approx, int sizeD, int nbTurns, double up, double down){
  ofstream summary("../resultsOfSteps/Summary_UnivariéSingulier");
  summary << "Welcome to the summary of Remez 1 for the singular univariate problem. \n";
  summary << "The parameters of Remez have been initialized as : \n";
  double error, newPoint;
  std::vector<double> errorStep = {};
  double errorStep1 = 0.0;
  int i = 1;
  bool itIsTheEnd = false;
  std::vector<double> D = step0Univarie(sizeD, up, down);
  
  summary << "-Size of the first discretized set : " << sizeD << "\n";
  summary << "-Approximation of x : " << approx << "\n";
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
    step1(D, sizeD);
    std::vector<string> P = parserResultStep1();
    summary << "[" << P[0] ;
    summary << "] \n" ;

      summary << "-a"  << "= " << P[0] << "\n";
    istringstream  istr(P[1]);
    istr >> errorStep1;
    //cout << "\n" << P[nbA+1] << "\n";
    P.pop_back();
    step2(P, up, down);
    std::vector<double> step2Result = parserResultStep2();
    error = step2Result[0];
    newPoint = step2Result[1];
    errorStep.push_back(errorStep1);
    summary << "\nError : " << error << "\n" ;
    summary << "Error for x in D" << i << ": " << errorStep1 << "\n" ;
    if (tooClose(D, sizeD, newPoint, approx) || i >= nbTurns){
      itIsTheEnd = true ;
      summary << "\n Last error is : " << errorStep[i-1] << " \n";
      summary << "\n END \n";
    } else {
        D.push_back(newPoint);
        sizeD = sizeD + 1;
        summary << "Added x is : " << newPoint << "\n";
    }
    i = i + 1;
  }
  summary << "\n Remez finished after " << i-1 << " turns. \n";
  return i;
}

//Tools : 

bool tooClose(std::vector<double> D, int sizeD, double newPoint, double approx){
  bool isClose = false;
  for(int i=0 ; i<sizeD; i++){
    if(D.at(i)>=newPoint-approx && D.at(i)<=newPoint+approx){ 
      isClose = true;
      cout << D.at(i) << " " << newPoint << "\n";
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

//step 1 using glpk
void step1(std::vector<double> D, int sizeD){
  createMPSFile(D, sizeD);
  glp_prob *mip;
  glp_tran *tran;
  int ret;
  mip = glp_create_prob();
  tran = glp_mpl_alloc_wksp();
  ret = glp_mpl_read_model(tran, "step1_uniSing.mod", 1);
  ret = glp_mpl_generate(tran, NULL);
  glp_mpl_build_prob(tran, mip);
  glp_simplex(mip, NULL);
  glp_intopt(mip, NULL);
  glp_print_sol(mip, "../resultsOfSteps/step1Result_UniSing.txt");
  glp_mpl_free_wksp(tran);
  glp_delete_prob(mip);   
}

//Creates the model to solve
void createMPSFile(std::vector<double> D, int sizeD){
  ofstream step1("step1_uniSing.mod");
  double x = 0.0;
  step1 << "# Parameters \n";
  step1 << "\n";
  
  step1 << "# Variable \n";
  step1 << "var a" << ";\n";
  step1 << "var y; \n";
  
  step1 << "# Objectif \n";  
  step1 << "minimize z : y; \n";
  step1 << "\n";
  
  step1 << "# Constraint \n";  
  for(int i = 0; i<sizeD; i++){
    step1 << "d" << i << " : ";
    step1 << "2*a*" << D[i] << "-" << pow(D[i], 2) << "<=y;\n";
  }
  /*  
  for(int i = 0; i<sizeD; i++){
    step1 << "d" << i+sizeD << " : ";
    step1 << "2*a*" << D[i] << "-" << pow(D[i], 2) << ">= -y;\n";
  }*/
  
  step1 << "\n solve ; \n";  

    step1 << "display a" << "; \n";
  step1.close();
}


//Parser recuperating the string of the polynomial 
//(no need to change them in double since it is to write them in step 2)
std::vector<string> parserResultStep1(){
  ifstream resultStep1("../resultsOfSteps/step1Result_UniSing.txt");
  std::vector<string> P = {};
  if(resultStep1) {
    bool varReached = false;
    string str;
    resultStep1 >> str;
    //To reach the first a
    while(str != "a"){
      resultStep1 >> str;
    }
      resultStep1 >> str;
      resultStep1 >> str;
      P.push_back(str);
      resultStep1 >> str;
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
  else {
    cout << "ERROR FILE NOT OPENING" << endl;
  }
  return P;
}

//STEP 2 :

void createModelStep2(std::vector<string> P, double up, double down){
  ofstream step2("step2_model_uniSing.txt");
  step2 << "Variables\n";
  step2 << "  x in [" << down <<" , " << up <<"];\n";
  step2 << "\nMinimize\n";
  
  step2 << "2*x*" << P[0] << "-x^2;\n";

  step2 << "\nConstraints\n";
  step2 << "end";
  step2.close();
}

void step2(std::vector<string> P, double up, double down){
  createModelStep2(P, up, down);
  ofstream step2R("../resultsOfSteps/step2_result_uniSing.txt");

  System sys("step2_model_uniSing.txt");
  DefaultOptimizer optimizer(sys,1e-07, 1e-10);
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

std::vector<double> parserResultStep2(){
  ifstream resultStep2("../resultsOfSteps/step2_result_uniSing.txt");  //Ouverture d'un fichier en lecture
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


