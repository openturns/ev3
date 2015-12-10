
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <iostream>

#include "expression.h"
#include "parser.h"
#include "muParser.h"

class Description : public std::vector<std::string>
{
public:

  explicit Description(const unsigned long Size = 0, const std::string & expr = "")
    : std::vector<std::string>(Size)
  { 
    for (unsigned long i = 0; i < Size; ++i)
    {
      (*this)[i] = expr;
    }
  }

  void add(const std::string & elt) {push_back(elt);}
    std::string toString() const
  {
    std::stringstream oss;
    oss << "[";
    for (unsigned long i = 0; i < size(); ++i)
    {
      if (i>0) oss << ",";
      oss << (*this)[i];
    }
    oss << "]";
    return oss.str();
  }
};


class NumericalPoint : public std::vector<double>
{
public:
  explicit NumericalPoint(const unsigned long Size = 0, const double expr = 0.)
    : std::vector<double>(Size)
  {
    for (unsigned long i = 0; i < Size; ++i)
    {
      (*this)[i] = expr;
    }
  }

  std::string toString() const
  {
    std::stringstream oss;
    oss << "[";
    for (unsigned long i = 0; i < size(); ++i)
    {
      if (i>0) oss << ",";
      oss << (*this)[i];
    }
    oss << "]";

    return oss.str();
  }
};


class AnalyticalParser : public mu::Parser
{

public:
  AnalyticalParser() : mu::Parser()
{
  DefineFun(_T("cotan"), Cotan); // modified
  DefineFun(_T("acotan"), ACotan); // modified
  DefineFun(_T("asinh"), ASinh); // modified
  DefineFun(_T("acosh"), ACosh); // modified
  DefineFun(_T("atanh"), ATanh); // modified
  DefineFun(_T("log"), Ln); // modified: assigned to log10 by default
  DefineFun(_T("log2"), Log2); // modified
  DefineFun(_T("lngamma"), LnGamma); // added
  DefineFun(_T("gamma"), Gamma); // added
  DefineFun(_T("erf"), Erf); // added
  DefineFun(_T("erfc"), Erfc); // added
  DefineFun(_T("abs"), Abs); // modified
  DefineFun(_T("cbrt"), Cbrt); // added
  DefineFun(_T("besselJ0"), J0); // added
  DefineFun(_T("besselJ1"), J1); // added
  DefineFun(_T("besselY0"), Y0); // added
  DefineFun(_T("besselY1"), Y1); // added
  DefineFun(_T("rint"), Rint); // modified
}

protected:

static double Cotan(double v)
{
  return 1.0 / tan(v);
}
static double ACotan(double v)
{
  if (v < 0.0) return -M_PI_2 - atan(v);
  return M_PI_2 - atan(v);
}
static double ASinh(double v)
{
  return asinh(v);
}
static double ACosh(double v)
{
  return acosh(v);
}
static double ATanh(double v)
{
  return atanh(v);
}
static double Ln(double v)
{
  return log(v);
}
static double Log2(double v)
{
  return log2(v);
}
static double LnGamma(double v)
{
  return lgamma(v);
}
static double Gamma(double v)
{
  return tgamma(v);
}
static double Erf(double v)
{
  return erf(v);
}
static double Erfc(double v)
{
  return erfc(v);
}
static double Abs(double v)
{
  return std::abs(v);
}
static double Cbrt(double v)
{
  return cbrt(v);
}
static double J0(double v)
{
  return j0(v);
}
static double J1(double v)
{
  return j1(v);
}
static double Y0(double v)
{
  return y0(v);
}
static double Y1(double v)
{
  return y1(v);
}
static double Rint(double v)
{
  return round(v);
}


};

class Function {
public:

  Function(const Description & inputVariables, const Description & evaluation)
  : inputVariables_(inputVariables), evaluation_(evaluation), gradient_(inputVariables.size())
  {
    unsigned long inputSize = inputVariables_.size();

    int nerr = 0;
    Ev3::Expression ev3Expression;
    Ev3::ExpressionParser ev3Parser;
//     try {
        for (unsigned long i = 0; i < inputSize; ++i)
          ev3Parser.SetVariableID(inputVariables_[i], i);
        ev3Expression = ev3Parser.Parse(evaluation_[0].c_str(), nerr);

        if(nerr!=0) throw std::exception();
        for (unsigned long i = 0; i < inputVariables_.size(); ++i) {
          Ev3::Expression derivative = Ev3::Diff(ev3Expression, i);
          std::stringstream oss;
          oss << std::setprecision(12) << derivative;
          gradient_[i] = oss.str();
        }
//       }
//       catch (Ev3::ErrBase & exc) {
//         throw std::exception();
//       }
      
    
  }

  
  NumericalPoint eval(const NumericalPoint & inP)
  {
    AnalyticalParser parser;
    unsigned long inputSize = inputVariables_.size();
    NumericalPoint buf(inP);
    for (unsigned long i = 0; i < inputSize; ++i) {
       parser.DefineVar(inputVariables_[i].c_str(), &buf[i]);
    }
    parser.SetExpr(evaluation_[0].c_str());
    return NumericalPoint(1, parser.Eval());
  }
  
  NumericalPoint grad(const NumericalPoint & inP)
  {
    unsigned long inputSize = inputVariables_.size();
    NumericalPoint buf(inP);

    NumericalPoint result(inputSize);
    for (unsigned long i = 0; i < inputSize; ++i) {
      AnalyticalParser parser;
      for (unsigned long j = 0; j < inputSize; ++j) {
        parser.DefineVar(inputVariables_[j].c_str(), &buf[j]);
      }

      parser.SetExpr(gradient_[i].c_str());
      result[i] =  parser.Eval();
      
    }
    return result;
  }
  
  NumericalPoint grad_fd(const NumericalPoint & inP)
  {
    unsigned long inputSize = inputVariables_.size();
    double eps = 1e-5;
    
    NumericalPoint result(inputSize);
    for (unsigned long i = 0; i < inputSize; ++i) {
      NumericalPoint x1(inP);
      x1[i] += eps;
      NumericalPoint x2(inP);
      x2[i] -= eps;
      double df = (eval(x1)[0]-eval(x2)[0])/(2.*eps);
      result[i] = df;
    }
    return result;
  }
  
  std::string toString() const {
    std::stringstream oss;
    oss << "f(";
    for (unsigned long i = 0; i < inputVariables_.size(); ++i) {
      if (i>0) oss << ",";
      oss << inputVariables_[i];
    }
      
    oss <<")=" << evaluation_[0]<< std::endl;
    for (unsigned long i = 0; i < inputVariables_.size(); ++i) {
      if (i>0) oss << std::endl;
      oss << "  df/d"<<inputVariables_[i]<<"="<<gradient_[i];
    }
    return oss.str();
  }
  
private:  
  Description inputVariables_;
  Description evaluation_;
  Description gradient_;
  

};

static std::string randFunc() {
  Description elementaryFunctions;
  elementaryFunctions.add("sin");
    elementaryFunctions.add("cos");
    elementaryFunctions.add("tan");
    elementaryFunctions.add("asin");
    elementaryFunctions.add("acos");
    elementaryFunctions.add("atan");
    elementaryFunctions.add("sinh");
    elementaryFunctions.add("cosh");
    elementaryFunctions.add("tanh");
    elementaryFunctions.add("asinh");
    elementaryFunctions.add("acosh");
    elementaryFunctions.add("atanh");
    elementaryFunctions.add("log2");
    elementaryFunctions.add("log10");
    elementaryFunctions.add("log");
    elementaryFunctions.add("ln");
//     elementaryFunctions.add("lngamma");
//     elementaryFunctions.add("gamma");
    elementaryFunctions.add("exp");
    elementaryFunctions.add("erf");
    elementaryFunctions.add("erfc");
    elementaryFunctions.add("sqrt");
    elementaryFunctions.add("cbrt");
//     elementaryFunctions.add("besselJ0");
//     elementaryFunctions.add("besselJ1");
//     elementaryFunctions.add("besselY0");
//     elementaryFunctions.add("besselY1");
//     elementaryFunctions.add("sign");
//     elementaryFunctions.add("rint");
    elementaryFunctions.add("abs");
    
  unsigned long index = rand() % elementaryFunctions.size();
  return elementaryFunctions[index];
}


static std::string randVar(int n)
{
  std::stringstream oss;
  oss << "x" << 1+(rand()%n);
  return oss.str();
}

static std::string randOp2()
{
  Description ops;
  ops.add("+");
  ops.add("-");
  ops.add("*");
  ops.add("/");
  ops.add("^");
  unsigned long index = rand() % ops.size();

  return ops[index];
}

static std::string randCoef()
{
  std::stringstream oss;
  oss << (rand()%12) - 5;
  return oss.str();
}

static std::string randExp(int n, int dim)
{
  if (n==0)
      return randVar(dim);
  else {
    int n1 = rand() % n;
    std::string exp1 = randExp(n1, dim);
    int t = rand() % 10;
    if (t<3)
        return "("+randCoef()+"*"+exp1+")";
    else if (t<7)
        return randFunc()+"("+exp1+")";
    else {
        int n2 = rand() %n;
        std::string exp2 = randExp(n2, dim);
        return "("+exp1+")"+randOp2()+"("+exp2+")";
    }
  }
}




int main()
{
  
  unsigned long dimension = 2;
  Description inputVars(dimension);
  NumericalPoint x(dimension, 0.4);
  
  for (unsigned long i = 0; i < dimension; ++i) {
    std::stringstream oss;
    oss << "x"<<i+1;
    inputVars[i] = oss.str();
  }
  
  for (unsigned long j = 0; j < 10000; ++ j) 
  {
    try {
      Description formulas;
      formulas.add("log2(((x1)+(x2))*(erf(x2)))");
      formulas.add("(x2)/(((x1)*(x1))+((-2*x1)))");
      formulas.add("sqrt(acosh((x1)-(x1)))");
      formulas.add("1./(1.-5*x1-1.)");
      formulas.add("x1*sinh(-3*x2/x2)");
      formulas.add("x1-tan(x1)+log(x1)");
      formulas.add("erf(x2-(x1-x2))");
      formulas.add("2*log(x1)*log(x1)");
      formulas.add("erf((x2)-((x2)+(x1)))");
      formulas.add("1/(((x2)-(x1))-(x2))");
      formulas.add("log(x1)+x1-log(x1)");
      formulas.add("x1*x2*log(x1*x1)");
      formulas.add("log2(x1)/abs(log(x1))");
      formulas.add("(4*x1)^(1/2)");
      formulas.add("(x1)^(1/log10(x1))");
      formulas.add("cos(x1)*sin(x1)");
      formulas.add("x1^(0.5*(x1*x2)/(x2))");
      formulas.add("exp(-x1^2/2)");
      formulas.add("((-1*tanh((x1)*(x1))))/((-1*x1))");
      formulas.add("(x2*x1)*(x1+x2-x2)");
      formulas.add("(4*(3*x1))/log10(exp(x1))");
      formulas.add("-(x1-1)^2");
      formulas.add("exp(-(x1-1)^2)");
      formulas.add("-x1+x1--3*x1+-7*x1");
       
      if (j < formulas.size())
        formulas = Description(1, formulas[j]);
      else
        formulas = Description(1, randExp(4,2));

      std::cout << "formula=" << formulas.toString() << std::endl;
      Function function(inputVars, formulas);
      std::cout << function.toString()<<std::endl;
      double df = function.grad(x)[0];
      if(formulas[0]== std::string("tanh(((x2)-(x2))^((x1)-(x2)))"))
        continue;
       if(formulas[0]== std::string("(((x1)-(x1))^((x1)-(x2)))*(x1)"))
        continue;
            if(std::isnan(df) || std::isinf(df))
        continue;

      double df2 = function.grad_fd(x)[0];
      if( std::isnan(df2)|| std::isinf(df2) || (df2==0))
        continue;
      double err_g = 0.;
      if (std::abs(df)>1e5||std::abs(df2)<1e-10)
        continue;
      if (std::abs(df) > 1e-5)
        err_g = std::abs(df2/df-1.);
      else
        err_g = std::abs(df - df2);

      if (err_g > 1e-2) {
        std::cout << "XXXXXXXXX df="<<df<<" df2="<<df2<<" err="<<err_g<<std::endl;
        throw std::exception();
      }
    } catch (mu::ParserError & ex)
    {
      continue;
    }    
    catch (Ev3::ErrBase & ex)
    {
      continue;
    }
  }
  return 0;
}

