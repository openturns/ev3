
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <iostream>

#include "expression.h"
#include "parser.h"
#include "exprtk.hpp"

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


class Point : public std::vector<double>
{
public:
  explicit Point(const unsigned long Size = 0, const double expr = 0.)
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


class Function {
public:

  Function(const Description & inputVariables, const Description & evaluation)
  : inputVariables_(inputVariables), evaluation_(evaluation), gradient_(inputVariables.size())
  {
    unsigned long inputSize = inputVariables_.size();

    int nerr = 0;
    Ev3::Expression ev3Expression;
    Ev3::ExpressionParser ev3Parser;

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
  }

  static double Function_sign(double v)
  {
    if (v > 0.0) return 1.0;
    else if (v < 0.0) return -1.0;
    else return 0.0;
  }
  
  static exprtk::symbol_table<double> get_symbol_table()
  {
    exprtk::symbol_table<double> symbol_table;
    symbol_table.add_constant("e_", 2.71828182845904523536028747135266249775724709369996);
    symbol_table.add_constant("pi_", 3.14159265358979323846264338327950288419716939937510);
    symbol_table.add_function("sign", Function_sign);
    symbol_table.add_function("ln", std::log);
    symbol_table.add_function("cbrt", cbrt);
    return symbol_table;
  }
  
  Point eval(const Point & inP)
  {
    exprtk::symbol_table<double> symbol_table = get_symbol_table();
    unsigned long inputSize = inputVariables_.size();
    Point buf(inputSize + 1);
    for (unsigned long i = 0; i < inputSize; ++i)
    {
      buf[i] = inP[i];
      if (!symbol_table.add_variable(inputVariables_[i], buf[i]))
        throw std::runtime_error(std::string("add_variable error"));
    }
    if (!symbol_table.add_variable("y", buf[inputSize]))
      throw std::runtime_error(std::string("add_variable error"));
    
    exprtk::parser<double> parser;
    exprtk::expression<double> expression;
    expression.register_symbol_table(symbol_table);
    if (!parser.compile(evaluation_[0], expression))
      throw std::runtime_error(std::string("eval compile error: ") + parser.get_error(0).diagnostic);
    expression.value();
    return Point(1, buf[inputSize]);
  }
  
  Point grad(const Point & inP)
  {
    exprtk::symbol_table<double> symbol_table = get_symbol_table();
    unsigned long inputSize = inputVariables_.size();
    Point buf(inputSize + 1);
    for (unsigned long i = 0; i < inputSize; ++i)
    {
      buf[i] = inP[i];
      if (!symbol_table.add_variable(inputVariables_[i], buf[i]))
        throw std::runtime_error(std::string("add_variable error"));
    }
    if (!symbol_table.add_variable("y", buf[inputSize]))
      throw std::runtime_error(std::string("add_variable error"));
    
    exprtk::parser<double> parser;
    exprtk::expression<double> expression;
    expression.register_symbol_table(symbol_table);
    Point grad(inputSize);
    for (unsigned long i = 0; i < inputSize; ++i)
    {
      if (!parser.compile(gradient_[i], expression))
        throw std::runtime_error(std::string("grad compile error: ") + parser.get_error(0).diagnostic);
      expression.value();
      grad[i] = buf[inputSize];
    }
    return grad;
  }

  Point grad_fd(const Point & inP)
  {
    unsigned long inputSize = inputVariables_.size();
    double eps = 1e-5;
    
    Point result(inputSize);
    for (unsigned long i = 0; i < inputSize; ++i) {
      Point x1(inP);
      x1[i] += eps;
      Point x2(inP);
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
  Point x(dimension, 0.4);
  
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
      formulas.add("(3*x1)^(-1)");
//       formulas.add("(-2*x1^2)^(-4*(x2)/(x2))");

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
        std::cout << "Wrong gradient! df="<<df<<" df2="<<df2<<" err="<<err_g<<std::endl;
        throw std::logic_error("wrong gradient");
      }
    }
    catch (std::exception & ex)
    {
      std::cout << ex.what() << std::endl;
      continue;
    }    
    catch (Ev3::ErrBase & ex)
    {
      continue;
    }
  }
  return 0;
}

