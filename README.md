[![Build Status](https://travis-ci.org/openturns/ev3.svg?branch=master)](https://travis-ci.org/openturns/ev3)

Ev3
===

Ev3 is a C++ library to compute symbolic derivatives written by [Leo Liberti] (http://www.lix.polytechnique.fr/~liberti/academic.html).

It is an alternative to [libmatheval](http://www.gnu.org/software/libmatheval/) but under the more permissive [CPL](http://en.wikipedia.org/wiki/Common_Public_License) license.

The library was originaly available [here](http://www.lix.polytechnique.fr/~liberti/Ev3-1.0.tar.gz), and used also [here] (https://projects.coin-or.org/ROSE/browser/rose/Ev3) with some modifications. This is the debugged version used in the [OpenTURNS](http://www.openturns.org) software.

Here a snippet which define two variables and derivates an expression:
```
#include <iostream>
#include "expression.h"
#include "parser.h"

int main()
{  
  Ev3::ExpressionParser parser;
  parser.SetVariableID("x1", 0);
  parser.SetVariableID("x2", 1);
  int nerr = 0;
  Ev3::Expression expr = parser.Parse("x1*sin(x2)", nerr);
  Ev3::Expression derivative = Ev3::Diff(expr, 1);
  std::cout << derivative->ToString() << std::endl;
}
```

