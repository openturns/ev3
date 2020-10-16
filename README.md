[![Build Status](https://travis-ci.org/openturns/ev3.svg?branch=master)](https://travis-ci.org/openturns/ev3)

Ev3
===

Ev3 is a C++ [library for symbolic computation](http://www.lix.polytechnique.fr/Labo/Leo.Liberti/Ev3.pdf) written by [Leo Liberti](http://www.lix.polytechnique.fr/~liberti/academic.html).

It is an alternative to [libmatheval](http://www.gnu.org/software/libmatheval/).

It was [originally published](http://www.lix.polytechnique.fr/~liberti/Ev3-1.0.tar.gz) by Leo Liberti under the [CPL](http://en.wikipedia.org/wiki/Common_Public_License) license.

The Computational Infrastructure for Operations Research project ([COIN-OR](https://www.coin-or.org/)) uses a modified version in [ROSE](https://github.com/coin-or/ROSE/).

Leo Liberti has since [republished the original Ev3](http://www.lix.polytechnique.fr/~liberti/Ev3-1.0.zip) under the LGPL license (see `COPYING` and `COPYING.LESSER`).

This repository stores a debugged version used in the [OpenTURNS](http://www.openturns.org) software also published under the LGPL.

Here is a snippet which defines two variables and differentiates an expression:
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

