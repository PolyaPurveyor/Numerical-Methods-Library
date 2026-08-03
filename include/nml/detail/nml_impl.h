/*
   Egwuchukwu Kalu
   Numerical Methods Library
   Template implementations
*/
#ifndef NML_DETAIL_NML_IMPL_H
#define NML_DETAIL_NML_IMPL_H


namespace NML
{
 namespace utils
 {
  //Generates an interval from two floating point numbers
  template <typename T>
  std::vector<T> interval_gen(T x0,T xn,std::size_t n)
  {
   //Divisor must be non-zero
   if (n == 0)
    throw std::invalid_argument("n must be greater than 0");

   //Swaps the positions if x0 > xn (the beginning of the interval is greater than the end)
   if (x0 > xn)
   {
    std::cout << "Error: the beginning of the interval is greater than the end of the "
     "interval. Swapping positions." << std::endl;
    std::swap(x0,xn);
   }

   //Push equally spaced increments from x0 to xn onto the vector
   T dx = (xn-x0)/static_cast<T>(n);
   std::vector<T> interval;
   interval.reserve(n+1);

   for (std::size_t i = 0; i <= n; ++i)
    interval.push_back(x0+static_cast<T>(i)*dx);

   interval.back() = xn;
   return interval;
  }
 }

 namespace functional
 {
  /*Some equation solving methods require a derivative,
   * so the finite difference is first*/
  namespace diff
  {
   /*Forward Difference Formula */
   template <typename T,typename Function>
   T fdiff(Function f,T x,T dx)
   {
    /*Make sure divisor is non-zero
     *it can be negative however to my surprise*/
    detail::validate_nonzero(static_cast<long double>(dx),"dx");
    return (f(x+dx)-f(x))/dx;
   }

   /*Backward Difference Formula, which is just forward differentiation with a negative dx lol*/
   template <typename T,typename Function>
   T bdiff(Function f,T x,T dx)
   {
    return fdiff<T>(f,x,-dx);
   }

   /*Central Difference Formula*/
   template <typename T,typename Function>
   T cdiff(Function f,T x,T dx)
   {
    detail::validate_nonzero(static_cast<long double>(dx),"dx");
    return (f(x+dx)-f(x-dx))/(static_cast<T>(2)*dx);
   }
  }

  namespace algebraic
  {
   /*Horner's Method for Synthetic Division*/
   template <typename T>
   std::pair<T,T> horner(const std::vector<T>& P,T x0)
   {
    if (P.empty())
     throw std::invalid_argument("Polynomial coefficient vector cannot be empty");

    T Px = P[0],Pprime = T{};
    if (P.size() > 1)
     Pprime = P[0];

    for(std::size_t i = 1; i < P.size(); ++i)
    {
     Px = x0*Px+P[i];
     if (i < P.size()-1)
      Pprime = x0*Pprime+Px;
    }

    return std::pair<T,T>(Px,Pprime);
   }

   /*Quadratic Synthetic-Division*/
   template <typename T>
   std::pair<T,T> qsd(const std::vector<T>& R,T B1,T B0)
   {
    if (R.size() < 2)
     throw std::invalid_argument("Quadratic synthetic division requires at least two coefficients");

    T Z = R[0],C = R[1],ZZ;
    for(std::size_t i = 1; i < R.size(); ++i)
    {
     ZZ = Z;
     Z = -B1*ZZ+C;
     if (i != R.size()-1)
      C = -B0*ZZ+R[i+1];
    }

    return std::pair<T,T>(Z,C);
   }
  }

  //For all these methods, the functions are assumed to be continuous
  //(or differentiable if necessary as with Newton-Raphson) or otherwise not pathological
  namespace eqsoln
  {
   //Bisection Root-Finding Method: Finds a solution to f(x) = 0 given the continuous function f on the interval [a,b] where f(a) and f(b)
   //have opposite signs through a binary search where x = a+(b-a)/2
   template <typename T,typename Function>
   T bisection(Function f,T a,T b,T tolerance,std::size_t max_iterations)
   {
    detail::validate_positive(static_cast<long double>(tolerance),"tolerance");
    detail::validate_iterations(max_iterations);

    if (a > b)
     std::swap(a,b);

    T fa = f(a),fb = f(b);
    if (fa == 0)
     return a;
    if (fb == 0)
     return b;
    if ((fa > 0) == (fb > 0))
     throw std::invalid_argument("Bisection requires endpoints with opposite signs");

    //While within the range of allowed iterations
    for(std::size_t i = 0; i < max_iterations; ++i)
    {
     /*Set x to be the midpoint of the interval [a,b]*/
     T x = a+(b-a)/static_cast<T>(2);
     T fx = f(x);

     /*If x is a root or x is within the tolerance range of the root, set as solved
      *If not, if f(a) and f(x) are the same sign, set a to be x. Set b to be x otherwise.*/
     if(fx == 0 || std::abs((b-a)/static_cast<T>(2)) <= tolerance)
      return x;

     if ((fa > 0) == (fx > 0))
     {
      a = x;
      fa = fx;
     }
     else
      b = x;
    }

    detail::throw_no_convergence("Bisection method");
   }

   //Fixed Point Iteration: To find a solution to x = f(x), give an initial solution x0 and iteratively set x = f(x0)
   //and set x0 = x if not within the proper tolerance
   template <typename T,typename Function>
   T fixed_point(Function f,T x0,T tolerance,std::size_t max_iterations)
   {
    detail::validate_positive(static_cast<long double>(tolerance),"tolerance");
    detail::validate_iterations(max_iterations);

    for(std::size_t i = 0; i < max_iterations; ++i)
    {
     //Set x to be f(x0). If within the tolerated range for a solution, set as solved.
     //Otherwise, set x0 equal to x (i.e. x0 = f(x0)), and iterate the action of f on
     //x until either a solution is found or the number of iterations exceeds the maximum
     T x = f(x0);
     if (std::abs(x-x0) <= tolerance)
      return x;
     x0 = x;
    }

    detail::throw_no_convergence("Fixed-point iteration");
   }

   //Newton-Raphson Method: f must be differentiable (twice for fast convergence)
   //Will converge given a sufficiently accurate estimate for a root of f
   //Useful for refining answers to other root-finding methods
   //The more accurate central finite difference method will be used to calculate the derivative
   template <typename T,typename Function>
   T newton_raphson(Function f,T x0,T tolerance,std::size_t max_iterations)
   {
    detail::validate_positive(static_cast<long double>(tolerance),"tolerance");
    detail::validate_iterations(max_iterations);

    for(std::size_t i = 0 ; i < max_iterations; ++i)
    {
     //In essence, solve for x from the equation for a finite difference:
     //f'(x0) = (f(x)-f(x0))/(x-x0) where x is a root of f (i.e. f(x) = 0).
     T f_x0 = f(x0),fprime_x0 = diff::cdiff<T>(f,x0);
     if (std::abs(f_x0) <= tolerance)
      return x0;

     //If at a critical point that isn't a root, kill the algorithm
     if (fprime_x0 == 0)
      detail::throw_division_error("Newton-Raphson method");

     T x = x0-(f_x0/fprime_x0);
     if (std::abs(x-x0) <= tolerance)
      return x;

     //Set x0 = x
     x0 = x;
    }

    detail::throw_no_convergence("Newton-Raphson method");
   }

   //Newton-Raphson using the Horner method (for polynomials)
   template <typename T>
   T newton_raphson_horner(const std::vector<T>& f,T x0,T tolerance,std::size_t max_iterations)
   {
    detail::validate_positive(static_cast<long double>(tolerance),"tolerance");
    detail::validate_iterations(max_iterations);

    for(std::size_t i = 0 ; i < max_iterations; ++i)
    {
     //In essence, solve for x from the equation for a finite difference:
     //f'(x0) = (f(x)-f(x0))/(x-x0) where x is a root of f (i.e. f(x) = 0).
     std::pair<T,T> f_and_fprime = algebraic::horner(f,x0);
     T f_x0 = f_and_fprime.first,fprime_x0 = f_and_fprime.second;
     if (std::abs(f_x0) <= tolerance)
      return x0;

     //If at a critical point that isn't a root, kill the algorithm
     if (fprime_x0 == 0)
      detail::throw_division_error("Newton-Raphson-Horner method");

     T x = x0-(f_x0/fprime_x0);
     if (std::abs(x-x0) <= tolerance)
      return x;

     //Set x0 = x
     x0 = x;
    }

    detail::throw_no_convergence("Newton-Raphson-Horner method");
   }

   //Secant method: A variant of NR that replaces the derivative with an approximation and uses iteration
   //on x to solve for the roots.
   //Easier to calculate due to not using a derivative.
   template <typename T,typename Function>
   T secant(Function f,T x0,T x1,T tolerance,std::size_t max_iterations)
   {
    detail::validate_positive(static_cast<long double>(tolerance),"tolerance");
    detail::validate_iterations(max_iterations);

    T y0 = f(x0),y1 = f(x1);
    for(std::size_t i = 0; i < max_iterations;++i)
    {
     if (std::abs(y1) <= tolerance)
      return x1;

     //Discrete analogue to the finite difference in NR
     T denominator = y1-y0;
     if (denominator == 0)
      detail::throw_division_error("Secant method");

     T x = x1-(y1*(x1-x0)/denominator);
     if (std::abs(x-x1) <= tolerance)
      return x;

     //Update the values of each variable:
     //x1->x0, y1->y0  x->x1,f(x)->y1
     x0 = x1; y0 = y1;
     x1 = x; y1 = f(x);
    }

    detail::throw_no_convergence("Secant method");
   }

   //Steffensen's Method: Application of Aitken's Delta-squared method to the fixed-point method
   //Ensures quadratic convergence of the sequence to a solution
   template <typename T,typename Function>
   T steffensen(Function f,T x0,T tolerance,std::size_t max_iterations)
   {
    detail::validate_positive(static_cast<long double>(tolerance),"tolerance");
    detail::validate_iterations(max_iterations);

    for(std::size_t i = 0; i < max_iterations; ++i)
    {
     T x1 = f(x0),x2 = f(x1);
     T denominator = x2-static_cast<T>(2)*x1+x0;
     if (denominator == 0)
      detail::throw_division_error("Steffensen's method");

     T delta = x1-x0;
     T x = x0-(delta*delta/denominator);
     if(std::abs(x-x0) <= tolerance)
      return x;

     x0 = x;
    }

    detail::throw_no_convergence("Steffensen's method");
   }
  }
 }
}

#endif
