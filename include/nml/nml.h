/*
   Egwuchukwu Kalu
   Numerical Methods Library
   A personal implementation of common numerical methods for differentiation, integration,
   ordinary and partial differential equations
*/
#ifndef NML_NML_H
#define NML_NML_H         // Library use-case
                          // -------------------------------------------------
#include <cmath>          //Mathematical operations
#include <csignal>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>         //Operating on a list of values in sequence
#include <utility>        //Move semantics
#include <algorithm>      //Iterating over the ranges of functions
#include <exception>      //Throwing divide-by-zero and other exceptions
#include <complex>        //If I feel cute, adding complex support to this
#include <Eigen/Eigen>    //Partial differentiation support
#include <limits>
#include <type_traits>


namespace NML
{
 namespace detail
 {
  void validate_nonzero(long double value,const char* name);
  void validate_positive(long double value,const char* name);
  void validate_iterations(std::size_t max_iterations);
  [[noreturn]] void throw_division_error(const char* method);
  [[noreturn]] void throw_no_convergence(const char* method);
 }

 namespace utils
 {
  //Generates an interval from two floating point numbers
  template <typename T>
  std::vector<T> interval_gen(T x0,T xn,std::size_t n);
 }

 namespace functional
 {
  /*Some equation solving methods require a derivative,
   * so the finite difference is first*/
  namespace diff
  {
   /*Forward Difference Formula */
   template <typename T,typename Function>
   T fdiff(Function f,T x,T dx = static_cast<T>(1e-6));

   /*Backward Difference Formula, which is just forward differentiation with a negative dx lol*/
   template <typename T,typename Function>
   T bdiff(Function f,T x,T dx = static_cast<T>(1e-6));

   /*Central Difference Formula*/
   template <typename T,typename Function>
   T cdiff(Function f,T x,T dx = static_cast<T>(1e-6));
  }

  namespace algebraic
  {
   /*Horner's Method for Synthetic Division*/
   template <typename T>
   std::pair<T,T> horner(const std::vector<T>& P,T x0);

   /*Quadratic Synthetic-Division*/
   template <typename T>
   std::pair<T,T> qsd(const std::vector<T>& R,T B1,T B0);
  }

  //For all these methods, the functions are assumed to be continuous
  //(or differentiable if necessary as with Newton-Raphson) or otherwise not pathological
  namespace eqsoln
  {
   //Bisection Root-Finding Method: Finds a solution to f(x) = 0 given the continuous function f on the interval [a,b] where f(a) and f(b)
   //have opposite signs through a binary search where x = a+(b-a)/2
   template <typename T,typename Function>
   T bisection(Function f,T a,T b,T tolerance = static_cast<T>(1e-6),std::size_t max_iterations = 20);

   //Fixed Point Iteration: To find a solution to x = f(x), give an initial solution x0 and iteratively set x = f(x0)
   //and set x0 = x if not within the proper tolerance
   template <typename T,typename Function>
   T fixed_point(Function f,T x0,T tolerance = static_cast<T>(1e-6),std::size_t max_iterations = 10);

   //Newton-Raphson Method: f must be differentiable (twice for fast convergence)
   //Will converge given a sufficiently accurate estimate for a root of f
   //Useful for refining answers to other root-finding methods
   //The more accurate central finite difference method will be used to calculate the derivative
   template <typename T,typename Function>
   T newton_raphson(Function f,T x0,T tolerance = static_cast<T>(1e-6),std::size_t max_iterations = 10);

   //Newton-Raphson using the Horner method (for polynomials)
   template <typename T>
   T newton_raphson_horner(const std::vector<T>& f,T x0,T tolerance = static_cast<T>(1e-6),std::size_t max_iterations = 10);

   //Secant method: A variant of NR that replaces the derivative with an approximation and uses iteration
   //on x to solve for the roots.
   //Easier to calculate due to not using a derivative.
   template <typename T,typename Function>
   T secant(Function f,T x0,T x1,T tolerance = static_cast<T>(1e-6),std::size_t max_iterations = 10);

   //Steffensen's Method: Application of Aitken's Delta-squared method to the fixed-point method
   //Ensures quadratic convergence of the sequence to a solution
   template <typename T,typename Function>
   T steffensen(Function f,T x0,T tolerance = static_cast<T>(1e-6),std::size_t max_iterations = 10);
  }

  namespace integral
  {

  }

  namespace ode
  {

  }
 }
}

#include "nml/detail/nml_impl.h"

#endif
