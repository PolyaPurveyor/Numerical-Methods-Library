#include "nml/nml.h"

#include <string>


namespace NML
{
 namespace detail
 {
  void validate_nonzero(long double value,const char* name)
  {
   if (value == 0.0L)
    throw std::invalid_argument(std::string(name)+" cannot equal 0");
  }

  void validate_positive(long double value,const char* name)
  {
   if (value <= 0.0L)
    throw std::invalid_argument(std::string(name)+" must be positive");
  }

  void validate_iterations(std::size_t max_iterations)
  {
   if (max_iterations == 0)
    throw std::invalid_argument("max_iterations must be greater than 0");
  }

  [[noreturn]] void throw_division_error(const char* method)
  {
   throw std::domain_error(std::string(method)+" encountered a zero-valued divisor");
  }

  [[noreturn]] void throw_no_convergence(const char* method)
  {
   throw std::runtime_error(std::string(method)+" failed after the maximum number of iterations");
  }
 }
}
