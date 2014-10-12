//
// SHORT DESCRIPTION
// =================
//
// The class "mpfr::real" is a C++ interface to the GNU MPFR library
// version 3.0.0 or later.
//
// COPYRIGHT/LICENSE
// =================
//
// Copyright 2010, 2011 Christian Schneider <software(at)chschneider(dot)eu>
//
// Version: 0.0.2-alpha
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 3 of the License, not any later
// version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef MPFR_REAL_HPP
#define MPFR_REAL_HPP 1

// //////////////////////////////////////////////////////////////////
// inclusion of headers
// //////////////////////////////////////////////////////////////////

#include <mpfr.h>
#include <cmath>    // for returned values of fpclassify() etc.
#include <cstring>  // for isdigit() and strcpy() etc.
#include <iostream>
#include <limits>
#include <sstream>
#include <string>



namespace mpfr {

  // //////////////////////////////////////////////////////////////////
  // type definitions
  // //////////////////////////////////////////////////////////////////

  typedef mpfr_prec_t real_prec_t;
  typedef mp_exp_t    real_exp_t;
  typedef mpfr_rnd_t  real_rnd_t;

  // exception type

  class exception_real: public std::exception {
    public:
      exception_real(const std::string& msg = "exception_real")
        throw(): msg_(msg) {}
      virtual ~exception_real() throw() {}
      // returns cause of error
      virtual const char* what() const throw() { return msg_.c_str(); };
    private:
      std::string msg_;
  };

}  // namespace mpfr



namespace mpfr {

  // //////////////////////////////////////////////////////////////////
  // helper functions
  // //////////////////////////////////////////////////////////////////

  inline int helper_set_stdstr(mpfr_ptr rop,
                               const std::string& op,
                               mpfr_rnd_t rnd) {
    const int err = mpfr_set_str(rop, op.c_str(), 0, rnd);
    if (err == -1)
      throw exception_real(
        std::string("in mpfr::helper_set_stdstr(mpfr_ptr, const std::string&, mpfr_rnd_t):\n  invalid input format ")
        + op);
    return err;
  }

  inline int helper_set_charptr(mpfr_ptr rop,
                                const char* op,
                                mpfr_rnd_t rnd) {
    const int err = mpfr_set_str(rop, op, 0, rnd);
    if (err == -1)
      throw exception_real(
        std::string("in mpfr::helper_set_charptr(mpfr_ptr, const char*, mpfr_rnd_t):\n  invalid input format ")
        + op);
    return err;
  }

}  // namespace mpfr



namespace mpfr {

  // //////////////////////////////////////////////////////////////
  // declaration of real
  // //////////////////////////////////////////////////////////////

  template <real_prec_t prec_, real_rnd_t rnd_>
  class real;

  // //////////////////////////////////////////////////////////////////
  // type traits
  // //////////////////////////////////////////////////////////////////

  template <class Tp1_, class Tp2_>
  struct type_traits {
    typedef Tp1_ real_type;
    typedef Tp2_ other_type;

    // support level in class real
    static const bool enable_impl_ctor  = false;  // implicit ctor (else expl.)
    static const bool enable_assign_op  = false;  // assignment operator
    static const bool enable_conv_func  = false;  // conversion member function
    static const bool enable_arithm_ops = false;  // arithmetic operators
    static const bool enable_compar_ops = false;  // comparison operators
    static const bool enable_math_funcs = false;  // mathematical functions

    // support in MPFR library
    // (Note: "has_get_a" beats "has_get_b", if both are "true".)
    static const bool has_set   = false;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = false;
    static const bool has_sub_a = false;
    static const bool has_sub_b = false;
    static const bool has_mul   = false;
    static const bool has_div_a = false;
    static const bool has_div_b = false;
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, unsigned long int> {
    typedef real<prec_, rnd_> real_type;
    typedef unsigned long int other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const unsigned long int op, mpfr_rnd_t rnd) {
      return mpfr_set_ui(rop, op, rnd);
    }
    inline static unsigned long int get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_ui(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const unsigned long int op2, mpfr_rnd_t rnd) {
      return mpfr_add_ui(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned long int op2, mpfr_rnd_t rnd) {
      return mpfr_sub_ui(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const unsigned long int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_ui_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const unsigned long int op2, mpfr_rnd_t rnd) {
      return mpfr_mul_ui(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned long int op2, mpfr_rnd_t rnd) {
      return mpfr_div_ui(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const unsigned long int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_ui_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, long int> {
    typedef real<prec_, rnd_> real_type;
    typedef long int other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const long int op, mpfr_rnd_t rnd) {
      return mpfr_set_si(rop, op, rnd);
    }
    inline static long int get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_si(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const long int op2, mpfr_rnd_t rnd) {
      return mpfr_add_si(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const long int op2, mpfr_rnd_t rnd) {
      return mpfr_sub_si(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const long int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const long int op2, mpfr_rnd_t rnd) {
      return mpfr_mul_si(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const long int op2, mpfr_rnd_t rnd) {
      return mpfr_div_si(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const long int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, unsigned int> {
    typedef real<prec_, rnd_> real_type;
    typedef unsigned int other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const unsigned int op, mpfr_rnd_t rnd) {
      return mpfr_set_ui(rop, op, rnd);
    }
    inline static unsigned int get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_ui(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const unsigned int op2, mpfr_rnd_t rnd) {
      return mpfr_add_ui(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned int op2, mpfr_rnd_t rnd) {
      return mpfr_sub_ui(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const unsigned int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_ui_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const unsigned int op2, mpfr_rnd_t rnd) {
      return mpfr_mul_ui(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned int op2, mpfr_rnd_t rnd) {
      return mpfr_div_ui(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const unsigned int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_ui_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, int> {
    typedef real<prec_, rnd_> real_type;
    typedef int other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const int op, mpfr_rnd_t rnd) {
      return mpfr_set_si(rop, op, rnd);
    }
    inline static int get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_si(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const int op2, mpfr_rnd_t rnd) {
      return mpfr_add_si(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const int op2, mpfr_rnd_t rnd) {
      return mpfr_sub_si(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const int op2, mpfr_rnd_t rnd) {
      return mpfr_mul_si(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const int op2, mpfr_rnd_t rnd) {
      return mpfr_div_si(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, unsigned short int> {
    typedef real<prec_, rnd_> real_type;
    typedef unsigned short int other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const unsigned short int op, mpfr_rnd_t rnd) {
      return mpfr_set_ui(rop, op, rnd);
    }
    inline static unsigned short int get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_ui(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const unsigned short int op2, mpfr_rnd_t rnd) {
      return mpfr_add_ui(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned short int op2, mpfr_rnd_t rnd) {
      return mpfr_sub_ui(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const unsigned short int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_ui_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const unsigned short int op2, mpfr_rnd_t rnd) {
      return mpfr_mul_ui(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned short int op2, mpfr_rnd_t rnd) {
      return mpfr_div_ui(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const unsigned short int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_ui_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, short int> {
    typedef real<prec_, rnd_> real_type;
    typedef short int other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const short int op, mpfr_rnd_t rnd) {
      return mpfr_set_si(rop, op, rnd);
    }
    inline static short int get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_si(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const short int op2, mpfr_rnd_t rnd) {
      return mpfr_add_si(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const short int op2, mpfr_rnd_t rnd) {
      return mpfr_sub_si(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const short int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const short int op2, mpfr_rnd_t rnd) {
      return mpfr_mul_si(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const short int op2, mpfr_rnd_t rnd) {
      return mpfr_div_si(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const short int op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, unsigned char> {
    typedef real<prec_, rnd_> real_type;
    typedef unsigned char other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const unsigned char op, mpfr_rnd_t rnd) {
      return mpfr_set_ui(rop, op, rnd);
    }
    inline static unsigned char get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_ui(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const unsigned char op2, mpfr_rnd_t rnd) {
      return mpfr_add_ui(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned char op2, mpfr_rnd_t rnd) {
      return mpfr_sub_ui(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const unsigned char op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_ui_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const unsigned char op2, mpfr_rnd_t rnd) {
      return mpfr_mul_ui(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned char op2, mpfr_rnd_t rnd) {
      return mpfr_div_ui(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const unsigned char op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_ui_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, signed char> {
    typedef real<prec_, rnd_> real_type;
    typedef signed char other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const signed char op, mpfr_rnd_t rnd) {
      return mpfr_set_si(rop, op, rnd);
    }
    inline static signed char get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_si(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const signed char op2, mpfr_rnd_t rnd) {
      return mpfr_add_si(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const signed char op2, mpfr_rnd_t rnd) {
      return mpfr_sub_si(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const signed char op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const signed char op2, mpfr_rnd_t rnd) {
      return mpfr_mul_si(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const signed char op2, mpfr_rnd_t rnd) {
      return mpfr_div_si(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const signed char op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, char> {
    typedef real<prec_, rnd_> real_type;
    typedef char other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const char op, mpfr_rnd_t rnd) {
      return mpfr_set_si(rop, op, rnd);
    }
    inline static char get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_si(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const char op2, mpfr_rnd_t rnd) {
      return mpfr_add_si(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const char op2, mpfr_rnd_t rnd) {
      return mpfr_sub_si(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const char op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const char op2, mpfr_rnd_t rnd) {
      return mpfr_mul_si(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const char op2, mpfr_rnd_t rnd) {
      return mpfr_div_si(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const char op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, wchar_t> {
    typedef real<prec_, rnd_> real_type;
    typedef wchar_t other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const wchar_t op, mpfr_rnd_t rnd) {
      return mpfr_set_si(rop, op, rnd);
    }
    inline static wchar_t get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_si(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const wchar_t op2, mpfr_rnd_t rnd) {
      return mpfr_add_si(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const wchar_t op2, mpfr_rnd_t rnd) {
      return mpfr_sub_si(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const wchar_t op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const wchar_t op2, mpfr_rnd_t rnd) {
      return mpfr_mul_si(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const wchar_t op2, mpfr_rnd_t rnd) {
      return mpfr_div_si(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const wchar_t op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_si_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, float> {
    typedef real<prec_, rnd_> real_type;
    typedef float other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const float op, mpfr_rnd_t rnd) {
      return mpfr_set_flt(rop, op, rnd);
    }
    inline static float get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_flt(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const float op2, mpfr_rnd_t rnd) {
      return mpfr_add_d(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const float op2, mpfr_rnd_t rnd) {
      return mpfr_sub_d(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const float op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_d_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const float op2, mpfr_rnd_t rnd) {
      return mpfr_mul_d(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const float op2, mpfr_rnd_t rnd) {
      return mpfr_div_d(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const float op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_d_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, double> {
    typedef real<prec_, rnd_> real_type;
    typedef double other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const double op, mpfr_rnd_t rnd) {
      return mpfr_set_d(rop, op, rnd);
    }
    inline static double get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_d(op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd) {
      return mpfr_add_d(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd) {
      return mpfr_sub_d(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, const double op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_d_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd) {
      return mpfr_mul_d(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd) {
      return mpfr_div_d(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, const double op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_d_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, long double> {
    typedef real<prec_, rnd_> real_type;
    typedef long double other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = true;
    static const bool has_get_b = false;
    static const bool has_add   = false;
    static const bool has_sub_a = false;
    static const bool has_sub_b = false;
    static const bool has_mul   = false;
    static const bool has_div_a = false;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const long double op, mpfr_rnd_t rnd) {
      return mpfr_set_ld(rop, op, rnd);
    }
    inline static long double get_a(mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_ld(op, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpz_t> {
    typedef real<prec_, rnd_> real_type;
    typedef mpz_t other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = true;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = false;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpz_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set_z(rop, op, rnd);
    }
    inline static int get_b(mpz_t rop, mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_z(rop, op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_add_z(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub_z(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_mul_z(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div_z(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpq_t> {
    typedef real<prec_, rnd_> real_type;
    typedef mpq_t other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = false;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpq_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set_q(rop, op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_add_q(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub_q(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_mul_q(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div_q(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpf_t> {
    typedef real<prec_, rnd_> real_type;
    typedef mpf_t other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = true;
    static const bool has_add   = false;
    static const bool has_sub_a = false;
    static const bool has_sub_b = false;
    static const bool has_mul   = false;
    static const bool has_div_a = false;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpf_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set_f(rop, op, rnd);
    }
    inline static int get_b(mpf_t rop, mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_f(rop, op, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpfr_t> {
    typedef real<prec_, rnd_> real_type;
    typedef mpfr_t other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = true;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set(rop, op, rnd);
    }
    inline static int get_b(mpfr_t rop, mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set(rop, op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_add(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_mul(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpz_ptr> {
    typedef real<prec_, rnd_> real_type;
    typedef mpz_ptr other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = true;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = false;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpz_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set_z(rop, op, rnd);
    }
    inline static int get_b(mpz_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_z(rop, op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_add_z(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub_z(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_mul_z(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div_z(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpq_ptr> {
    typedef real<prec_, rnd_> real_type;
    typedef mpq_ptr other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = false;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpq_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set_q(rop, op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_add_q(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub_q(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_mul_q(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div_q(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpf_ptr> {
    typedef real<prec_, rnd_> real_type;
    typedef mpf_ptr other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = true;
    static const bool has_add   = false;
    static const bool has_sub_a = false;
    static const bool has_sub_b = false;
    static const bool has_mul   = false;
    static const bool has_div_a = false;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpf_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set_f(rop, op, rnd);
    }
    inline static int get_b(mpf_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_get_f(rop, op, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpfr_ptr> {
    typedef real<prec_, rnd_> real_type;
    typedef mpfr_ptr other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = true;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set(rop, op, rnd);
    }
    inline static int get_b(mpfr_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set(rop, op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_add(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_mul(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpz_srcptr> {
    typedef real<prec_, rnd_> real_type;
    typedef mpz_srcptr other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = false;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpz_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set_z(rop, op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_add_z(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub_z(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_mul_z(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div_z(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpq_srcptr> {
    typedef real<prec_, rnd_> real_type;
    typedef mpq_srcptr other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = false;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpq_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set_q(rop, op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_add_q(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub_q(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_mul_q(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div_q(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpf_srcptr> {
    typedef real<prec_, rnd_> real_type;
    typedef mpf_srcptr other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = false;
    static const bool has_sub_a = false;
    static const bool has_sub_b = false;
    static const bool has_mul   = false;
    static const bool has_div_a = false;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpf_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set_f(rop, op, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, mpfr_srcptr> {
    typedef real<prec_, rnd_> real_type;
    typedef mpfr_srcptr other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = true;
    static const bool has_sub_a = true;
    static const bool has_sub_b = true;
    static const bool has_mul   = true;
    static const bool has_div_a = true;
    static const bool has_div_b = true;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd) {
      return mpfr_set(rop, op, rnd);
    }
    inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_add(rop, op1, op2, rnd);
    }
    inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub(rop, op1, op2, rnd);
    }
    inline static int sub_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_sub(rop, op1, op2, rnd);
    }
    inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_mul(rop, op1, op2, rnd);
    }
    inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div(rop, op1, op2, rnd);
    }
    inline static int div_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd) {
      return mpfr_div(rop, op1, op2, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, std::string> {
    typedef real<prec_, rnd_> real_type;
    typedef std::string other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = false;
    static const bool has_sub_a = false;
    static const bool has_sub_b = false;
    static const bool has_mul   = false;
    static const bool has_div_a = false;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const std::string& op, mpfr_rnd_t rnd) {
      return helper_set_stdstr(rop, op, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, char*> {
    typedef real<prec_, rnd_> real_type;
    typedef char* other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = false;
    static const bool has_sub_a = false;
    static const bool has_sub_b = false;
    static const bool has_mul   = false;
    static const bool has_div_a = false;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const char* op, mpfr_rnd_t rnd) {
      return helper_set_charptr(rop, op, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_>
  struct type_traits<real<prec_, rnd_>, const char*> {
    typedef real<prec_, rnd_> real_type;
    typedef const char* other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = false;
    static const bool has_sub_a = false;
    static const bool has_sub_b = false;
    static const bool has_mul   = false;
    static const bool has_div_a = false;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const char* op, mpfr_rnd_t rnd) {
      return helper_set_charptr(rop, op, rnd);
    }
  };

  template <real_prec_t prec_, real_rnd_t rnd_, int i_>
  struct type_traits<real<prec_, rnd_>, char[i_]> {
    typedef real<prec_, rnd_> real_type;
    typedef char other_type[i_];

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    static const bool has_set   = true;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = false;
    static const bool has_sub_a = false;
    static const bool has_sub_b = false;
    static const bool has_mul   = false;
    static const bool has_div_a = false;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (must be defined if corresponding "has_..." boolean is set to "true")
    inline static int set(mpfr_ptr rop, const char* op, mpfr_rnd_t rnd) {
      return helper_set_charptr(rop, op, rnd);
    }
  };

  template <real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd_>
  struct type_traits<real<prec1_, rnd_>, real<prec2_, rnd_> > {
    typedef real<prec1_, rnd_> real_type;
    typedef real<prec2_, rnd_> other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    // (all "has_..." booleans *MUST* be set to "false" for "real<prec_, rnd_>")
    static const bool has_set   = false;
    static const bool has_get_a = false;
    static const bool has_get_b = false;
    static const bool has_add   = false;
    static const bool has_sub_a = false;
    static const bool has_sub_b = false;
    static const bool has_mul   = false;
    static const bool has_div_a = false;
    static const bool has_div_b = false;

    // functions in MPFR library
    // (there should be no function definitions for "real<prec_, rnd_>")
  };

}  // namespace mpfr



namespace mpfr {

  // //////////////////////////////////////////////////////////////
  // declaration of real
  // //////////////////////////////////////////////////////////////

  template <real_prec_t prec_, real_rnd_t rnd_>
  class real;

  // //////////////////////////////////////////////////////////////////
  // basic meta-programming
  // //////////////////////////////////////////////////////////////////

  // enable_if

  template<bool, class Tp_ = void>
  struct enable_if {};

  template<class Tp_>
  struct enable_if<true, Tp_> {
    typedef Tp_ type;
  };

  // //////////////////////////////////////////////////////////////////
  // calculated result type for two arguments
  // //////////////////////////////////////////////////////////////////

  // At least one argument must be of type "real" and all "real"s must have
  // the same rounding ("rnd_").  The result type is the "real" with highest
  // precision.

  template<class Tp1_, class Tp2_>
  struct result_type2 {
  };

  template <real_prec_t prec_, real_rnd_t rnd_, class Tp_>
  struct result_type2<real<prec_, rnd_>, Tp_> {
    typedef real<prec_, rnd_> type;
    static const real_rnd_t  rnd  = rnd_;
    static const real_prec_t prec = prec_;
  };

  template <real_prec_t prec_, real_rnd_t rnd_, class Tp_>
  struct result_type2<Tp_, real<prec_, rnd_> > {
    typedef real<prec_, rnd_> type;
    static const real_rnd_t  rnd  = rnd_;
    static const real_prec_t prec = prec_;
  };

  template <real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd_>
  struct result_type2<real<prec1_, rnd_>, real<prec2_, rnd_> > {
    typedef real<((prec1_ < prec2_) ? prec2_ : prec1_), rnd_> type;
    static const real_rnd_t  rnd  = rnd_;
    static const real_prec_t prec = ((prec1_ < prec2_) ? prec2_ : prec1_);
  };

  template <real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd1_, real_rnd_t rnd2_>
  struct result_type2<real<prec1_, rnd1_>, real<prec2_, rnd2_> > {
  };

  // //////////////////////////////////////////////////////////////////
  // calculated result type for three arguments
  // //////////////////////////////////////////////////////////////////

  // At least one argument must be of type "real" and all "real"s must have
  // the same rounding ("rnd_").  The result type is the "real" with highest
  // precision.

  template<class Tp1_, class Tp2_, class Tp3_>
  struct result_type3 {
  };

  template<real_prec_t prec_, real_rnd_t rnd_, class Tp1_, class Tp2_>
  struct result_type3<real<prec_, rnd_>, Tp1_, Tp2_> {
    typedef real<prec_, rnd_> type;
    static const real_rnd_t  rnd  = rnd_;
    static const real_prec_t prec = prec_;
  };

  template<real_prec_t prec_, real_rnd_t rnd_, class Tp1_, class Tp2_>
  struct result_type3<Tp2_, real<prec_, rnd_>, Tp1_> {
    typedef real<prec_, rnd_> type;
    static const real_rnd_t  rnd  = rnd_;
    static const real_prec_t prec = prec_;
  };

  template<real_prec_t prec_, real_rnd_t rnd_, class Tp1_, class Tp2_>
  struct result_type3<Tp1_, Tp2_, real<prec_, rnd_> > {
    typedef real<prec_, rnd_> type;
    static const real_rnd_t  rnd  = rnd_;
    static const real_prec_t prec = prec_;
  };

  template<real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd_, class Tp_>
  struct result_type3<real<prec1_, rnd_>, real<prec2_, rnd_>, Tp_> {
    typedef real<((prec1_ < prec2_) ? prec2_ : prec1_), rnd_> type;
    static const real_rnd_t  rnd  = rnd_;
    static const real_prec_t prec = ((prec1_ < prec2_) ? prec2_ : prec1_);
  };

  template<real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd1_, real_rnd_t rnd2_, class Tp_>
  struct result_type3<real<prec1_, rnd1_>, real<prec2_, rnd2_>, Tp_> {
  };

  template<real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd_, class Tp_>
  struct result_type3<Tp_, real<prec1_, rnd_>, real<prec2_, rnd_> > {
    typedef real<((prec1_ < prec2_) ? prec2_ : prec1_), rnd_> type;
    static const real_rnd_t  rnd  = rnd_;
    static const real_prec_t prec = ((prec1_ < prec2_) ? prec2_ : prec1_);
  };

  template<real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd1_, real_rnd_t rnd2_, class Tp_>
  struct result_type3<Tp_, real<prec1_, rnd1_>, real<prec2_, rnd2_> > {
  };

  template<real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd_, class Tp_>
  struct result_type3<real<prec2_, rnd_>, Tp_, real<prec1_, rnd_> > {
    typedef real<((prec1_ < prec2_) ? prec2_ : prec1_), rnd_> type;
    static const real_rnd_t  rnd  = rnd_;
    static const real_prec_t prec = ((prec1_ < prec2_) ? prec2_ : prec1_);
  };

  template<real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd1_, real_rnd_t rnd2_, class Tp_>
  struct result_type3<real<prec2_, rnd2_>, Tp_, real<prec1_, rnd1_> > {
  };

  template<real_prec_t prec1_, real_prec_t prec2_, real_prec_t prec3_, real_rnd_t rnd_>
  struct result_type3<real<prec1_, rnd_>, real<prec2_, rnd_>, real<prec3_, rnd_> > {
    typedef real<((((prec3_ < prec2_) ? prec2_ : prec3_) < prec1_) ? prec1_ :
                    ((prec3_ < prec2_) ? prec2_ : prec3_)), rnd_> type;
    static const real_rnd_t  rnd  = rnd_;
    static const real_prec_t prec = ((((prec3_ < prec2_) ? prec2_ : prec3_)
      < prec1_) ? prec1_ : ((prec3_ < prec2_) ? prec2_ : prec3_));
  };

  template<real_prec_t prec1_, real_prec_t prec2_, real_prec_t prec3_,
    real_rnd_t rnd1_, real_rnd_t rnd2_, real_rnd_t rnd3_>
  struct result_type3<real<prec1_, rnd1_>, real<prec2_, rnd2_>, real<prec3_, rnd3_> > {
  };

  // //////////////////////////////////////////////////////////////////
  // promotion to real
  // //////////////////////////////////////////////////////////////////

  template <class Tp1_, class Tp2_>
  struct promote {
  };

  template <real_prec_t prec_, real_rnd_t rnd_, class Tp_>
  struct promote<real<prec_, rnd_>, Tp_> {
    typedef const real<prec_, rnd_> type;
  };

  template <real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd_>
  struct promote<real<prec1_, rnd_>, real<prec2_, rnd_> > {
    typedef const real<prec2_, rnd_>& type;
  };

  // //////////////////////////////////////////////////////////////////
  // check for equal types
  // //////////////////////////////////////////////////////////////////

  template <class Tp1_, class Tp2_>
  struct equal_types2 {
    static const bool val = false;
  };

  template <class Tp_>
  struct equal_types2<Tp_, Tp_> {
    static const bool val = true;
  };

  // //////////////////////////////////////////////////////////////////
  // check for type "real"
  // //////////////////////////////////////////////////////////////////

  template <class Tp1_, class Tp2_>
  struct has_real2 {
    static const bool val = false;
  };

  template <real_prec_t prec_, real_rnd_t rnd_, class Tp_>
  struct has_real2<real<prec_, rnd_>, Tp_> {
    static const bool val = true;
  };

  template <real_prec_t prec_, real_rnd_t rnd_, class Tp_>
  struct has_real2<Tp_, real<prec_, rnd_> > {
    static const bool val = true;
  };

  template <real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd_>
  struct has_real2<real<prec1_, rnd_>, real<prec2_, rnd_> > {
    static const bool val = true;
  };

  template <real_prec_t prec1_, real_prec_t prec2_, real_rnd_t rnd1_, real_rnd_t rnd2_>
  struct has_real2<real<prec1_, rnd1_>, real<prec2_, rnd2_> > {
    static const bool val = false;
  };

}  // namespace mpfr




namespace mpfr {

  // //////////////////////////////////////////////////////////////
  // class declaration
  // //////////////////////////////////////////////////////////////

  template <real_prec_t prec_, real_rnd_t rnd_>
  class real;

  // //////////////////////////////////////////////////////////////
  // generic operators (definitions of binary operators)
  // //////////////////////////////////////////////////////////////

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_arithm_ops &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_arithm_ops,
    const typename result_type2<Tp1_, Tp2_>::type>::type
  operator +(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_add(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_arithm_ops &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_arithm_ops,
    const typename result_type2<Tp1_, Tp2_>::type>::type
  operator -(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_sub(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_arithm_ops &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_arithm_ops,
    const typename result_type2<Tp1_, Tp2_>::type>::type
  operator *(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_mul(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_arithm_ops &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_arithm_ops,
    const typename result_type2<Tp1_, Tp2_>::type>::type
  operator /(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_div(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_compar_ops &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_compar_ops,
    const bool>::type
  operator ==(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_equal_p(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_compar_ops &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_compar_ops,
    const bool>::type
  operator !=(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_lessgreater_p(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_compar_ops &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_compar_ops,
    const bool>::type
  operator <(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_less_p(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_compar_ops &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_compar_ops,
    const bool>::type
  operator <=(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_lessequal_p(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_compar_ops &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_compar_ops,
    const bool>::type
  operator >(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_greater_p(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_compar_ops &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_compar_ops,
    const bool>::type
  operator >=(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_greaterequal_p(temp1.x_, temp2.x_);
  }

  // //////////////////////////////////////////////////////////////
  // mathematical functions (definitions for multiple "real" arguments)
  // //////////////////////////////////////////////////////////////

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    const int>::type
  isgreater(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_greater_p(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    const int>::type
  isgreaterequal(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_greaterequal_p(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    const int>::type
  isless(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_less_p(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    const int>::type
  islessequal(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_lessequal_p(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    const int>::type
  islessgreater(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_lessgreater_p(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    const int>::type
  isunordered(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_unordered_p(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  atan2(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_atan2(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  copysign(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_copysign(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  fdim(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_dim(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_, class Tp3_>
  inline typename enable_if<
    type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
      Tp2_>::enable_math_funcs &&
    type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
      Tp3_>::enable_math_funcs,
    typename result_type3<Tp1_, Tp2_, Tp3_>::type>::type
  fma(const Tp1_& r1, const Tp2_& r2, const Tp3_& r3) {
    typedef typename result_type3<Tp1_, Tp2_, Tp3_>::type temp_type;
    const real_rnd_t rnd = result_type3<Tp1_, Tp2_, Tp3_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    typename promote<temp_type, Tp3_>::type temp3(r3);
    temp_type temp;
    mpfr_fma(temp.x_, temp1.x_, temp2.x_, temp3.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  fmax(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_max(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  fmin(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_min(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  fmod(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_fmod(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  hypot(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_hypot(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  pow(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_pow(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  remainder(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_remainder(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  agm(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_agm(temp.x_, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    const int>::type
  cmpabs(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    return mpfr_cmpabs(temp1.x_, temp2.x_);
  }

  template <class Tp1_, class Tp2_, class Tp3_>
  inline typename enable_if<
    type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
      Tp2_>::enable_math_funcs &&
    type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
      Tp3_>::enable_math_funcs,
    typename result_type3<Tp1_, Tp2_, Tp3_>::type>::type
  fms(const Tp1_& r1, const Tp2_& r2, const Tp3_& r3) {
    typedef typename result_type3<Tp1_, Tp2_, Tp3_>::type temp_type;
    const real_rnd_t rnd = result_type3<Tp1_, Tp2_, Tp3_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    typename promote<temp_type, Tp3_>::type temp3(r3);
    temp_type temp;
    mpfr_fms(temp.x_, temp1.x_, temp2.x_, temp3.x_, rnd);
    return temp;
  }

  template <real_prec_t prec_, real_rnd_t rnd_, class Tp_>
  inline typename enable_if<
    type_traits<typename result_type2<real<prec_, rnd_>, Tp_>::type,
      Tp_>::enable_math_funcs &&
    type_traits<typename result_type2<real<prec_, rnd_>, Tp_>::type,
      real<prec_, rnd_> >::enable_math_funcs,
    typename result_type2<real<prec_, rnd_>, Tp_>::type>::type
  modf(const Tp_& r, real<prec_, rnd_>* iptr) {
    typedef typename result_type2<real<prec_, rnd_>, Tp_>::type temp_type;
    const real_rnd_t rnd = result_type2<real<prec_, rnd_>, Tp_>::rnd;
    temp_type temp;
    mpfr_modf(iptr->x_, temp.x_, r.x_, rnd);
    return temp;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  nextafter(const Tp1_& r1, const Tp2_& r2) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    temp_type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    mpfr_nexttoward(temp1.x_, temp2.x_);
    return temp1;
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  nexttoward(const Tp1_& r1, const Tp2_& r2) {
    return nextafter(r1, r2);
  }

  template <class Tp1_, class Tp2_>
  inline typename enable_if<
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp1_>::enable_math_funcs &&
    type_traits<typename result_type2<Tp1_, Tp2_>::type,
      Tp2_>::enable_math_funcs,
    typename result_type2<Tp1_, Tp2_>::type>::type
  remquo(const Tp1_& r1, const Tp2_& r2, long* quo) {
    typedef typename result_type2<Tp1_, Tp2_>::type temp_type;
    const real_rnd_t rnd = result_type2<Tp1_, Tp2_>::rnd;
    typename promote<temp_type, Tp1_>::type temp1(r1);
    typename promote<temp_type, Tp2_>::type temp2(r2);
    temp_type temp;
    mpfr_remquo(temp.x_, quo, temp1.x_, temp2.x_, rnd);
    return temp;
  }

  // //////////////////////////////////////////////////////////////////
  // class definition
  // //////////////////////////////////////////////////////////////////

  template <real_prec_t prec_ = 53, real_rnd_t rnd_ = MPFR_RNDN>
  class real {
    // private:
    public:
      mpfr_t x_;

    public:
      // //////////////////////////////////////////////////////////////
      // default and copy constructors, default assignment operator, destructor
      // //////////////////////////////////////////////////////////////

      // default and copy constructor

      inline real() {
        mpfr_init2(x_, prec_);
        mpfr_set_zero(x_, +1);
      }

      inline real(const real& o) {
        mpfr_init2(x_, prec_);
        mpfr_set(x_, o.x_, rnd_);
      }

      // default assignment operator

      inline real& operator =(const real& o) {
        if (&o != this)
          mpfr_set(x_, o.x_, rnd_);
        return *this;
      }

      // destructor

      inline ~real() {
        mpfr_clear(x_);
      }

      // //////////////////////////////////////////////////////////////
      // converting constructors and converting assignment operators
      // //////////////////////////////////////////////////////////////

      // friend of other reals

      template <real_prec_t prec1_, real_rnd_t rnd1_>
      friend class real;

      // implicit conversion constructors

      template <class Tp_>
      inline real(const Tp_& o,
          typename enable_if<type_traits<real, Tp_>::has_set &&
          type_traits<real, Tp_>::enable_impl_ctor>::type* = 0) {
        mpfr_init2(x_, prec_);
        type_traits<real, Tp_>::set(x_, o, rnd_);
      }

      template <real_prec_t prec1_, real_rnd_t rnd1_>
      inline real(const real<prec1_, rnd1_>& o,
          typename enable_if<
          type_traits<real, real<prec1_, rnd1_> >::enable_impl_ctor>::type* = 0) {
        mpfr_init2(x_, prec1_);
        mpfr_set(x_, o.x_, rnd1_);
      }

      // explicit conversion constructors

      template <class Tp_>
      inline explicit real(const Tp_& o,
          typename enable_if<
          type_traits<real, Tp_>::has_set &&
          (! type_traits<real, Tp_>::enable_impl_ctor)>::type* = 0) {
        mpfr_init2(x_, prec_);
        type_traits<real, Tp_>::set(x_, o, rnd_);
      }

      template <real_prec_t prec1_, real_rnd_t rnd1_>
      inline explicit real(const real<prec1_, rnd1_>& o,
          typename enable_if<
          (! type_traits<real, real<prec1_, rnd1_> >::enable_impl_ctor)>::type* = 0) {
        mpfr_init2(x_, prec1_);
        mpfr_set(x_, o.x_, rnd1_);
      }

      // converting assignment operators

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_assign_op &&
        type_traits<real, Tp_>::has_set,
        real&>::type
      operator =(const Tp_& o) {
        type_traits<real, Tp_>::set(x_, o, rnd_);
        return *this;
      }

      template <real_prec_t prec1_, real_rnd_t rnd1_>
      inline typename enable_if<
        type_traits<real, real<prec1_, rnd1_> >::enable_assign_op,
        real&>::type
      operator =(const real<prec1_, rnd1_>& o) {
        mpfr_set(x_, o.x_, rnd1_);
        return *this;
      }

      // //////////////////////////////////////////////////////////////
      // generic operators
      // //////////////////////////////////////////////////////////////

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        (! type_traits<real, Tp_>::has_add),
        real&>::type
      operator +=(const Tp_& o) {
        typename promote<real, Tp_>::type temp(o);
        mpfr_add(x_, x_, temp.x_, rnd_);
        return *this;
      }

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        (! type_traits<real, Tp_>::has_sub_a),
        real&>::type
      operator -=(const Tp_& o) {
        typename promote<real, Tp_>::type temp(o);
        mpfr_sub(x_, x_, temp.x_, rnd_);
        return *this;
      }

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        (! type_traits<real, Tp_>::has_mul),
        real&>::type
      operator *=(const Tp_& o) {
        typename promote<real, Tp_>::type temp(o);
        mpfr_mul(x_, x_, temp.x_, rnd_);
        return *this;
      }

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        (! type_traits<real, Tp_>::has_div_a),
        real&>::type
      operator /=(const Tp_& o) {
        typename promote<real, Tp_>::type temp(o);
        mpfr_div(x_, x_, temp.x_, rnd_);
        return *this;
      }

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_arithm_ops &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_arithm_ops,
        const typename result_type2<Tp1_, Tp2_>::type>::type
      operator +(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_arithm_ops &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_arithm_ops,
        const typename result_type2<Tp1_, Tp2_>::type>::type
      operator -(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_arithm_ops &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_arithm_ops,
        const typename result_type2<Tp1_, Tp2_>::type>::type
      operator *(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_arithm_ops &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_arithm_ops,
        const typename result_type2<Tp1_, Tp2_>::type>::type
      operator /(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_compar_ops &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_compar_ops,
        const bool>::type
      operator ==(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_compar_ops &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_compar_ops,
        const bool>::type
      operator !=(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_compar_ops &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_compar_ops,
        const bool>::type
      operator <(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_compar_ops &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_compar_ops,
        const bool>::type
      operator <=(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_compar_ops &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_compar_ops,
        const bool>::type
      operator >(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_compar_ops &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_compar_ops,
        const bool>::type
      operator >=(const Tp1_& r1, const Tp2_& r2);

      // //////////////////////////////////////////////////////////////
      // optimized operators
      // //////////////////////////////////////////////////////////////

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_add,
        real&>::type
      operator +=(const Tp_& o) {
        type_traits<real, Tp_>::add(x_, x_, o, rnd_);
        return *this;
      }

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_sub_a,
        real&>::type
      operator -=(const Tp_& o) {
        type_traits<real, Tp_>::sub_a(x_, x_, o, rnd_);
        return *this;
      }

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_mul,
        real&>::type
      operator *=(const Tp_& o) {
        type_traits<real, Tp_>::mul(x_, x_, o, rnd_);
        return *this;
      }

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_div_a,
        real&>::type
      operator /=(const Tp_& o) {
        type_traits<real, Tp_>::div_a(x_, x_, o, rnd_);
        return *this;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_add,
        const real>::type
      operator +(const real& r1, const Tp_& r2) {
        real temp;
        type_traits<real, Tp_>::add(temp.x_, r1.x_, r2, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_add,
        const real>::type
      operator +(const Tp_& r1, const real& r2) {
        real temp;
        type_traits<real, Tp_>::add(temp.x_, r2.x_, r1, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_sub_a,
        const real>::type
      operator -(const real& r1, const Tp_& r2) {
        real temp;
        type_traits<real, Tp_>::sub_a(temp.x_, r1.x_, r2, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_sub_b,
        const real>::type
      operator -(const Tp_& r1, const real& r2) {
        real temp;
        type_traits<real, Tp_>::sub_b(temp.x_, r1, r2.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_mul,
        const real>::type
      operator *(const real& r1, const Tp_& r2) {
        real temp;
        type_traits<real, Tp_>::mul(temp.x_, r1.x_, r2, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_mul,
        const real>::type
      operator *(const Tp_& r1, const real& r2) {
        real temp;
        type_traits<real, Tp_>::mul(temp.x_, r2.x_, r1, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_div_a,
        const real>::type
      operator /(const real& r1, const Tp_& r2) {
        real temp;
        type_traits<real, Tp_>::div_a(temp.x_, r1.x_, r2, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_arithm_ops &&
        type_traits<real, Tp_>::has_div_b,
        const real>::type
      operator /(const Tp_& r1, const real& r2) {
        real temp;
        type_traits<real, Tp_>::div_b(temp.x_, r1, r2.x_, rnd_);
        return temp;
      }

      // //////////////////////////////////////////////////////////////
      // conversion operators and functions
      // //////////////////////////////////////////////////////////////

      // conversion operators
      // (can be enabled with preprocessor macro)

      #ifdef REAL_ENABLE_CONVERSION_OPERATORS
      inline operator unsigned long int() const {
        return mpfr_get_ui(x_, rnd_);
      }

      inline operator long int() const {
        return mpfr_get_si(x_, rnd_);
      }

      inline operator unsigned int() const {
        return mpfr_get_ui(x_, rnd_);
      }

      inline operator int() const {
        return mpfr_get_si(x_, rnd_);
      }

      inline operator unsigned short int() const {
        return mpfr_get_ui(x_, rnd_);
      }

      inline operator short int() const {
        return mpfr_get_si(x_, rnd_);
      }

      inline operator unsigned char() const {
        return mpfr_get_ui(x_, rnd_);
      }

      inline operator signed char() const {
        return mpfr_get_si(x_, rnd_);
      }

      inline operator char() const {
        return mpfr_get_si(x_, rnd_);
      }

      inline operator wchar_t() const {
        return mpfr_get_si(x_, rnd_);
      }

      inline operator float() const {
        return mpfr_get_flt(x_, rnd_);
      }

      inline operator double() const {
        return mpfr_get_d(x_, rnd_);
      }

      inline operator long double() const {
        return mpfr_get_ld(x_, rnd_);
      }

      inline operator std::string() const {
        std::stringstream temp;
        temp.precision(-1);
        try {
          temp << *this;
        }
        catch (...) {
          throw exception_real(
            "in real<prec_, rnd_>& real<prec_, rnd_>::operator std::string() const:\n  conversion failed");
        }
        return temp.str();
      }
      #endif  // REAL_ENABLE_CONVERSION_OPERATORS

      // conversion functions

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_conv_func &&
        type_traits<real, Tp_>::has_get_a,
        void>::type
      conv(Tp_& o) const {
        o = type_traits<real, Tp_>::get_a(x_, rnd_);
      }

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_conv_func &&
        (! type_traits<real, Tp_>::has_get_a) &&
        type_traits<real, Tp_>::has_get_b,
        void>::type
      conv(Tp_ const o) const {
        type_traits<real, Tp_>::get_b(o, x_, rnd_);
      }

      template <class Tp_>
      inline typename enable_if<
        type_traits<real, Tp_>::enable_conv_func &&
        equal_types2<char*, Tp_>::val,
        void>::type
      conv(Tp_ const o) const {
        std::stringstream temp;
        temp.precision(-1);
        try {
          temp << *this;
        }
        catch (...) {
          throw exception_real(
            "in const char* real<prec_, rnd_>::c_str() const:\n  conversion failed");
        }
        strcpy(o, temp.str().c_str());
      }

      // //////////////////////////////////////////////////////////////
      // increment, decrement, and negation operators
      // //////////////////////////////////////////////////////////////

      // increment operators

      inline real& operator ++() {
        static const real<prec_, rnd_> one_(1);
        *this += one_;
        return *this;
      }

      inline const real operator ++(int) {
        real<prec_, rnd_> temp = *this;
        ++(*this);
        return temp;
      }

      // decrement operators

      inline real& operator --() {
        static const real<prec_, rnd_> one_(1);
        *this -= one_;
        return *this;
      }

      inline const real operator --(int) {
        real<prec_, rnd_> temp = *this;
        --(*this);
        return temp;
      }

      // NOTE: The unary member operator- is declared after any template
      // binary friend operator-, because the latter may be unqualified
      // in the code above.  This way we make sure that binary - operations
      // do not match the unary member operator- (in any case).

      inline const real operator -() const {
        real<prec_, rnd_> temp;
        mpfr_neg(temp.x_, x_, rnd_);
        return temp;
      }

      // //////////////////////////////////////////////////////////////
      // std::istream and std::ostream operators
      // //////////////////////////////////////////////////////////////

      friend inline std::istream& operator >>(std::istream& s, real& r) {
        char c;
        std::string t;

        s >> std::ws;
        while ((c = s.get())) {
          if(c >= 0 && c <= 32) {  // a bit too much
            s.putback(c);
            break;
          }
          t += c;
        }

        try {
          r = t;
        }
        catch (...) {
          throw exception_real(
            std::string("in std::istream& operator >>(std::istream& s, real<prec_, rnd_>& r):\n  invalid input format ")
            + t);
        }

        return s;
      }

      // there might be some room for improvements for the next
      // missing: handling of ios_base::fixed/ios_base::scientific and
      // ios_base::showpoint

      friend inline std::ostream& operator <<(std::ostream& s, const real& r) {
        real_exp_t exp;
        char* ch = mpfr_get_str(0, &exp, 10, s.precision() + 1, r.x_, rnd_);
        if (! ch)
          throw exception_real(
            "in std::ostream& operator <<(std::ostream& s, const real<prec_, rnd_>& r):\n  conversion failed");
        std::string t = ch;
        mpfr_free_str(ch);

        const std::ios_base::fmtflags flags = s.flags();
        std::string::iterator t_iter = t.begin();

        if (*t_iter == '-')
          t_iter++;

        if (isdigit(*t_iter)) {  // normal number?
          // positive sign
          if ((t_iter == t.begin()) && (flags & std::ios_base::showpos)) {
            t_iter = t.insert(t_iter, '+');
            t_iter++;
          }

          // decimal point
          t_iter++;
          t.insert(t_iter, '.');

          // fixing exponent after insertion of decimal point
          // why must life be so difficult? (any suggestions for improvements?)
          if (! mpfr_zero_p(r.x_)) {
            const real_exp_t exp_prev = exp;
            volatile real_exp_t* exp_ptr = &exp;
            exp--;
            if (*exp_ptr > exp_prev)
              throw exception_real(
                "in std::ostream& operator <<(std::ostream& s, const real<prec_, rnd_>& r):\n  exponent out of range");
          }

          // composing of the exponent
          if (flags & std::ios_base::uppercase)
            t += 'E';
          else
            t += 'e';
          if (exp >= 0)
            t += '+';
          else {
            t += '-';
            exp = -exp;
          }
          if (exp >= -9 && exp <= 9)
            t += '0';
          std::stringstream temp;
          temp << exp;
          t += temp.str();
        }

        // width and adjustment
        if (s.width() > 0 && static_cast<unsigned int>(s.width()) > t.size()) {
          if (flags & std::ios_base::left)
            t_iter = t.end();
          else if (flags & std::ios_base::internal) {
            t_iter = t.begin();
            if (*t_iter == '+' || *t_iter == '-')
              t_iter++;
          }
          else
            t_iter = t.begin();
          while (t.size() < static_cast<unsigned int>(s.width()))
            t_iter = t.insert(t_iter, s.fill());
        }

        s << t;

        return s;
      }

      // //////////////////////////////////////////////////////////////
      // mathematical functions with zero or one "real" argument
      // //////////////////////////////////////////////////////////////

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const int>::type
      isfinite(const Tp_& r) {
        return mpfr_number_p(r.x_);
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const int>::type
      isinf(const Tp_& r) {
        return mpfr_inf_p(r.x_);
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const int>::type
      isnan(const Tp_& r) {
        return mpfr_nan_p(r.x_);
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const int>::type
      isnormal(const Tp_& r) {
        return mpfr_regular_p(r.x_);
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const int>::type
      signbit(const Tp_& r) {
        return mpfr_signbit(r.x_);
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      acos(const Tp_& r) {
        real temp;
        mpfr_acos(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      acosh(const Tp_& r) {
        real temp;
        mpfr_acosh(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      asin(const Tp_& r) {
        real temp;
        mpfr_asin(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      asinh(const Tp_& r) {
        real temp;
        mpfr_asinh(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      atan(const Tp_& r) {
        real temp;
        mpfr_atan(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      atanh(const Tp_& r) {
        real temp;
        mpfr_atanh(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      cbrt(const Tp_& r) {
        real temp;
        mpfr_cbrt(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      ceil(const Tp_& r) {
        real temp;
        mpfr_ceil(temp.x_, r.x_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      cos(const Tp_& r) {
        real temp;
        mpfr_cos(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      cosh(const Tp_& r) {
        real temp;
        mpfr_cosh(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      erf(const Tp_& r) {
        real temp;
        mpfr_erf(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      erfc(const Tp_& r) {
        real temp;
        mpfr_erfc(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      exp(const Tp_& r) {
        real temp;
        mpfr_exp(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      exp2(const Tp_& r) {
        real temp;
        mpfr_exp2(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      expm1(const Tp_& r) {
        real temp;
        mpfr_expm1(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      fabs(const Tp_& r) {
        real temp;
        mpfr_abs(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      abs(const Tp_& r) {
        real temp;
        mpfr_abs(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      floor(const Tp_& r) {
        real temp;
        mpfr_floor(temp.x_, r.x_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      log(const Tp_& r) {
        real temp;
        mpfr_log(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      log10(const Tp_& r) {
        real temp;
        mpfr_log10(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      log1p(const Tp_& r) {
        real temp;
        mpfr_log1p(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      log2(const Tp_& r) {
        real temp;
        mpfr_log2(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      nearbyint(const Tp_& r) {
        real temp;
        mpfr_rint(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      rint(const Tp_& r) {
        real temp;
        mpfr_rint(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      round(const Tp_& r) {
        real temp;
        mpfr_round(temp.x_, r.x_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      sin(const Tp_& r) {
        real temp;
        mpfr_sin(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      sinh(const Tp_& r) {
        real temp;
        mpfr_sinh(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      sqrt(const Tp_& r) {
        real temp;
        mpfr_sqrt(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      tan(const Tp_& r) {
        real temp;
        mpfr_tan(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      tanh(const Tp_& r) {
        real temp;
        mpfr_tanh(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      tgamma(const Tp_& r) {
        real temp;
        mpfr_gamma(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      trunc(const Tp_& r) {
        real temp;
        mpfr_trunc(temp.x_, r.x_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      j0(const Tp_& r) {
        real temp;
        mpfr_j0(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      j1(const Tp_& r) {
        real temp;
        mpfr_j1(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      y0(const Tp_& r) {
        real temp;
        mpfr_y0(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      y1(const Tp_& r) {
        real temp;
        mpfr_y1(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      ai(const Tp_& r) {
        real temp;
        mpfr_ai(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      cot(const Tp_& r) {
        real temp;
        mpfr_cot(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      coth(const Tp_& r) {
        real temp;
        mpfr_coth(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      csc(const Tp_& r) {
        real temp;
        mpfr_csc(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      csch(const Tp_& r) {
        real temp;
        mpfr_csch(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      digamma(const Tp_& r) {
        real temp;
        mpfr_digamma(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      exp10(const Tp_& r) {
        real temp;
        mpfr_exp10(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      expint(const Tp_& r) {
        real temp;
        mpfr_eint(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      frac(const Tp_& r) {
        real temp;
        mpfr_frac(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const int>::type
      isinteger(const Tp_& r) {
        return mpfr_integer_p(r.x_);
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const int>::type
      iszero(const Tp_& r) {
        return mpfr_zero_p(r.x_);
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      li2(const Tp_& r) {
        real temp;
        mpfr_li2(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      rec_sqrt(const Tp_& r) {
        real temp;
        mpfr_rec_sqrt(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      sec(const Tp_& r) {
        real temp;
        mpfr_sec(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      sech(const Tp_& r) {
        real temp;
        mpfr_sech(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const int>::type
      sgn(const Tp_& r) {
        return mpfr_sgn(r.x_);
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real>::type
      zeta(const Tp_& r) {
        real temp;
        mpfr_zeta(temp.x_, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const int>::type
      fpclassify(const Tp_& r) {
        if (mpfr_nan_p(r.x_))
          return FP_NAN;
        else if (mpfr_inf_p(r.x_))
          return FP_INFINITE;
        else if (mpfr_zero_p(r.x_))
          return FP_ZERO;
        else
          return FP_NORMAL;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      frexp(const Tp_& r, real_exp_t* exp) {
        if (mpfr_zero_p(r.x_)) {
          *exp = 0;
          return r;
        }
        else if (mpfr_inf_p(r.x_) || mpfr_nan_p(r.x_)) {
          //*exp = 0;
          return r;
        }
        else {
          Tp_ temp = r;
          *exp = mpfr_get_exp(r.x_);
          mpfr_set_exp(temp.x_, 0);
          return temp;
        }
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const real_exp_t>::type
      ilogb(const Tp_& r) {
        if (mpfr_zero_p(r.x_) || mpfr_nan_p(r.x_))
          return std::numeric_limits<real_exp_t>::min();
        else if (mpfr_inf_p(r.x_))
          return std::numeric_limits<real_exp_t>::max();
        else {
          real_exp_t temp = mpfr_get_exp(r.x_);
          if (temp != std::numeric_limits<real_exp_t>::min())
            temp--;
          return temp;
        }
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_ >::type
      ldexp(const Tp_& r, const long exp) {
        Tp_ temp;
        mpfr_mul_2si(temp.x_, r.x_, exp, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      lgamma(const Tp_& r) {
        Tp_ temp;
        int signp;
        mpfr_lgamma(temp.x_, &signp, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      logb(const Tp_& r) {
        Tp_ temp;
        if (mpfr_zero_p(r.x_))
          mpfr_set_inf(temp.x_, -1);
        else if (mpfr_nan_p(r.x_))
          mpfr_set_nan(temp.x_);
        else if (mpfr_inf_p(r.x_))
          mpfr_set_inf(temp.x_, 1);
        else {
          temp = mpfr_get_exp(r.x_);
          temp--;
        }
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      nan(const char*) {
        Tp_ temp;
        mpfr_set_nan(temp.x_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      scalbln(const Tp_& r, const long exp) {
        return ldexp(r, exp);  // FLT_RADIX == 2???
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      scalbn(const Tp_& r, const int exp) {
        return ldexp(r, exp);  // FLT_RADIX == 2???
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      jn(const long n, const Tp_& r) {
        Tp_ temp;
        mpfr_jn(temp.x_, n, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      yn(const long n, const Tp_& r) {
        Tp_ temp;
        mpfr_yn(temp.x_, n, r.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      cyl_bessel_j(const long n, const Tp_& r) {
        return jn(n, r);
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      cyl_neumann(const long n, const Tp_& r) {
        return yn(n, r);
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      factorial(const long n) {
        Tp_ temp;
        mpfr_fac_ui(temp.x_, n, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      root(const Tp_& r, const unsigned long n) {
        Tp_ temp;
        mpfr_root(temp.x_, r.x_, n, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      inf(const int n) {
        Tp_ temp;
        mpfr_set_inf(temp.x_, n);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      zero(const int n) {
        Tp_ temp;
        mpfr_set_zero(temp.x_, n);
        return temp;
      }

      // //////////////////////////////////////////////////////////////
      // mathematical functions with multiple "real" arguments
      // //////////////////////////////////////////////////////////////

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        const int>::type
      isgreater(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        const int>::type
      isgreaterequal(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        const int>::type
      isless(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        const int>::type
      islessequal(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        const int>::type
      islessgreater(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        const int>::type
      isunordered(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      atan2(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      copysign(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      fdim(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_, class Tp3_>
      friend typename enable_if<
        type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
          Tp2_>::enable_math_funcs &&
        type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
          Tp3_>::enable_math_funcs,
        typename result_type3<Tp1_, Tp2_, Tp3_>::type>::type
      fma(const Tp1_& r1, const Tp2_& r2, const Tp3_& r3);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      fmax(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      fmin(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      fmod(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      hypot(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      pow(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      remainder(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      agm(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        const int>::type
      cmpabs(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_, class Tp3_>
      friend typename enable_if<
        type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
          Tp2_>::enable_math_funcs &&
        type_traits<typename result_type3<Tp1_, Tp2_, Tp3_>::type,
          Tp3_>::enable_math_funcs,
        typename result_type3<Tp1_, Tp2_, Tp3_>::type>::type
      fms(const Tp1_& r1, const Tp2_& r2, const Tp3_& r3);

      template <real_prec_t prec1_, real_rnd_t rnd1_, class Tp_>
      friend typename enable_if<
        type_traits<typename result_type2<real<prec1_, rnd1_>, Tp_>::type,
          Tp_>::enable_math_funcs &&
        type_traits<typename result_type2<real<prec1_, rnd1_>, Tp_>::type,
          real<prec1_, rnd1_> >::enable_math_funcs,
        typename result_type2<real<prec1_, rnd1_>, Tp_>::type>::type
      modf(const Tp_& r, real<prec1_, rnd1_>* iptr);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      nextafter(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      nexttoward(const Tp1_& r1, const Tp2_& r2);

      template <class Tp1_, class Tp2_>
      friend typename enable_if<
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp1_>::enable_math_funcs &&
        type_traits<typename result_type2<Tp1_, Tp2_>::type,
          Tp2_>::enable_math_funcs,
        typename result_type2<Tp1_, Tp2_>::type>::type
      remquo(const Tp1_& r1, const Tp2_& r2, long* quo);

      // //////////////////////////////////////////////////////////////
      // mathematical constants
      // //////////////////////////////////////////////////////////////

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      const_log2() {
        Tp_ temp;
        mpfr_const_log2(temp.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      const_pi() {
        Tp_ temp;
        mpfr_const_pi(temp.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      const_euler() {
        Tp_ temp;
        mpfr_const_euler(temp.x_, rnd_);
        return temp;
      }

      template <class Tp_>
      friend inline typename enable_if<
        type_traits<real, Tp_>::enable_math_funcs &&
        equal_types2<real, Tp_>::val,
        const Tp_>::type
      const_catalan() {
        Tp_ temp;
        mpfr_const_catalan(temp.x_, rnd_);
        return temp;
      }
  };  // class real

}  // namespace mpfr

#endif  // MPFR_REAL_HPP
