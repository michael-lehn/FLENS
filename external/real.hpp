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
        throw(): _msg(msg) {}
      virtual ~exception_real() throw() {}
      // returns cause of error
      virtual const char* what() const throw() { return _msg.c_str(); };
    private:
      std::string _msg;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  class real;

  // //////////////////////////////////////////////////////////////////
  // type traits
  // //////////////////////////////////////////////////////////////////

  template <class _Tp1, class _Tp2>
  struct type_traits {
    typedef _Tp1 real_type;
    typedef _Tp2 other_type;

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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, unsigned long int> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, long int> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, unsigned int> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, int> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, unsigned short int> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, short int> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, unsigned char> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, signed char> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, char> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, wchar_t> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, float> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, double> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, long double> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpz_t> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpq_t> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpf_t> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpfr_t> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpz_ptr> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpq_ptr> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpf_ptr> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpfr_ptr> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpz_srcptr> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpq_srcptr> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpf_srcptr> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, mpfr_srcptr> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, std::string> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, char*> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd>
  struct type_traits<real<_prec, _rnd>, const char*> {
    typedef real<_prec, _rnd> real_type;
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

  template <real_prec_t _prec, real_rnd_t _rnd, int _i>
  struct type_traits<real<_prec, _rnd>, char[_i]> {
    typedef real<_prec, _rnd> real_type;
    typedef char other_type[_i];

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

  template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd>
  struct type_traits<real<_prec1, _rnd>, real<_prec2, _rnd> > {
    typedef real<_prec1, _rnd> real_type;
    typedef real<_prec2, _rnd> other_type;

    // support level in class real
    static const bool enable_impl_ctor  = true;
    static const bool enable_assign_op  = true;
    static const bool enable_conv_func  = true;
    static const bool enable_arithm_ops = true;
    static const bool enable_compar_ops = true;
    static const bool enable_math_funcs = true;

    // support in MPFR library
    // (all "has_..." booleans *MUST* be set to "false" for "real<_prec, _rnd>")
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
    // (there should be no function definitions for "real<_prec, _rnd>")
  };

}  // namespace mpfr



namespace mpfr {

  // //////////////////////////////////////////////////////////////
  // declaration of real
  // //////////////////////////////////////////////////////////////

  template <real_prec_t _prec, real_rnd_t _rnd>
  class real;

  // //////////////////////////////////////////////////////////////////
  // basic meta-programming
  // //////////////////////////////////////////////////////////////////

  // enable_if

  template<bool, class _Tp = void>
  struct enable_if {};

  template<class _Tp>
  struct enable_if<true, _Tp> {
    typedef _Tp type;
  };

  // //////////////////////////////////////////////////////////////////
  // calculated result type for two arguments
  // //////////////////////////////////////////////////////////////////

  // At least one argument must be of type "real" and all "real"s must have
  // the same rounding ("_rnd").  The result type is the "real" with highest
  // precision.

  template<class _Tp1, class _Tp2>
  struct result_type2 {
  };

  template <real_prec_t _prec, real_rnd_t _rnd, class _Tp>
  struct result_type2<real<_prec, _rnd>, _Tp> {
    typedef real<_prec, _rnd> type;
    static const real_rnd_t  rnd  = _rnd;
    static const real_prec_t prec = _prec;
  };

  template <real_prec_t _prec, real_rnd_t _rnd, class _Tp>
  struct result_type2<_Tp, real<_prec, _rnd> > {
    typedef real<_prec, _rnd> type;
    static const real_rnd_t  rnd  = _rnd;
    static const real_prec_t prec = _prec;
  };

  template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd>
  struct result_type2<real<_prec1, _rnd>, real<_prec2, _rnd> > {
    typedef real<((_prec1 < _prec2) ? _prec2 : _prec1), _rnd> type;
    static const real_rnd_t  rnd  = _rnd;
    static const real_prec_t prec = ((_prec1 < _prec2) ? _prec2 : _prec1);
  };

  template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd1, real_rnd_t _rnd2>
  struct result_type2<real<_prec1, _rnd1>, real<_prec2, _rnd2> > {
  };

  // //////////////////////////////////////////////////////////////////
  // calculated result type for three arguments
  // //////////////////////////////////////////////////////////////////

  // At least one argument must be of type "real" and all "real"s must have
  // the same rounding ("_rnd").  The result type is the "real" with highest
  // precision.

  template<class _Tp1, class _Tp2, class _Tp3>
  struct result_type3 {
  };

  template<real_prec_t _prec, real_rnd_t _rnd, class _Tp1, class _Tp2>
  struct result_type3<real<_prec, _rnd>, _Tp1, _Tp2> {
    typedef real<_prec, _rnd> type;
    static const real_rnd_t  rnd  = _rnd;
    static const real_prec_t prec = _prec;
  };

  template<real_prec_t _prec, real_rnd_t _rnd, class _Tp1, class _Tp2>
  struct result_type3<_Tp2, real<_prec, _rnd>, _Tp1> {
    typedef real<_prec, _rnd> type;
    static const real_rnd_t  rnd  = _rnd;
    static const real_prec_t prec = _prec;
  };

  template<real_prec_t _prec, real_rnd_t _rnd, class _Tp1, class _Tp2>
  struct result_type3<_Tp1, _Tp2, real<_prec, _rnd> > {
    typedef real<_prec, _rnd> type;
    static const real_rnd_t  rnd  = _rnd;
    static const real_prec_t prec = _prec;
  };

  template<real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd, class _Tp>
  struct result_type3<real<_prec1, _rnd>, real<_prec2, _rnd>, _Tp> {
    typedef real<((_prec1 < _prec2) ? _prec2 : _prec1), _rnd> type;
    static const real_rnd_t  rnd  = _rnd;
    static const real_prec_t prec = ((_prec1 < _prec2) ? _prec2 : _prec1);
  };

  template<real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd1, real_rnd_t _rnd2, class _Tp>
  struct result_type3<real<_prec1, _rnd1>, real<_prec2, _rnd2>, _Tp> {
  };

  template<real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd, class _Tp>
  struct result_type3<_Tp, real<_prec1, _rnd>, real<_prec2, _rnd> > {
    typedef real<((_prec1 < _prec2) ? _prec2 : _prec1), _rnd> type;
    static const real_rnd_t  rnd  = _rnd;
    static const real_prec_t prec = ((_prec1 < _prec2) ? _prec2 : _prec1);
  };

  template<real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd1, real_rnd_t _rnd2, class _Tp>
  struct result_type3<_Tp, real<_prec1, _rnd1>, real<_prec2, _rnd2> > {
  };

  template<real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd, class _Tp>
  struct result_type3<real<_prec2, _rnd>, _Tp, real<_prec1, _rnd> > {
    typedef real<((_prec1 < _prec2) ? _prec2 : _prec1), _rnd> type;
    static const real_rnd_t  rnd  = _rnd;
    static const real_prec_t prec = ((_prec1 < _prec2) ? _prec2 : _prec1);
  };

  template<real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd1, real_rnd_t _rnd2, class _Tp>
  struct result_type3<real<_prec2, _rnd2>, _Tp, real<_prec1, _rnd1> > {
  };

  template<real_prec_t _prec1, real_prec_t _prec2, real_prec_t _prec3, real_rnd_t _rnd>
  struct result_type3<real<_prec1, _rnd>, real<_prec2, _rnd>, real<_prec3, _rnd> > {
    typedef real<((((_prec3 < _prec2) ? _prec2 : _prec3) < _prec1) ? _prec1 :
                    ((_prec3 < _prec2) ? _prec2 : _prec3)), _rnd> type;
    static const real_rnd_t  rnd  = _rnd;
    static const real_prec_t prec = ((((_prec3 < _prec2) ? _prec2 : _prec3)
      < _prec1) ? _prec1 : ((_prec3 < _prec2) ? _prec2 : _prec3));
  };

  template<real_prec_t _prec1, real_prec_t _prec2, real_prec_t _prec3,
    real_rnd_t _rnd1, real_rnd_t _rnd2, real_rnd_t _rnd3>
  struct result_type3<real<_prec1, _rnd1>, real<_prec2, _rnd2>, real<_prec3, _rnd3> > {
  };

  // //////////////////////////////////////////////////////////////////
  // promotion to real
  // //////////////////////////////////////////////////////////////////

  template <class _Tp1, class _Tp2>
  struct promote {
  };

  template <real_prec_t _prec, real_rnd_t _rnd, class _Tp>
  struct promote<real<_prec, _rnd>, _Tp> {
    typedef const real<_prec, _rnd> type;
  };

  template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd>
  struct promote<real<_prec1, _rnd>, real<_prec2, _rnd> > {
    typedef const real<_prec2, _rnd>& type;
  };

  // //////////////////////////////////////////////////////////////////
  // check for equal types
  // //////////////////////////////////////////////////////////////////

  template <class _Tp1, class _Tp2>
  struct equal_types2 {
    static const bool val = false;
  };

  template <class _Tp>
  struct equal_types2<_Tp, _Tp> {
    static const bool val = true;
  };

  // //////////////////////////////////////////////////////////////////
  // check for type "real"
  // //////////////////////////////////////////////////////////////////

  template <class _Tp1, class _Tp2>
  struct has_real2 {
    static const bool val = false;
  };

  template <real_prec_t _prec, real_rnd_t _rnd, class _Tp>
  struct has_real2<real<_prec, _rnd>, _Tp> {
    static const bool val = true;
  };

  template <real_prec_t _prec, real_rnd_t _rnd, class _Tp>
  struct has_real2<_Tp, real<_prec, _rnd> > {
    static const bool val = true;
  };

  template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd>
  struct has_real2<real<_prec1, _rnd>, real<_prec2, _rnd> > {
    static const bool val = true;
  };

  template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd1, real_rnd_t _rnd2>
  struct has_real2<real<_prec1, _rnd1>, real<_prec2, _rnd2> > {
    static const bool val = false;
  };

}  // namespace mpfr




namespace mpfr {

  // //////////////////////////////////////////////////////////////
  // class declaration
  // //////////////////////////////////////////////////////////////

  template <real_prec_t _prec, real_rnd_t _rnd>
  class real;

  // //////////////////////////////////////////////////////////////
  // generic operators (definitions of binary operators)
  // //////////////////////////////////////////////////////////////

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_arithm_ops &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_arithm_ops,
    const typename result_type2<_Tp1, _Tp2>::type>::type
  operator +(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_add(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_arithm_ops &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_arithm_ops,
    const typename result_type2<_Tp1, _Tp2>::type>::type
  operator -(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_sub(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_arithm_ops &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_arithm_ops,
    const typename result_type2<_Tp1, _Tp2>::type>::type
  operator *(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_mul(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_arithm_ops &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_arithm_ops,
    const typename result_type2<_Tp1, _Tp2>::type>::type
  operator /(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_div(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_compar_ops &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_compar_ops,
    const bool>::type
  operator ==(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_equal_p(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_compar_ops &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_compar_ops,
    const bool>::type
  operator !=(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_lessgreater_p(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_compar_ops &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_compar_ops,
    const bool>::type
  operator <(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_less_p(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_compar_ops &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_compar_ops,
    const bool>::type
  operator <=(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_lessequal_p(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_compar_ops &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_compar_ops,
    const bool>::type
  operator >(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_greater_p(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_compar_ops &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_compar_ops,
    const bool>::type
  operator >=(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_greaterequal_p(temp1._x, temp2._x);
  }

  // //////////////////////////////////////////////////////////////
  // mathematical functions (definitions for multiple "real" arguments)
  // //////////////////////////////////////////////////////////////

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    const int>::type
  isgreater(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_greater_p(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    const int>::type
  isgreaterequal(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_greaterequal_p(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    const int>::type
  isless(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_less_p(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    const int>::type
  islessequal(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_lessequal_p(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    const int>::type
  islessgreater(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_lessgreater_p(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    const int>::type
  isunordered(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_unordered_p(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  atan2(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_atan2(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  copysign(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_copysign(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  fdim(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_dim(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2, class _Tp3>
  inline typename enable_if<
    type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
      _Tp2>::enable_math_funcs &&
    type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
      _Tp3>::enable_math_funcs,
    typename result_type3<_Tp1, _Tp2, _Tp3>::type>::type
  fma(const _Tp1& r1, const _Tp2& r2, const _Tp3& r3) {
    typedef typename result_type3<_Tp1, _Tp2, _Tp3>::type temp_type;
    const real_rnd_t rnd = result_type3<_Tp1, _Tp2, _Tp3>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    typename promote<temp_type, _Tp3>::type temp3(r3);
    temp_type temp;
    mpfr_fma(temp._x, temp1._x, temp2._x, temp3._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  fmax(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_max(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  fmin(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_min(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  fmod(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_fmod(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  hypot(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_hypot(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  pow(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_pow(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  remainder(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_remainder(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  agm(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_agm(temp._x, temp1._x, temp2._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    const int>::type
  cmpabs(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    return mpfr_cmpabs(temp1._x, temp2._x);
  }

  template <class _Tp1, class _Tp2, class _Tp3>
  inline typename enable_if<
    type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
      _Tp2>::enable_math_funcs &&
    type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
      _Tp3>::enable_math_funcs,
    typename result_type3<_Tp1, _Tp2, _Tp3>::type>::type
  fms(const _Tp1& r1, const _Tp2& r2, const _Tp3& r3) {
    typedef typename result_type3<_Tp1, _Tp2, _Tp3>::type temp_type;
    const real_rnd_t rnd = result_type3<_Tp1, _Tp2, _Tp3>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    typename promote<temp_type, _Tp3>::type temp3(r3);
    temp_type temp;
    mpfr_fms(temp._x, temp1._x, temp2._x, temp3._x, rnd);
    return temp;
  }

  template <real_prec_t _prec, real_rnd_t _rnd, class _Tp>
  inline typename enable_if<
    type_traits<typename result_type2<real<_prec, _rnd>, _Tp>::type,
      _Tp>::enable_math_funcs &&
    type_traits<typename result_type2<real<_prec, _rnd>, _Tp>::type,
      real<_prec, _rnd> >::enable_math_funcs,
    typename result_type2<real<_prec, _rnd>, _Tp>::type>::type
  modf(const _Tp& r, real<_prec, _rnd>* iptr) {
    typedef typename result_type2<real<_prec, _rnd>, _Tp>::type temp_type;
    const real_rnd_t rnd = result_type2<real<_prec, _rnd>, _Tp>::rnd;
    temp_type temp;
    mpfr_modf(iptr->_x, temp._x, r._x, rnd);
    return temp;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  nextafter(const _Tp1& r1, const _Tp2& r2) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    temp_type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    mpfr_nexttoward(temp1._x, temp2._x);
    return temp1;
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  nexttoward(const _Tp1& r1, const _Tp2& r2) {
    return nextafter(r1, r2);
  }

  template <class _Tp1, class _Tp2>
  inline typename enable_if<
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp1>::enable_math_funcs &&
    type_traits<typename result_type2<_Tp1, _Tp2>::type,
      _Tp2>::enable_math_funcs,
    typename result_type2<_Tp1, _Tp2>::type>::type
  remquo(const _Tp1& r1, const _Tp2& r2, long* quo) {
    typedef typename result_type2<_Tp1, _Tp2>::type temp_type;
    const real_rnd_t rnd = result_type2<_Tp1, _Tp2>::rnd;
    typename promote<temp_type, _Tp1>::type temp1(r1);
    typename promote<temp_type, _Tp2>::type temp2(r2);
    temp_type temp;
    mpfr_remquo(temp._x, quo, temp1._x, temp2._x, rnd);
    return temp;
  }

  // //////////////////////////////////////////////////////////////////
  // class definition
  // //////////////////////////////////////////////////////////////////

  template <real_prec_t _prec = 53, real_rnd_t _rnd = MPFR_RNDN>
  class real {
    // private:
    public:
      mpfr_t _x;

    public:
      // //////////////////////////////////////////////////////////////
      // default and copy constructors, default assignment operator, destructor
      // //////////////////////////////////////////////////////////////

      // default and copy constructor

      inline real() {
        mpfr_init2(_x, _prec);
        mpfr_set_zero(_x, +1);
      }

      inline real(const real& o) {
        mpfr_init2(_x, _prec);
        mpfr_set(_x, o._x, _rnd);
      }

      // default assignment operator

      inline real& operator =(const real& o) {
        if (&o != this)
          mpfr_set(_x, o._x, _rnd);
        return *this;
      }

      // destructor

      inline ~real() {
        mpfr_clear(_x);
      }

      // //////////////////////////////////////////////////////////////
      // converting constructors and converting assignment operators
      // //////////////////////////////////////////////////////////////

      // friend of other reals

      template <real_prec_t _prec1, real_rnd_t _rnd1>
      friend class real;

      // implicit conversion constructors

      template <class _Tp>
      inline real(const _Tp& o,
          typename enable_if<type_traits<real, _Tp>::has_set &&
          type_traits<real, _Tp>::enable_impl_ctor>::type* dummy = 0) {
        mpfr_init2(_x, _prec);
        type_traits<real, _Tp>::set(_x, o, _rnd);
      }

      template <real_prec_t _prec1, real_rnd_t _rnd1>
      inline real(const real<_prec1, _rnd1>& o,
          typename enable_if<
          type_traits<real, real<_prec1, _rnd1> >::enable_impl_ctor>::type*
          dummy = 0) {
        mpfr_init2(_x, _prec1);
        mpfr_set(_x, o._x, _rnd1);
      }

      // explicit conversion constructors

      template <class _Tp>
      inline explicit real(const _Tp& o,
          typename enable_if<
          type_traits<real, _Tp>::has_set &&
          (! type_traits<real, _Tp>::enable_impl_ctor)>::type* dummy = 0) {
        mpfr_init2(_x, _prec);
        type_traits<real, _Tp>::set(_x, o, _rnd);
      }

      template <real_prec_t _prec1, real_rnd_t _rnd1>
      inline explicit real(const real<_prec1, _rnd1>& o,
          typename enable_if<
          (! type_traits<real, real<_prec1, _rnd1> >::enable_impl_ctor)>::type*
          dummy = 0) {
        mpfr_init2(_x, _prec1);
        mpfr_set(_x, o._x, _rnd1);
      }

      // converting assignment operators

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_assign_op &&
        type_traits<real, _Tp>::has_set,
        real&>::type
      operator =(const _Tp& o) {
        type_traits<real, _Tp>::set(_x, o, _rnd);
        return *this;
      }

      template <real_prec_t _prec1, real_rnd_t _rnd1>
      inline typename enable_if<
        type_traits<real, real<_prec1, _rnd1> >::enable_assign_op,
        real&>::type
      operator =(const real<_prec1, _rnd1>& o) {
        mpfr_set(_x, o._x, _rnd1);
        return *this;
      }

      // //////////////////////////////////////////////////////////////
      // generic operators
      // //////////////////////////////////////////////////////////////

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        (! type_traits<real, _Tp>::has_add),
        real&>::type
      operator +=(const _Tp& o) {
        typename promote<real, _Tp>::type temp(o);
        mpfr_add(_x, _x, temp._x, _rnd);
        return *this;
      }

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        (! type_traits<real, _Tp>::has_sub_a),
        real&>::type
      operator -=(const _Tp& o) {
        typename promote<real, _Tp>::type temp(o);
        mpfr_sub(_x, _x, temp._x, _rnd);
        return *this;
      }

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        (! type_traits<real, _Tp>::has_mul),
        real&>::type
      operator *=(const _Tp& o) {
        typename promote<real, _Tp>::type temp(o);
        mpfr_mul(_x, _x, temp._x, _rnd);
        return *this;
      }

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        (! type_traits<real, _Tp>::has_div_a),
        real&>::type
      operator /=(const _Tp& o) {
        typename promote<real, _Tp>::type temp(o);
        mpfr_div(_x, _x, temp._x, _rnd);
        return *this;
      }

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_arithm_ops &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_arithm_ops,
        const typename result_type2<_Tp1, _Tp2>::type>::type
      operator +(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_arithm_ops &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_arithm_ops,
        const typename result_type2<_Tp1, _Tp2>::type>::type
      operator -(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_arithm_ops &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_arithm_ops,
        const typename result_type2<_Tp1, _Tp2>::type>::type
      operator *(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_arithm_ops &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_arithm_ops,
        const typename result_type2<_Tp1, _Tp2>::type>::type
      operator /(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_compar_ops &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_compar_ops,
        const bool>::type
      operator ==(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_compar_ops &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_compar_ops,
        const bool>::type
      operator !=(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_compar_ops &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_compar_ops,
        const bool>::type
      operator <(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_compar_ops &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_compar_ops,
        const bool>::type
      operator <=(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_compar_ops &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_compar_ops,
        const bool>::type
      operator >(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_compar_ops &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_compar_ops,
        const bool>::type
      operator >=(const _Tp1& r1, const _Tp2& r2);

      // //////////////////////////////////////////////////////////////
      // optimized operators
      // //////////////////////////////////////////////////////////////

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_add,
        real&>::type
      operator +=(const _Tp& o) {
        type_traits<real, _Tp>::add(_x, _x, o, _rnd);
        return *this;
      }

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_sub_a,
        real&>::type
      operator -=(const _Tp& o) {
        type_traits<real, _Tp>::sub_a(_x, _x, o, _rnd);
        return *this;
      }

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_mul,
        real&>::type
      operator *=(const _Tp& o) {
        type_traits<real, _Tp>::mul(_x, _x, o, _rnd);
        return *this;
      }

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_div_a,
        real&>::type
      operator /=(const _Tp& o) {
        type_traits<real, _Tp>::div_a(_x, _x, o, _rnd);
        return *this;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_add,
        const real>::type
      operator +(const real& r1, const _Tp& r2) {
        real temp;
        type_traits<real, _Tp>::add(temp._x, r1._x, r2, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_add,
        const real>::type
      operator +(const _Tp& r1, const real& r2) {
        real temp;
        type_traits<real, _Tp>::add(temp._x, r2._x, r1, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_sub_a,
        const real>::type
      operator -(const real& r1, const _Tp& r2) {
        real temp;
        type_traits<real, _Tp>::sub_a(temp._x, r1._x, r2, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_sub_b,
        const real>::type
      operator -(const _Tp& r1, const real& r2) {
        real temp;
        type_traits<real, _Tp>::sub_b(temp._x, r1, r2._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_mul,
        const real>::type
      operator *(const real& r1, const _Tp& r2) {
        real temp;
        type_traits<real, _Tp>::mul(temp._x, r1._x, r2, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_mul,
        const real>::type
      operator *(const _Tp& r1, const real& r2) {
        real temp;
        type_traits<real, _Tp>::mul(temp._x, r2._x, r1, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_div_a,
        const real>::type
      operator /(const real& r1, const _Tp& r2) {
        real temp;
        type_traits<real, _Tp>::div_a(temp._x, r1._x, r2, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_arithm_ops &&
        type_traits<real, _Tp>::has_div_b,
        const real>::type
      operator /(const _Tp& r1, const real& r2) {
        real temp;
        type_traits<real, _Tp>::div_b(temp._x, r1, r2._x, _rnd);
        return temp;
      }

      // //////////////////////////////////////////////////////////////
      // conversion operators and functions
      // //////////////////////////////////////////////////////////////

      // conversion operators
      // (can be enabled with preprocessor macro)

      #ifdef REAL_ENABLE_CONVERSION_OPERATORS
      inline operator unsigned long int() const {
        return mpfr_get_ui(_x, _rnd);
      }

      inline operator long int() const {
        return mpfr_get_si(_x, _rnd);
      }

      inline operator unsigned int() const {
        return mpfr_get_ui(_x, _rnd);
      }

      inline operator int() const {
        return mpfr_get_si(_x, _rnd);
      }

      inline operator unsigned short int() const {
        return mpfr_get_ui(_x, _rnd);
      }

      inline operator short int() const {
        return mpfr_get_si(_x, _rnd);
      }

      inline operator unsigned char() const {
        return mpfr_get_ui(_x, _rnd);
      }

      inline operator signed char() const {
        return mpfr_get_si(_x, _rnd);
      }

      inline operator char() const {
        return mpfr_get_si(_x, _rnd);
      }

      inline operator wchar_t() const {
        return mpfr_get_si(_x, _rnd);
      }

      inline operator float() const {
        return mpfr_get_flt(_x, _rnd);
      }

      inline operator double() const {
        return mpfr_get_d(_x, _rnd);
      }

      inline operator long double() const {
        return mpfr_get_ld(_x, _rnd);
      }

      inline operator std::string() const {
        std::stringstream temp;
        temp.precision(-1);
        try {
          temp << *this;
        }
        catch (...) {
          throw exception_real(
            "in real<_prec, _rnd>& real<_prec, _rnd>::operator std::string() const:\n  conversion failed");
        }
        return temp.str();
      }
      #endif  // REAL_ENABLE_CONVERSION_OPERATORS

      // conversion functions

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_conv_func &&
        type_traits<real, _Tp>::has_get_a,
        void>::type
      conv(_Tp& o) const {
        o = type_traits<real, _Tp>::get_a(_x, _rnd);
      }

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_conv_func &&
        (! type_traits<real, _Tp>::has_get_a) &&
        type_traits<real, _Tp>::has_get_b,
        void>::type
      conv(_Tp const o) const {
        type_traits<real, _Tp>::get_b(o, _x, _rnd);
      }

      template <class _Tp>
      inline typename enable_if<
        type_traits<real, _Tp>::enable_conv_func &&
        equal_types2<char*, _Tp>::val,
        void>::type
      conv(_Tp const o) const {
        std::stringstream temp;
        temp.precision(-1);
        try {
          temp << *this;
        }
        catch (...) {
          throw exception_real(
            "in const char* real<_prec, _rnd>::c_str() const:\n  conversion failed");
        }
        strcpy(o, temp.str().c_str());
      }

      // //////////////////////////////////////////////////////////////
      // increment, decrement, and negation operators
      // //////////////////////////////////////////////////////////////

      // increment operators

      inline real& operator ++() {
        static const real<_prec, _rnd> _one(1);
        *this += _one;
        return *this;
      }

      inline const real operator ++(int) {
        real<_prec, _rnd> temp = *this;
        ++(*this);
        return temp;
      }

      // decrement operators

      inline real& operator --() {
        static const real<_prec, _rnd> _one(1);
        *this -= _one;
        return *this;
      }

      inline const real operator --(int) {
        real<_prec, _rnd> temp = *this;
        --(*this);
        return temp;
      }

      // NOTE: The unary member operator- is declared after any template
      // binary friend operator-, because the latter may be unqualified
      // in the code above.  This way we make sure that binary - operations
      // do not match the unary member operator- (in any case).

      inline const real operator -() const {
        real<_prec, _rnd> temp;
        mpfr_neg(temp._x, _x, _rnd);
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
            std::string("in std::istream& operator >>(std::istream& s, real<_prec, _rnd>& r):\n  invalid input format ")
            + t);
        }

        return s;
      }

      // there might be some room for improvements for the next
      // missing: handling of ios_base::fixed/ios_base::scientific and
      // ios_base::showpoint

      friend inline std::ostream& operator <<(std::ostream& s, const real& r) {
        real_exp_t exp;
        char* ch = mpfr_get_str(0, &exp, 10, s.precision() + 1, r._x, _rnd);
        if (! ch)
          throw exception_real(
            "in std::ostream& operator <<(std::ostream& s, const real<_prec, _rnd>& r):\n  conversion failed");
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
          if (! mpfr_zero_p(r._x)) {
            const real_exp_t exp_prev = exp;
            volatile real_exp_t* exp_ptr = &exp;
            exp--;
            if (*exp_ptr > exp_prev)
              throw exception_real(
                "in std::ostream& operator <<(std::ostream& s, const real<_prec, _rnd>& r):\n  exponent out of range");
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

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const int>::type
      isfinite(const _Tp& r) {
        return mpfr_number_p(r._x);
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const int>::type
      isinf(const _Tp& r) {
        return mpfr_inf_p(r._x);
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const int>::type
      isnan(const _Tp& r) {
        return mpfr_nan_p(r._x);
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const int>::type
      isnormal(const _Tp& r) {
        return mpfr_regular_p(r._x);
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const int>::type
      signbit(const _Tp& r) {
        return mpfr_signbit(r._x);
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      acos(const _Tp& r) {
        real temp;
        mpfr_acos(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      acosh(const _Tp& r) {
        real temp;
        mpfr_acosh(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      asin(const _Tp& r) {
        real temp;
        mpfr_asin(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      asinh(const _Tp& r) {
        real temp;
        mpfr_asinh(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      atan(const _Tp& r) {
        real temp;
        mpfr_atan(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      atanh(const _Tp& r) {
        real temp;
        mpfr_atanh(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      cbrt(const _Tp& r) {
        real temp;
        mpfr_cbrt(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      ceil(const _Tp& r) {
        real temp;
        mpfr_ceil(temp._x, r._x);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      cos(const _Tp& r) {
        real temp;
        mpfr_cos(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      cosh(const _Tp& r) {
        real temp;
        mpfr_cosh(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      erf(const _Tp& r) {
        real temp;
        mpfr_erf(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      erfc(const _Tp& r) {
        real temp;
        mpfr_erfc(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      exp(const _Tp& r) {
        real temp;
        mpfr_exp(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      exp2(const _Tp& r) {
        real temp;
        mpfr_exp2(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      expm1(const _Tp& r) {
        real temp;
        mpfr_expm1(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      fabs(const _Tp& r) {
        real temp;
        mpfr_abs(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      abs(const _Tp& r) {
        real temp;
        mpfr_abs(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      floor(const _Tp& r) {
        real temp;
        mpfr_floor(temp._x, r._x);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      log(const _Tp& r) {
        real temp;
        mpfr_log(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      log10(const _Tp& r) {
        real temp;
        mpfr_log10(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      log1p(const _Tp& r) {
        real temp;
        mpfr_log1p(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      log2(const _Tp& r) {
        real temp;
        mpfr_log2(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      nearbyint(const _Tp& r) {
        real temp;
        mpfr_rint(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      rint(const _Tp& r) {
        real temp;
        mpfr_rint(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      round(const _Tp& r) {
        real temp;
        mpfr_round(temp._x, r._x);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      sin(const _Tp& r) {
        real temp;
        mpfr_sin(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      sinh(const _Tp& r) {
        real temp;
        mpfr_sinh(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      sqrt(const _Tp& r) {
        real temp;
        mpfr_sqrt(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      tan(const _Tp& r) {
        real temp;
        mpfr_tan(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      tanh(const _Tp& r) {
        real temp;
        mpfr_tanh(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      tgamma(const _Tp& r) {
        real temp;
        mpfr_gamma(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      trunc(const _Tp& r) {
        real temp;
        mpfr_trunc(temp._x, r._x);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      j0(const _Tp& r) {
        real temp;
        mpfr_j0(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      j1(const _Tp& r) {
        real temp;
        mpfr_j1(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      y0(const _Tp& r) {
        real temp;
        mpfr_y0(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      y1(const _Tp& r) {
        real temp;
        mpfr_y1(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      ai(const _Tp& r) {
        real temp;
        mpfr_ai(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      cot(const _Tp& r) {
        real temp;
        mpfr_cot(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      coth(const _Tp& r) {
        real temp;
        mpfr_coth(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      csc(const _Tp& r) {
        real temp;
        mpfr_csc(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      csch(const _Tp& r) {
        real temp;
        mpfr_csch(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      digamma(const _Tp& r) {
        real temp;
        mpfr_digamma(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      exp10(const _Tp& r) {
        real temp;
        mpfr_exp10(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      expint(const _Tp& r) {
        real temp;
        mpfr_eint(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      frac(const _Tp& r) {
        real temp;
        mpfr_frac(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const int>::type
      isinteger(const _Tp& r) {
        return mpfr_integer_p(r._x);
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const int>::type
      iszero(const _Tp& r) {
        return mpfr_zero_p(r._x);
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      li2(const _Tp& r) {
        real temp;
        mpfr_li2(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      rec_sqrt(const _Tp& r) {
        real temp;
        mpfr_rec_sqrt(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      sec(const _Tp& r) {
        real temp;
        mpfr_sec(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      sech(const _Tp& r) {
        real temp;
        mpfr_sech(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const int>::type
      sgn(const _Tp& r) {
        return mpfr_sgn(r._x);
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real>::type
      zeta(const _Tp& r) {
        real temp;
        mpfr_zeta(temp._x, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const int>::type
      fpclassify(const _Tp& r) {
        if (mpfr_nan_p(r._x))
          return FP_NAN;
        else if (mpfr_inf_p(r._x))
          return FP_INFINITE;
        else if (mpfr_zero_p(r._x))
          return FP_ZERO;
        else
          return FP_NORMAL;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      frexp(const _Tp& r, real_exp_t* exp) {
        if (mpfr_zero_p(r._x)) {
          *exp = 0;
          return r;
        }
        else if (mpfr_inf_p(r._x) || mpfr_nan_p(r._x)) {
          //*exp = 0;
          return r;
        }
        else {
          _Tp temp = r;
          *exp = mpfr_get_exp(r._x);
          mpfr_set_exp(temp._x, 0);
          return temp;
        }
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const real_exp_t>::type
      ilogb(const _Tp& r) {
        if (mpfr_zero_p(r._x) || mpfr_nan_p(r._x))
          return std::numeric_limits<real_exp_t>::min();
        else if (mpfr_inf_p(r._x))
          return std::numeric_limits<real_exp_t>::max();
        else {
          real_exp_t temp = mpfr_get_exp(r._x);
          if (temp != std::numeric_limits<real_exp_t>::min())
            temp--;
          return temp;
        }
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp >::type
      ldexp(const _Tp& r, const long exp) {
        _Tp temp;
        mpfr_mul_2si(temp._x, r._x, exp, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      lgamma(const _Tp& r) {
        _Tp temp;
        int signp;
        mpfr_lgamma(temp._x, &signp, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      logb(const _Tp& r) {
        _Tp temp;
        if (mpfr_zero_p(r._x))
          mpfr_set_inf(temp._x, -1);
        else if (mpfr_nan_p(r._x))
          mpfr_set_nan(temp._x);
        else if (mpfr_inf_p(r._x))
          mpfr_set_inf(temp._x, 1);
        else {
          temp = mpfr_get_exp(r._x);
          temp--;
        }
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      nan(const char* tagp) {
        _Tp temp;
        mpfr_set_nan(temp._x);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      scalbln(const _Tp& r, const long exp) {
        return ldexp(r, exp);  // FLT_RADIX == 2???
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      scalbn(const _Tp& r, const int exp) {
        return ldexp(r, exp);  // FLT_RADIX == 2???
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      jn(const long n, const _Tp& r) {
        _Tp temp;
        mpfr_jn(temp._x, n, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      yn(const long n, const _Tp& r) {
        _Tp temp;
        mpfr_yn(temp._x, n, r._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      cyl_bessel_j(const long n, const _Tp& r) {
        return jn(n, r);
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      cyl_neumann(const long n, const _Tp& r) {
        return yn(n, r);
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      factorial(const long n) {
        _Tp temp;
        mpfr_fac_ui(temp._x, n, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      root(const _Tp& r, const unsigned long n) {
        _Tp temp;
        mpfr_root(temp._x, r._x, n, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      inf(const int n) {
        _Tp temp;
        mpfr_set_inf(temp._x, n);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      zero(const int n) {
        _Tp temp;
        mpfr_set_zero(temp._x, n);
        return temp;
      }

      // //////////////////////////////////////////////////////////////
      // mathematical functions with multiple "real" arguments
      // //////////////////////////////////////////////////////////////

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        const int>::type
      isgreater(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        const int>::type
      isgreaterequal(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        const int>::type
      isless(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        const int>::type
      islessequal(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        const int>::type
      islessgreater(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        const int>::type
      isunordered(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      atan2(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      copysign(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      fdim(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2, class _Tp3>
      friend typename enable_if<
        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
          _Tp2>::enable_math_funcs &&
        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
          _Tp3>::enable_math_funcs,
        typename result_type3<_Tp1, _Tp2, _Tp3>::type>::type
      fma(const _Tp1& r1, const _Tp2& r2, const _Tp3& r3);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      fmax(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      fmin(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      fmod(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      hypot(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      pow(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      remainder(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      agm(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        const int>::type
      cmpabs(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2, class _Tp3>
      friend typename enable_if<
        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
          _Tp2>::enable_math_funcs &&
        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3>::type,
          _Tp3>::enable_math_funcs,
        typename result_type3<_Tp1, _Tp2, _Tp3>::type>::type
      fms(const _Tp1& r1, const _Tp2& r2, const _Tp3& r3);

      template <real_prec_t _prec1, real_rnd_t _rnd1, class _Tp>
      friend typename enable_if<
        type_traits<typename result_type2<real<_prec1, _rnd1>, _Tp>::type,
          _Tp>::enable_math_funcs &&
        type_traits<typename result_type2<real<_prec1, _rnd1>, _Tp>::type,
          real<_prec1, _rnd1> >::enable_math_funcs,
        typename result_type2<real<_prec1, _rnd1>, _Tp>::type>::type
      modf(const _Tp& r, real<_prec1, _rnd1>* iptr);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      nextafter(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      nexttoward(const _Tp1& r1, const _Tp2& r2);

      template <class _Tp1, class _Tp2>
      friend typename enable_if<
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp1>::enable_math_funcs &&
        type_traits<typename result_type2<_Tp1, _Tp2>::type,
          _Tp2>::enable_math_funcs,
        typename result_type2<_Tp1, _Tp2>::type>::type
      remquo(const _Tp1& r1, const _Tp2& r2, long* quo);

      // //////////////////////////////////////////////////////////////
      // mathematical constants
      // //////////////////////////////////////////////////////////////

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      const_log2() {
        _Tp temp;
        mpfr_const_log2(temp._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      const_pi() {
        _Tp temp;
        mpfr_const_pi(temp._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      const_euler() {
        _Tp temp;
        mpfr_const_euler(temp._x, _rnd);
        return temp;
      }

      template <class _Tp>
      friend inline typename enable_if<
        type_traits<real, _Tp>::enable_math_funcs &&
        equal_types2<real, _Tp>::val,
        const _Tp>::type
      const_catalan() {
        _Tp temp;
        mpfr_const_catalan(temp._x, _rnd);
        return temp;
      }
  };  // class real

}  // namespace mpfr

//
// import math functions into std
//

namespace std {

using mpfr::pow;

} // namespace std


#endif  // MPFR_REAL_HPP
