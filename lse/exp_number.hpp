/*****************************************************************************
*
* Copyright (C) 1997-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef LSE_EXP_NUMBER_HPP
#define LSE_EXP_NUMBER_HPP

#include <boost/throw_exception.hpp>
#include <boost/cstdint.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef ALPS_INDEP_SOURCE
# include <alps/osiris/dump.h>
#endif

namespace lse {

namespace detail {
    
template<typename T, typename U>
struct fp_promotion_traits;

#define FP_PROMOTION_TRAITS(T, U, V) \
template<> \
struct fp_promotion_traits<T, U> { \
  typedef V type; \
};

FP_PROMOTION_TRAITS(float, float, float)
FP_PROMOTION_TRAITS(float, double, double)
FP_PROMOTION_TRAITS(float, long double, long double)
FP_PROMOTION_TRAITS(double, float, double)
FP_PROMOTION_TRAITS(double, double, double)
FP_PROMOTION_TRAITS(double, long double, long double)
FP_PROMOTION_TRAITS(long double, float, long double)
FP_PROMOTION_TRAITS(long double, double, long double)
FP_PROMOTION_TRAITS(long double, long double, long double)

#undef FP_PROMOTION_TRAITS

} // end namespace detail

template<typename T>
class exp_number {
private:
  static const int32_t is_zero = 0;
  static const int32_t is_positive = 1;
  static const int32_t is_negative = -1;

public:
  typedef T value_type;
  typedef exp_number<value_type> self_;

  template<typename U> friend class exp_number;

  exp_number() : log_(0), sign_(is_zero) {}
  exp_number(self_ const& v) : log_(v.log_), sign_(v.sign_) {}
  exp_number(value_type v) {
    if (v > 0) {
      sign_ = is_positive;
      log_ = std::log(v);
    } else if (v < 0) {
      sign_ = is_negative;
      log_ = std::log(-v);
    } else {
      sign_ = is_zero;
      log_ = 0;
    }
  }
  template<typename U>
  exp_number(exp_number<U> const& v) : log_(static_cast<value_type>(v.log_)), sign_(v.sign_) {}

  void set_log(value_type v) {
    log_ = v;
    sign_ = is_positive;
  }

  operator value_type() const {
    if (sign_ == is_positive)
      return std::exp(log_);
    else if (sign_ == is_negative)
      return -std::exp(log_);
    else
      return 0;
  }
  value_type log() const {
    if (sign_ != is_positive) {
      std::cerr << (*this) << ' ' << sign_ << ' ' << log_ << std::endl;
      boost::throw_exception(std::range_error("exp_number::log()"));
    }
    return log_;
  }
  self_ operator-() const {
    self_ res(*this);
    res.sign_ *= is_negative;
    return res;
  }
  self_ pow(value_type p) const {
    self_ res(*this);
    res.log_ *= p;
    return res;
  }
  self_ sqrt() const {
    if (sign_ == is_negative)
      boost::throw_exception(std::range_error("exp_number::sqrt()"));
    self_ res(*this);
    res.log_ *= 0.5;
    return res;
  }

  bool operator>(int v) const { return this->operator>(static_cast<value_type>(v)); }
  bool operator>(value_type v) const {
    if (v >= 0)
      return (sign_ == is_positive) ? (log_ > std::log(v)) : false;
    else
      return (sign_ == is_negative) ? (log_ < std::log(-v)) : true;
  }
  template<typename U>
  bool operator>(exp_number<U> const& rhs) const {
    if (rhs.sign_ == is_positive)
      return (sign_ == is_positive) ? (log_ > rhs.log_) : false;
    else if (rhs.sign_ == is_negative)
      return (sign_ == is_negative) ? (log_ < rhs.log_) : true;
    else
      return (sign_ == is_positive);
  }

  bool operator==(int v) const { return this->operator==(static_cast<value_type>(v)); }
  bool operator==(value_type v) const {
    if (v > 0)
      return (sign_ == is_positive) ? (log_ == std::log(v)) : false;
    else if (v < 0)
      return (sign_ == is_negative) ? (log_ == std::log(-v)) : false;
    else
      return (sign_ == is_zero);
  }
  template<typename U>
  bool operator==(exp_number<U> const& rhs) const {
    if (rhs.sign_ == is_positive)
      return (sign_ == is_positive) ? (log_ == rhs.log_) : false;
    else if (rhs.sign_ == is_negative)
      return (sign_ == is_negative) ? (log_ == rhs.log_) : false;
    else
      return (sign_ == is_zero);
  }

  self_& operator+=(int v) { return operator+=(static_cast<value_type>(v)); }
  self_& operator+=(value_type v) { return operator+=(self_(v)); }
  template<typename U>
  self_& operator+=(exp_number<U> const& rhs) {
    if (sign_ == is_positive) {
      if (rhs.sign_ == is_positive) {
        // pos + pos
        log_ = std::max(log_, rhs.log_) + std::log(1 + std::exp(-std::abs(log_ - rhs.log_)));
      } else if (rhs.sign_ == is_negative) {
        if (log_ > rhs.log_) {
          // pos + neg = pos
          log_ = log_ + std::log(1 - std::exp(rhs.log_ - log_));
        } else if (log_ < rhs.log_) {
          // pos + neg = neg
          log_ = rhs.log_ + std::log(1 - std::exp(log_ - rhs.log_));
          sign_ = is_negative;
        } else {
          // pos + neg = zero
          sign_ = is_zero;
        }
      }
    } else if (sign_ == is_negative) {
      if (rhs.sign_ == is_positive) {
        if (log_ > rhs.log_) {
          // neg + pos = neg
          log_ = log_ + std::log(1 - std::exp(rhs.log_ - log_));
        } else if (log_ < rhs.log_) {
          // neg + pos = pos
          log_ = rhs.log_ + std::log(1 - std::exp(log_ - rhs.log_));
          sign_ = is_positive;
        } else {
          // neg + pos = zero
          sign_ = is_zero;
        }
      } else if (rhs.sign_ == is_negative)
        // neg + neg
        log_ = std::max(log_, rhs.log_) + std::log(1 + std::exp(-std::abs(log_ - rhs.log_)));
    } else {
      log_ = rhs.log_;
      sign_ = rhs.sign_;
    }
    return *this;
  }

  self_& operator-=(int v) { return this->operator+=(-v); }
  self_& operator-=(value_type v) { return this->operator+=(-v); }
  template<typename U>
  self_& operator-=(exp_number<U> const& rhs) { return this->operator+=(-rhs); }

  self_& operator*=(int v) { return this->operator*=(self_(v)); }
  self_& operator*=(value_type v) { return this->operator*=(self_(v)); }
  self_& operator*=(self_ const& rhs) {
    log_ += rhs.log_;
    sign_ *= rhs.sign_;
    if (std::abs(log_/rhs.log_) < 1.0e-10) log_ = 0;
    return *this;
  }

  self_& operator/=(int v) { return this->operator/=(self_(v)); }
  self_& operator/=(value_type v) { return this->operator/=(self_(v)); }
  self_& operator/=(self_ const& rhs) {
    if (sign_ != is_zero) {
      if (rhs.sign_ != is_zero) {
        log_ -= rhs.log_;
        sign_ *= rhs.sign_;
        if (std::abs(log_/rhs.log_) < 1.0e-10) log_ = 0;
      } else {
        boost::throw_exception(std::range_error("exp_number::operator/=()"));
      }
    }
    return *this;
  }

#ifndef ALPS_INDEP_SOURCE  
  alps::ODump& save(alps::ODump& dp) const { dp << log_ << sign_; return dp; }
  alps::IDump& load(alps::IDump& dp) { dp >> log_ >> sign_; return dp; }
#endif

private:
  value_type log_;
  int32_t sign_;
};

typedef exp_number<float> exp_float;
typedef exp_number<double> exp_double;
typedef exp_number<long double> exp_long_double;

template<typename T>
T log(exp_number<T> const& x) {
  return x.log();
}

template<typename T>
exp_number<T> exp(exp_number<T> const& x) {
  exp_number<T> res;
  res.set_log(static_cast<double>(x));
  return res;
}

template<typename T>
exp_number<T> pow(exp_number<T> const& x, int p) {
  return x.pow(p);
}

template<typename T>
exp_number<T> pow(exp_number<T> const& x, double p) {
  return x.pow(p);
}

template<typename T>
exp_number<T> sqrt(exp_number<T> const& x) {
  return x.sqrt();
}

inline exp_number<double> exp_value(double v) {
  exp_number<double> x;
  x.set_log(v);
  return x;
}

} // end namespace exp_number

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace lse {
#endif

//
// opertor>
//

template<typename T>
bool operator>(T x, lse::exp_number<T> const& y) {
  return (lse::exp_number<T>(x) > y);
}

//
// opertor==
//

template<typename T, typename U>
bool operator==(T x, lse::exp_number<T> const& y) {
  return (lse::exp_number<T>(x) == y);
}

//
// opertor>=
//

template<typename T, typename U>
bool operator>=(lse::exp_number<T> const& x, lse::exp_number<U> const& y) {
  return (x > y) || (x == y);
}
template<typename T, typename U>
bool operator>=(lse::exp_number<T> const& x, U y) {
  return (x > y) || (x == y);
}
template<typename T, typename U>
bool operator>=(U x, lse::exp_number<T> const& y) {
  lse::exp_number<T> lhs(x);
  return (lhs > y) || (lhs == y);
}

//
// opertor<
//

template<typename T, typename U>
bool operator<(lse::exp_number<T> const& x, lse::exp_number<U> const& y) {
  return (y > x);
}
template<typename T, typename U>
bool operator<(lse::exp_number<T> const& x, U y) {
  return (y > x);
}
template<typename T, typename U>
bool operator<(U x, lse::exp_number<T> const& y) {
  return (y > x);
}

//
// opertor<=
//

template<typename T, typename U>
bool operator<=(lse::exp_number<T> const& x, lse::exp_number<U> const& y) {
  return (y >= x);
}
template<typename T, typename U>
bool operator<=(lse::exp_number<T> const& x, U y) {
  lse::exp_number<T> z = y;
  return (z >= x);
}
template<typename T, typename U>
bool operator<=(U x, lse::exp_number<T> const& y) {
  return (y >= x);
}

//
// opertor+
//

template<typename T, typename U>
lse::exp_number<typename lse::detail::fp_promotion_traits<T, U>::type>
operator+(lse::exp_number<T> const& x, lse::exp_number<U> const& y) {
  typedef typename lse::detail::fp_promotion_traits<T, U>::type value_type;
  lse::exp_number<value_type> res = x;
  res += y;
  return res;
}
template<typename T, typename U>
lse::exp_number<T> operator+(lse::exp_number<T> const& x, U y) {
  lse::exp_number<T> res = x;
  res += y;
  return res;
}
template<typename T, typename U>
lse::exp_number<T> operator+(U x, lse::exp_number<T> const& y) {
  lse::exp_number<T> res = x;
  res += y;
  return res;
}

//
// opertor-
//

template<typename T, typename U>
lse::exp_number<typename lse::detail::fp_promotion_traits<T, U>::type>
operator-(lse::exp_number<T> const& x, lse::exp_number<U> const& y) {
  typedef typename lse::detail::fp_promotion_traits<T, U>::type value_type;
  lse::exp_number<value_type> res = x;
  res -= y;
  return res;
}
template<typename T, typename U>
lse::exp_number<T> operator-(lse::exp_number<T> const& x, U y) {
  lse::exp_number<T> res = x;
  res -= y;
  return res;
}
template<typename T, typename U>
lse::exp_number<T> operator-(U x, lse::exp_number<T> const& y) {
  lse::exp_number<T> res = x;
  res -= y;
  return res;
}

//
// operator*
//

template<typename T, typename U>
lse::exp_number<typename lse::detail::fp_promotion_traits<T, U>::type>
operator*(lse::exp_number<T> const& x, lse::exp_number<U> const& y) {
  typedef typename lse::detail::fp_promotion_traits<T, U>::type value_type;
  lse::exp_number<value_type> res = x;
  res *= y;
  return res;
}
template<typename T, typename U>
lse::exp_number<T> operator*(lse::exp_number<T> const& x, U y) {
  lse::exp_number<T> res = x;
  res *= y;
  return res;
}
template<typename T, typename U>
lse::exp_number<T> operator*(U x, lse::exp_number<T> const& y) {
  lse::exp_number<T> res = y;
  res *= x;
  return res;
}

//
// opertor/
//

template<typename T, typename U>
lse::exp_number<typename lse::detail::fp_promotion_traits<T, U>::type>
operator/(lse::exp_number<T> const& x, lse::exp_number<U> const& y) {
  typedef typename lse::detail::fp_promotion_traits<T, U>::type value_type;
  lse::exp_number<value_type> res = x;
  res /= y;
  return res;
}
template<typename T, typename U>
lse::exp_number<T> operator/(lse::exp_number<T> const& x, U y) {
  lse::exp_number<T> res = x;
  res /= y;
  return res;
}
template<typename T, typename U>
lse::exp_number<T> operator/(U x, lse::exp_number<T> const& y) {
  lse::exp_number<T> res = x;
  res /= y;
  return res;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, lse::exp_number<T> const& x) {
  if (x > 0)
    os << "exp(" << log(x) << ')';
  else if (x < 0)
    os << "-exp(" << log(-x) << ')';
  else
    os << "0";
  return os;
}

#ifndef ALPS_INDEP_SOURCE

template<typename T>
alps::ODump& operator<<(alps::ODump& dp, lse::exp_number<T> const& x) {
  x.save(dp);
  return dp;
}

template<typename T>
alps::IDump& operator>>(alps::IDump& dp, lse::exp_number<T>& x) {
  x.load(dp);
  return dp;
}

#endif

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // namespace lse
#endif

#endif // LSE_EXP_NUMBER_H
