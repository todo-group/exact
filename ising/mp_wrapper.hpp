/*****************************************************************************
*
* Copyright (C) 2015-2021 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Wrapper for multi-precision libraries

#pragma once

#include <limits>
#include <string>
#include <boost/math/policies/policy.hpp>

using std::enable_if;
using std::is_arithmetic;

template<class T>
struct mp_wrapper {
  typedef T base_type;
  mp_wrapper() {}
  mp_wrapper(const base_type& t) : base(t) {}
  template<class U>
  mp_wrapper(const U& u, typename enable_if<is_arithmetic<U>::value>::type* = 0) : base(u) {}
  operator int() const { return static_cast<int>(base); }
  template<class U>
  mp_wrapper& assign(const U& u) {
    base.assign(u);
    return *this;
  }
  base_type base;
};

template<class T>
mp_wrapper<T> operator-(const mp_wrapper<T>& t) { return mp_wrapper<T>(-t.base); }

template<class T>
mp_wrapper<T>& operator+=(mp_wrapper<T>& t, const mp_wrapper<T>& u) {
  t.base += u.base;
  return t;
}
template<class T, class V>
typename enable_if<is_arithmetic<V>::value, mp_wrapper<T>&>::type operator+=(mp_wrapper<T>& t, V v) {
  t.base += v;
  return t;
}

template<class T>
mp_wrapper<T>& operator-=(mp_wrapper<T>& t, const mp_wrapper<T>& u) {
  t.base -= u.base;
  return t;
}
template<class T, class V>
typename enable_if<is_arithmetic<V>::value, mp_wrapper<T>&>::type operator-=(mp_wrapper<T>& t, V v) {
  t.base -= v;
  return t;
}

template<class T>
mp_wrapper<T>& operator*=(mp_wrapper<T>& t, const mp_wrapper<T>& u) {
  t.base *= u.base;
  return t;
}
template<class T, class V>
typename enable_if<is_arithmetic<V>::value, mp_wrapper<T>&>::type operator*=(mp_wrapper<T>& t, V v) {
  t.base *= v;
  return t;
}

template<class T>
mp_wrapper<T>& operator/=(mp_wrapper<T>& t, const mp_wrapper<T>& u) {
  t.base /= u.base;
  return t;
}
template<class T, class V>
typename enable_if<is_arithmetic<V>::value, mp_wrapper<T>&>::type operator/=(mp_wrapper<T>& t, V v) {
  t.base /= v;
  return t;
}

template<class T>
mp_wrapper<T> operator+(const mp_wrapper<T>& t, const mp_wrapper<T>& u) {
  return mp_wrapper<T>(t.base + u.base);
}
template<class T, class V>
typename enable_if<is_arithmetic<V>::value, mp_wrapper<T>>::type operator+(const mp_wrapper<T>& t, const V& v) {
  return mp_wrapper<T>(t.base + v);
}
template<class T, class V>
typename enable_if<is_arithmetic<V>::value, mp_wrapper<T>>::type operator+(const V& v, const mp_wrapper<T>& t) {
  return mp_wrapper<T>(v + t.base);
}

template<class T>
mp_wrapper<T> operator-(const mp_wrapper<T>& t, const mp_wrapper<T>& u) {
  return mp_wrapper<T>(t.base - u.base);
}
template<class T, class V>
typename enable_if<is_arithmetic<V>::value, mp_wrapper<T>>::type operator-(const mp_wrapper<T>& t, V v) {
  return mp_wrapper<T>(t.base - v);
}
template<class T, class V>
typename enable_if<is_arithmetic<V>::value, mp_wrapper<T>>::type operator-(V v, const mp_wrapper<T>& t) {
  return mp_wrapper<T>(v - t.base);
}

template<class T>
mp_wrapper<T> operator*(const mp_wrapper<T>& t, const mp_wrapper<T>& u) {
  return mp_wrapper<T>(t.base * u.base);
}
template<class T, class V>
typename enable_if<is_arithmetic<V>::value, mp_wrapper<T>>::type operator*(const mp_wrapper<T>& t, V v) {
  return mp_wrapper<T>(t.base * v);
}
template<class T, class V>
typename enable_if<is_arithmetic<V>::value, mp_wrapper<T>>::type operator*(V v, const mp_wrapper<T>& t) {
  return mp_wrapper<T>(v * t.base);
}

template<class T>
mp_wrapper<T> operator/(const mp_wrapper<T>& t, const mp_wrapper<T>& u) {
  return mp_wrapper<T>(t.base / u.base);
}
template<class T, class V>
mp_wrapper<T> operator/(const mp_wrapper<T>& t, V v) { return mp_wrapper<T>(t.base / v); }
template<class T, class V>
mp_wrapper<T> operator/(V v, const mp_wrapper<T>& t) { return mp_wrapper<T>(v / t.base); }

template<class T>
bool operator==(const mp_wrapper<T>& t, const mp_wrapper<T>& u) { return t.base == u.base; }
template<class T, class V>
bool operator==(const mp_wrapper<T>& t, const V& v) { return t.base == v; }
template<class T, class V>
bool operator==(const V& v, const mp_wrapper<T>& t) { return v == t.base; }

template<class T>
bool operator!=(const mp_wrapper<T>& t, const mp_wrapper<T>& u) { return t.base != u.base; }
template<class T, class V>
bool operator!=(const mp_wrapper<T>& t, const V& v) { return t.base != v; }
template<class T, class V>
bool operator!=(const V& v, const mp_wrapper<T>& t) { return v != t.base; }

template<class T>
bool operator<(const mp_wrapper<T>& t, const mp_wrapper<T>& u) { return t.base < u.base; }
template<class T, class V>
bool operator<(const mp_wrapper<T>& t, const V& v) { return t.base < v; }
template<class T, class V>
bool operator<(const V& v, const mp_wrapper<T>& t) { return v < t.base; }

template<class T>
bool operator>(const mp_wrapper<T>& t, const mp_wrapper<T>& u) { return t.base > u.base; }
template<class T, class V>
bool operator>(const mp_wrapper<T>& t, const V& v) { return t.base > v; }
template<class T, class V>
bool operator>(const V& v, const mp_wrapper<T>& t) { return v > t.base; }

template<class T>
bool operator<=(const mp_wrapper<T>& t, const mp_wrapper<T>& u) { return t.base <= u.base; }
template<class T, class V>
bool operator<=(const mp_wrapper<T>& t, const V& v) { return t.base <= v; }
template<class T, class V>
bool operator<=(const V& v, const mp_wrapper<T>& t) { return v <= t.base; }

template<class T>
bool operator>=(const mp_wrapper<T>& t, const mp_wrapper<T>& u) { return t.base >= u.base; }
template<class T, class V>
bool operator>=(const mp_wrapper<T>& t, const V& v) { return t.base >= v; }
template<class T, class V>
bool operator>=(const V& v, const mp_wrapper<T>& t) { return v >= t.base; }

// mathematical functions

template<class T>
mp_wrapper<T> abs(const mp_wrapper<T>& t) { return mp_wrapper<T>(abs(t.base)); }

template<class T>
mp_wrapper<T> acos(const mp_wrapper<T>& t) { return mp_wrapper<T>(acos(t.base)); }

template<class T>
mp_wrapper<T> asin(const mp_wrapper<T>& t) { return mp_wrapper<T>(asin(t.base)); }

template<class T>
mp_wrapper<T> atan(const mp_wrapper<T>& t) { return mp_wrapper<T>(atan(t.base)); }

template<class T>
mp_wrapper<T> atan2(const mp_wrapper<T>& t) { return mp_wrapper<T>(atan2(t.base)); }

template<class T>
mp_wrapper<T> ceil(const mp_wrapper<T>& t) { return mp_wrapper<T>(ceil(t.base)); }

template<class T>
mp_wrapper<T> cos(const mp_wrapper<T>& t) { return mp_wrapper<T>(cos(t.base)); }

template<class T>
mp_wrapper<T> cosh(const mp_wrapper<T>& t) { return mp_wrapper<T>(cosh(t.base)); }

template<class T>
mp_wrapper<T> exp(const mp_wrapper<T>& t) { return mp_wrapper<T>(exp(t.base)); }

template<class T>
mp_wrapper<T> fabs(const mp_wrapper<T>& t) { return mp_wrapper<T>(fabs(t.base)); }

template<class T>
mp_wrapper<T> floor(const mp_wrapper<T>& t) { return mp_wrapper<T>(floor(t.base)); }

template<class T>
mp_wrapper<T> fmod(const mp_wrapper<T>& t, const mp_wrapper<T>& u) {
  return mp_wrapper<T>(fmod(t.base, u.base));
}

template<class T>
mp_wrapper<T> frexp(const mp_wrapper<T>& t, int* u) { return mp_wrapper<T>(frexp(t.base, u)); }

template<class T>
mp_wrapper<T> ldexp(const mp_wrapper<T>& t, int u) { return mp_wrapper<T>(ldexp(t.base, u)); }

template<class T>
mp_wrapper<T> log(const mp_wrapper<T>& t) { return mp_wrapper<T>(log(t.base)); }

template<class T>
mp_wrapper<T> log10(const mp_wrapper<T>& t) { return mp_wrapper<T>(log10(t.base)); }

template<class T>
mp_wrapper<T> pow(const mp_wrapper<T>& t, const mp_wrapper<T>& u) {
  return mp_wrapper<T>(pow(t.base, u.base));
}

template<class T>
mp_wrapper<T> round(const mp_wrapper<T>& t) { return mp_wrapper<T>(round(t.base)); }

template<class T>
mp_wrapper<T> sin(const mp_wrapper<T>& t) { return mp_wrapper<T>(sin(t.base)); }

template<class T>
mp_wrapper<T> sinh(const mp_wrapper<T>& t) { return mp_wrapper<T>(sinh(t.base)); }

template<class T>
mp_wrapper<T> sqrt(const mp_wrapper<T>& t) { return mp_wrapper<T>(sqrt(t.base)); }

template<class T>
mp_wrapper<T> tan(const mp_wrapper<T>& t) { return mp_wrapper<T>(tan(t.base)); }

template<class T>
mp_wrapper<T> tanh(const mp_wrapper<T>& t) { return mp_wrapper<T>(tanh(t.base)); }

template<class T>
mp_wrapper<T> trunc(const mp_wrapper<T>& t) { return mp_wrapper<T>(trunc(t.base)); }

template<class T>
inline std::istream& operator>>(std::istream& is, mp_wrapper<T>& t) { return is >> t.base; }

template<class T>
inline std::ostream& operator<<(std::ostream& os, const mp_wrapper<T>& t) { return os << t.base; }

namespace std {

template<class T>
class numeric_limits<mp_wrapper<T>> {
public:
  static constexpr bool is_specialized = numeric_limits<T>::is_specialized;

  static constexpr mp_wrapper<T> epsilon() { return numeric_limits<T>::epsilon(); }
  static constexpr mp_wrapper<T> infinity() { return numeric_limits<T>::infinity(); }
  static constexpr mp_wrapper<T> lowest() { return numeric_limits<T>::lowest(); }
  static constexpr mp_wrapper<T> max() { return numeric_limits<T>::max(); }
  static constexpr mp_wrapper<T> min() { return numeric_limits<T>::min(); }
  static constexpr mp_wrapper<T> quiet_NaN() { return numeric_limits<T>::quiet_NaN(); }
  static constexpr mp_wrapper<T> round_error() { return numeric_limits<T>::round_error(); }

  static constexpr int digits = numeric_limits<T>::digits;
  static constexpr int digits10 = numeric_limits<T>::digits10;
  static constexpr int max_digits10 = numeric_limits<T>::max_digits10;
  static constexpr int max_exponent = numeric_limits<T>::max_exponent;
  static constexpr int max_exponent10 = numeric_limits<T>::max_exponent10;
  static constexpr int min_exponent = numeric_limits<T>::min_exponent;
  static constexpr int min_exponent10 = numeric_limits<T>::min_exponent10;
  static constexpr int radix = numeric_limits<T>::radix;

  static constexpr bool is_exact = numeric_limits<T>::is_exact;
  static constexpr bool is_integer = numeric_limits<T>::is_integer;
  static constexpr bool is_signed = numeric_limits<T>::is_signed;

  static constexpr bool has_denorm = numeric_limits<T>::has_denorm;
  static constexpr bool has_denorm_loss = numeric_limits<T>::has_denorm_loss;
  static constexpr bool has_infinity = numeric_limits<T>::has_infinity;
  static constexpr bool has_quiet_NaN = numeric_limits<T>::has_quiet_NaN;
  static constexpr bool has_signaling_NaN = numeric_limits<T>::has_signaling_NaN;
};

} // end namespace std

namespace boost {
namespace math {
namespace policies {

template <class T, class Policy>
struct precision<mp_wrapper<T>, Policy> {
  typedef typename Policy::precision_type precision_type;
  typedef typename precision<T, Policy>::digits_2 digits_2;
  typedef typename precision<T, Policy>::type type;
};

} // end namespace policies
} // end namespace math
} // end namespace boost

// convert from string

template<class T>
inline T convert(const char* str) {
  T r;
  r.assign(str);
  return r;
}

template<class T>
inline T convert(const std::string& str) {
  T r;
  r.assign(str);
  return r;
}

template<>
inline float convert<float>(const char* str) { return std::stof(str); }

template<>
inline float convert<float>(const std::string& str) { return std::stof(str); }

template<>
inline double convert<double>(const char* str) { return std::stod(str); }

template<>
inline double convert<double>(const std::string& str) { return std::stod(str); }

template<>
inline long double convert<long double>(const char* str) { return std::stold(str); }

template<>
inline long double convert<long double>(const std::string& str) { return std::stold(str); }
