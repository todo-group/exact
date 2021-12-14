/*
   Copyright (C) 2015-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

#include <string>

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
