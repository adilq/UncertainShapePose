// g2o - General Graph Optimization
// Copyright (C) 2011 H. Strasdat
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef G2O_MY_TYPES_H_
#define G2O_MY_TYPES_H_

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace g2o {
  using namespace Eigen;
  
  using number_t = double;
  
  inline constexpr number_t cst(long double v) { return (number_t)v; }

  typedef Matrix<number_t, 3, 1, ColMajor> Vector3;
  typedef Matrix<number_t, 4, 1, ColMajor> Vector4;
  typedef Matrix<number_t, 6, 1, ColMajor> Vector6;
  typedef Matrix<number_t, 7, 1, ColMajor> Vector7;
  
  typedef Matrix<number_t,3,3,ColMajor> Matrix3;
  
  typedef Transform<number_t,2,Isometry,ColMajor> Isometry2;
  typedef Transform<number_t,3,Isometry,ColMajor> Isometry3;
  
  typedef Eigen::Quaternion<number_t> Quaternion;
} // end namespace

#endif
