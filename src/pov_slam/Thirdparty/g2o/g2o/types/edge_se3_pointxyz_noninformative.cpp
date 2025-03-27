// g2o - General Graph Optimization
// Copyright (C) 2011 R. Kuemmerle, G. Grisetti, W. Burgard
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

#include "edge_se3_pointxyz_noninformative.h"

#include <iostream>


namespace g2o {
  using namespace std;

  // point to camera projection, monocular
  EdgeSE3PointXYZNonInformative::EdgeSE3PointXYZNonInformative() : EdgeSE3PointXYZ() {}


  void EdgeSE3PointXYZNonInformative::computeError(double weight) {
	 VertexPointXYZ *point = static_cast<VertexPointXYZ*>(_vertices[1]);
	 Vector3 perr = cache->w2n() * point->estimate();
     _error = sqrt(weight) * (perr - _measurement);
  }

  void EdgeSE3PointXYZNonInformative::computeErrorNonInformative() {
    _error = Vector3::Zero();
  }


  void EdgeSE3PointXYZNonInformative::linearizeOplus() {
	VertexPointXYZ *vp = static_cast<VertexPointXYZ*>(_vertices[1]);

	Vector3 Zcam = cache->w2l() * vp->estimate();

	J(0, 4) = -2 * Zcam(2);
	J(0, 5) = 2 * Zcam(1);
	J(1, 3) = 2 * Zcam(2);
	J(1, 5) = -2 * Zcam(0);
	J(2, 3) = -2 * Zcam(1);
	J(2, 4) = 2 * Zcam(0);

	J.block<3, 3>(0, 6) = cache->w2l().rotation();

	Eigen::Matrix < number_t, 3, 9, Eigen::ColMajor > Jhom =
			offsetParam->inverseOffset().rotation() * J;

	_jacobianOplusXi = Jhom.block<3, 6>(0, 0);
	_jacobianOplusXj = Jhom.block<3, 3>(0, 6);
  }

  void EdgeSE3PointXYZNonInformative::linearizeOplusNonInformative() {
    _jacobianOplusXi = Eigen::Matrix<number_t,3,6,Eigen::ColMajor>::Zero();
    _jacobianOplusXj = Eigen::Matrix<number_t,3,3,Eigen::ColMajor>::Zero();
  }

} // end namespace
