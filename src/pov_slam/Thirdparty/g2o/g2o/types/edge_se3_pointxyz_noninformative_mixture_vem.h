// Max-Mixture plugin for g2o
// Copyright (C) 2012 P. Agarwal, E. Olson, W. Burgard
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

/*
 * Pratik Agarwal -- Max Mixture, the key thing here is that the graph can connect
 * multiple pair of nodes with only the required number of pairs being used in the optimization
 * The implementation becomes a little complicated when we try to model multiple data association
 * using a single edge
 */

#ifndef EDGE_SE3_POINTXYZ_NONINFORMATIVE_MIXTURE_VEM
#define EDGE_SE3_POINTXYZ_NONINFORMATIVE_MIXTURE_VEM

#include <boost/math/special_functions/digamma.hpp>

#include "edge_se3_pointxyz_noninformative.h"
#include "vertex_se3.h"
#include "vertex_pointxyz.h"
#include "g2o/config.h"

using namespace std;
using namespace Eigen;
using namespace g2o;

class EdgeSE3PointXYZNonInformativeMixtureELBO: public g2o::EdgeSE3PointXYZNonInformative {
public:

	EdgeSE3PointXYZNonInformativeMixtureELBO();
	EdgeSE3PointXYZNonInformativeMixtureELBO(std::pair<g2o::EdgeSE3PointXYZ*, g2o::EdgeSE3PointXYZNonInformative*> &_edges,
			const std::vector<double> &_weights, const std::vector<double> &beta_params);
	virtual bool read(std::istream &is);
	virtual bool write(std::ostream &os) const;
	virtual ~EdgeSE3PointXYZNonInformativeMixtureELBO();

	void initializeComponents(std::pair<g2o::EdgeSE3PointXYZ*, g2o::EdgeSE3PointXYZNonInformative*> &edges,
			const std::vector<double> &weights, const std::vector<double> &beta_params);
	void UpdateBelief(int i);
	void computeError();
	void linearizeOplus();

	double getEv() {
		return weights.at(0);
	}

	double getEpi() {
		return ELBOweights.at(0);
	}

	void computeBestEdge();
	int getBestComponent() const;

	const std::vector<double>& getELBOWeights() {
		return ELBOweights;
	}

	int numberComponents;

	EdgeSE3PointXYZ* getGaussianComponent();
	EdgeSE3PointXYZNonInformative* getNonInformativeComponent();

private:

	void getELBOParams();
	double getInlierLogProb();
	double getOutlierLogProb();
	double getNegLogProb(unsigned int c);

	//out of the components which is the best one (max-probability)

	int bestComponent;
	std::pair<g2o::EdgeSE3PointXYZ*, g2o::EdgeSE3PointXYZNonInformative*> allEdges;
	std::vector<double> weights;
	std::vector<double> ELBOweights;
	std::vector<double> betaParams;
	std::vector<double> determinants;
	bool verticesChanged;
	//to read vertices for subedges
	std::vector<std::pair<int, int> > vertexPairs;

	//bool onlyCorruptedGaussian = true; //if only corruptedGaussian then we dont need to allocate hessian memory
};

#endif
