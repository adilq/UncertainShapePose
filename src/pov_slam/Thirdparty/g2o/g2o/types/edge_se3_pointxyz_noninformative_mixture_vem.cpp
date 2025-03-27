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
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE

#include <limits.h>
#include <math.h>

#include <g2o/core/factory.h>

#include "edge_se3_pointxyz_noninformative_mixture_vem.h"

EdgeSE3PointXYZNonInformativeMixtureELBO::EdgeSE3PointXYZNonInformativeMixtureELBO() :
		g2o::EdgeSE3PointXYZNonInformative::EdgeSE3PointXYZNonInformative() {
	numberComponents = 0;
	bestComponent = -1;
	ELBOweights = std::vector<double> {0.0, 0.0};
}

EdgeSE3PointXYZNonInformativeMixtureELBO::~EdgeSE3PointXYZNonInformativeMixtureELBO() {
	delete allEdges.first;
	delete allEdges.second;
}

EdgeSE3PointXYZNonInformativeMixtureELBO::EdgeSE3PointXYZNonInformativeMixtureELBO(
		std::pair<g2o::EdgeSE3PointXYZ*, g2o::EdgeSE3PointXYZNonInformative*> &_edges,
		const std::vector<double> &_weights, const std::vector<double> &_beta_params) :
		g2o::EdgeSE3PointXYZNonInformative::EdgeSE3PointXYZNonInformative() {
	initializeComponents(_edges, _weights, _beta_params);
}

void EdgeSE3PointXYZNonInformativeMixtureELBO::initializeComponents(
		std::pair<g2o::EdgeSE3PointXYZ*, g2o::EdgeSE3PointXYZNonInformative*> &_edges,
		const std::vector<double> &_weights,
		const std::vector<double> &_beta_params) {
	this->allEdges = _edges;
	this->weights = _weights;
	this->betaParams = _beta_params;
	numberComponents = 2;
	determinants.push_back(allEdges.first->information().inverse().determinant());
	determinants.push_back(allEdges.second->information().inverse().determinant());
	bestComponent = 0;
	UpdateBelief(bestComponent);
}

void EdgeSE3PointXYZNonInformativeMixtureELBO::getELBOParams() {
	double log_pi_tilda = boost::math::digamma(betaParams.at(0)) - boost::math::digamma(betaParams.at(0) + betaParams.at(1));
	double log_one_minus_pi_tilda = boost::math::digamma(betaParams.at(1)) - boost::math::digamma(betaParams.at(0) + betaParams.at(1));

	double p_in = exp(log_pi_tilda + getInlierLogProb());
	double p_out = exp(log_one_minus_pi_tilda + getOutlierLogProb());

	p_in = p_in / (p_in + p_out);
	p_out = p_out / (p_in + p_out);

	ELBOweights = std::vector<double> {p_in, p_out};
}

void EdgeSE3PointXYZNonInformativeMixtureELBO::UpdateBelief(int i) {
	//required for multimodal max-mixtures, this part may be optimized for speed,
	//do this only when component changed.
	if (i == 0) {
		this->setVertex(0, allEdges.first->vertex(0));
		this->setVertex(1, allEdges.first->vertex(1));

		double p[3];
		allEdges.first->getMeasurementData(p);
		this->setMeasurementData(p);
		this->setInformation(allEdges.first->information());
	} else {
		this->setVertex(0, allEdges.second->vertex(0));
		this->setVertex(1, allEdges.second->vertex(1));

		double p[3];
		allEdges.second->getMeasurementData(p);
		this->setMeasurementData(p);
		this->setInformation(allEdges.second->information());
	}
}

//get the edge with max probability and then create initialize the hessian
//memory if the nodes have changed
void EdgeSE3PointXYZNonInformativeMixtureELBO::computeError() {
	int best = -1;
	double minError = numeric_limits<double>::max();

	allEdges.first->setCaches(getCaches());
	allEdges.second->setCaches(getCaches());

	if (ELBOweights.at(0) == 0.0) getELBOParams();

	double thisNegLogProb = getNegLogProb(0);
	if (minError > thisNegLogProb) {
		best = 0;
		minError = thisNegLogProb;
	}

	thisNegLogProb = getNegLogProb(1);
	if (minError > thisNegLogProb) {
		best = 1;
		minError = thisNegLogProb;
	}

	bestComponent = best;
	UpdateBelief(bestComponent);
	if (bestComponent == 0) g2o::EdgeSE3PointXYZNonInformative::computeError(ELBOweights.at(0));
	else g2o::EdgeSE3PointXYZNonInformative::computeErrorNonInformative();
}


double EdgeSE3PointXYZNonInformativeMixtureELBO::getInlierLogProb() {
	allEdges.first->computeError();
	auto error = allEdges.first->error();
	return -1.5 * log(2 * M_PI) - 0.5 * log(determinants[0])
				- 0.5 * error.transpose() * allEdges.first->information() * error;
}

double EdgeSE3PointXYZNonInformativeMixtureELBO::getOutlierLogProb() {
	return log(allEdges.second -> getDensity());
}

double EdgeSE3PointXYZNonInformativeMixtureELBO::getNegLogProb(unsigned int c) {
	if (c == 0) {
		allEdges.first->computeError();
		auto error = allEdges.first->error();
		return -log(weights.at(0)) - getInlierLogProb();
	} else {
		return -log(weights.at(1)) - getOutlierLogProb();
	}
}

void EdgeSE3PointXYZNonInformativeMixtureELBO::linearizeOplus() {
	EdgeSE3PointXYZNonInformativeMixtureELBO::computeError();
	if (bestComponent == 0) {
		g2o::EdgeSE3PointXYZNonInformative::linearizeOplus();
	} else {
		g2o::EdgeSE3PointXYZNonInformative::linearizeOplusNonInformative();
	}
}

int EdgeSE3PointXYZNonInformativeMixtureELBO::getBestComponent() const {
	return bestComponent;
}

void EdgeSE3PointXYZNonInformativeMixtureELBO::computeBestEdge() {
	this->computeError();
}

//note g2o takes care of creating the first two vertices
//its problematic since all edges are taken care in this way
//EDGE_SE3_MIXTURE na nb numComponents Edgetype_i w_i na_i nb_i
bool EdgeSE3PointXYZNonInformativeMixtureELBO::read(std::istream &is) {
	int pId;
	is >> pId;
	setParameterId(0, pId);

	is >> bestComponent;

	is >> numberComponents;

	weights.reserve(numberComponents);
	determinants.reserve(numberComponents);

	double p1[3], p2[3];
	double w1, w2;

	//might throw error if first vertex is landmark
	VertexSE3 *va = static_cast<VertexSE3*>(this->vertex(0));
	assert(va != NULL);


		EdgeSE3PointXYZ *e1 = new EdgeSE3PointXYZ;
		EdgeSE3PointXYZNonInformative *e2 = new EdgeSE3PointXYZNonInformative;
		allEdges = std::make_pair(e1, e2);
		
		std::string buf;
		is >> buf;
		is >> w1;
		weights.push_back(w1);
		int na1, nb1;
		is >> na1;
		is >> nb1;

		VertexSE3 *v0 = static_cast<VertexSE3*>(va->graph()->vertex(na1));
		VertexPointXYZ *v1 = static_cast<VertexPointXYZ*>(va->graph()->vertex(nb1));
		assert(v0 != NULL);
		assert(v1 != NULL);
		allEdges.first->setVertex(0, v0);
		allEdges.first->setVertex(1, v1);
		allEdges.first->setParameterId(0, pId);

		is >> p1[0] >> p1[1] >> p1[2];
		allEdges.first->setMeasurementData(p1);

		InformationType inf1;
		for (int i = 0; i < information().rows(); ++i)
			for (int j = i; j < information().cols(); ++j) {
				is >> inf1(i, j);
				if (i != j)
					inf1(j, i) = inf1(i, j);
			}
		allEdges.first->setInformation(inf1);
		
		is >> buf;
		is >> w2;
		weights.push_back(w2);
		int na2, nb2;
		is >> na2;
		is >> nb2;

		v0 = static_cast<VertexSE3*>(va->graph()->vertex(na2));
		v1 = static_cast<VertexPointXYZ*>(va->graph()->vertex(nb2));
		assert(v0 != NULL);
		assert(v1 != NULL);
		allEdges.second->setVertex(0, v0);
		allEdges.second->setVertex(1, v1);
		allEdges.second->setParameterId(0, pId);

		is >> p2[0] >> p2[1] >> p2[2];
		allEdges.second->setMeasurementData(p2);

		InformationType inf2;
		for (int i = 0; i < information().rows(); ++i)
			for (int j = i; j < information().cols(); ++j) {
				is >> inf2(i, j);
				if (i != j)
					inf2(j, i) = inf2(i, j);
			}
		allEdges.second->setInformation(inf2);
	

		determinants.push_back(allEdges.first->information().inverse().determinant());
		determinants.push_back(allEdges.second->information().inverse().determinant());

	UpdateBelief(bestComponent);
	return is.good();
}

EdgeSE3PointXYZ* EdgeSE3PointXYZNonInformativeMixtureELBO::getGaussianComponent() { return allEdges.first; }

EdgeSE3PointXYZNonInformative* EdgeSE3PointXYZNonInformativeMixtureELBO::getNonInformativeComponent() { return allEdges.second; }

bool EdgeSE3PointXYZNonInformativeMixtureELBO::write(std::ostream &os) const {
	os << offsetParam->id() << " ";
	os << bestComponent << " ";
	os << numberComponents << " ";
	g2o::Factory *factory = g2o::Factory::instance();

	double p1[3], p2[3];
	

		os << factory->tag(allEdges.first) << " ";
		os << weights[0] << " ";
		os << allEdges.first->vertex(0)->id() << " ";
		os << allEdges.first->vertex(1)->id() << " ";

		allEdges.first->getMeasurementData(p1);
		os << p1[0] << " " << p1[1] << " " << p1[2] << " ";
		for (int i = 0; i < information().rows(); ++i)
			for (int j = i; j < information().cols(); ++j)
				os << " " << allEdges.first->information()(i, j);

		os << " ";
		
		os << factory->tag(allEdges.second) << " ";
		os << weights[1] << " ";
		os << allEdges.second->vertex(0)->id() << " ";
		os << allEdges.second->vertex(1)->id() << " ";

		allEdges.second->getMeasurementData(p2);
		os << p2[0] << " " << p2[1] << " " << p2[2] << " ";
		for (int i = 0; i < information().rows(); ++i)
			for (int j = i; j < information().cols(); ++j)
				os << " " << allEdges.second->information()(i, j);

		os << " ";

	return os.good();
}

