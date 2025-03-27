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

#include <limits.h>
#include <math.h>

#include <g2o/core/factory.h>

#include "edge_se3_mixture.h"


EdgeSE3Mixture::EdgeSE3Mixture() :
		g2o::EdgeSE3Orig::EdgeSE3Orig() {
	numberComponents = 0;
	int bestComponent = -1;
}

EdgeSE3Mixture::~EdgeSE3Mixture() {
	for (int i(0); i < numberComponents; i++)
		delete allEdges[i];
}

EdgeSE3Mixture::EdgeSE3Mixture(std::vector<g2o::EdgeSE3Orig*> &_edges,
		std::vector<double> &_weights) :
		g2o::EdgeSE3Orig::EdgeSE3Orig() {
	initializeComponents(_edges, _weights);
}

void EdgeSE3Mixture::initializeComponents(std::vector<g2o::EdgeSE3Orig*> &_edges,
		std::vector<double> &_weights) {
	this->allEdges = _edges;
	this->weights = _weights;
	//UpdateBelief(0);
	numberComponents = _edges.size();
	for (unsigned int i = 0; i < numberComponents; i++)
		determinants.push_back(
				allEdges[i]->information().inverse().determinant());

	computeBestEdge();
}

void EdgeSE3Mixture::UpdateBelief(int i) {
	//required for multimodal max-mixtures, this part may be optimized for speed,
	//do this only when component changed.
	this->setVertex(0, allEdges[i]->vertex(0));
	this->setVertex(1, allEdges[i]->vertex(1));

	double p[7];
	allEdges[i]->getMeasurementData(p);
	this->setMeasurementData(p);
	this->setInformation(allEdges[i]->information());
}

void EdgeSE3Mixture::computeError() {
	int best = -1;
	double minError = numeric_limits<double>::max();
	for (unsigned int i = 0; i < numberComponents; i++) {
		//allEdges[i]->computeError();
		double thisNegLogProb = getNegLogProb(i);
		if (minError > thisNegLogProb) {
			best = i;
			minError = thisNegLogProb;
		}
	}

	//if(best!=bestComponent){
	bestComponent = best;
	UpdateBelief(bestComponent);
	//}
	g2o::EdgeSE3Orig::computeError();
}

double EdgeSE3Mixture::getNegLogProb(unsigned int c) {
	allEdges[c]->computeError();
	//cerr << "\nchi2 for " <<c<<" "<< allEdges[c]->chi2();
	return -log(weights[c]) + 0.5 * log(determinants[c])
			+ 0.5 * allEdges[c]->chi2();
}

void EdgeSE3Mixture::linearizeOplus() {
	EdgeSE3Mixture::computeError();
	g2o::EdgeSE3Orig::linearizeOplus();
}

int EdgeSE3Mixture::getBestComponent() const {
	return bestComponent;
}

void EdgeSE3Mixture::computeBestEdge() {
	computeError();
}

bool EdgeSE3Mixture::read(std::istream &is) {

	is >> numberComponents;

	allEdges.reserve(numberComponents);
	weights.reserve(numberComponents);
	determinants.reserve(numberComponents);

	Vector7 p;
	double w;

	VertexSE3 *va = static_cast<VertexSE3*>(this->vertex(0));

	for (int c = 0; c < numberComponents; c++) {
		EdgeSE3Orig *e = new EdgeSE3Orig;
		allEdges.push_back(e);
		std::string buf;
		is >> buf;
		is >> w;
		weights.push_back(w);
		int na, nb;
		is >> na;
		is >> nb;

		VertexSE3 *v0 = static_cast<VertexSE3*>(va->graph()->vertex(na));
		VertexSE3 *v1 = static_cast<VertexSE3*>(va->graph()->vertex(nb));
		assert(v0 != NULL);
		assert(v1 != NULL);
		allEdges[c]->setVertex(0, v0);
		allEdges[c]->setVertex(1, v1);

		for (int i = 0; i < 7; i++)
			is >> p[i];
		Vector4d::MapType(p.data() + 3).normalize();

		allEdges[c]->setMeasurementData(p.data());
		InformationType inf;
		for (int i = 0; i < 6; ++i)
			for (int j = i; j < 6; ++j) {
				is >> inf(i, j);
				if (i != j)
					inf(j, i) = inf(i, j);
			}
		allEdges[c]->setInformation(inf);
	}

	for (unsigned int i = 0; i < numberComponents; i++)
		determinants.push_back(
				allEdges[i]->information().inverse().determinant());

	computeBestEdge();
	return is.good();
}

bool EdgeSE3Mixture::write(std::ostream &os) const {
	os << numberComponents << " ";
	g2o::Factory *factory = g2o::Factory::instance();

	double p[7];
	//cerr <<"\n";
	for (int c = 0; c < numberComponents; c++) {
		os << factory->tag(allEdges[c]) << " ";
		os << weights[c] << " ";
		//cerr<< "Writing component "<< c << " ";
		os << allEdges[c]->vertex(0)->id() << " ";
		os << allEdges[c]->vertex(1)->id() << " ";
		//os >> this->_vertices[0]->id()<<" ";
		//os << this->_vertices[1]->id()<< " ";

		allEdges[c]->getMeasurementData(p);
		for (int i = 0; i < 7; i++)
			os << p[i] << " ";
		//os<<p[0] <<" " <<p[1]<<" "<< p[2] <<" ";
		InformationType inf;
		for (int i = 0; i < 6; ++i)
			for (int j = i; j < 6; ++j)
				os << " " << allEdges[c]->information()(i, j);

		os << " ";
	}
	return os.good();
}

EdgeSE3Orig* EdgeSE3Mixture::getComponent(int i) {
	return allEdges[i];
}

