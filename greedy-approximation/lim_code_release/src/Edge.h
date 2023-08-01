#include <string>
#include <set>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

class Edge{
public:
	Edge(int src, int dst, double value){
		this->src = src;
		this->dst = dst;
		this->value = value;
	}
	int src;//cited paper
	int dst;
	double value;
	bool operator < (const Edge& edge) const{
		if(edge.src!=this->src){
			return this->src < edge.src;
		}
		else{
			return this->dst < edge.dst;
		}
	}

	bool operator > (const Edge& edge) const{
		if(edge.src!=this->src){
			return this->src > edge.src;
		}
		else{
			return this->dst > edge.dst;
		}
	}

	bool operator != (const Edge& edge) const{
		return this->src != edge.src || this->dst != edge.dst;
	}

	bool operator == (const Edge& edge) const{
		return this->src == edge.src && this->dst == edge.dst;
	}

	void print() const{
		cout<<src<<"  "<<dst<<"  "<<value<<endl;
	}
};

class EdgeValue{
public:
	EdgeValue(int src, int dst, double value){
		this->src = src;
		this->dst = dst;
		this->value = value;
	}
	int src;//cited paper
	int dst;
	double value;
	bool operator < (const EdgeValue& edge) const{
		if(edge.value != this->value){
			return this->value < edge.value;
		}
		else if(edge.src!=this->src){
			return this->src < edge.src;
		}
		else{
			return this->dst < edge.dst;
		}
	}

	bool operator > (const EdgeValue& edge) const{
		if(edge.value != this->value){
			return this->value > edge.value;
		}
		else if(edge.src!=this->src){
			return this->src > edge.src;
		}
		else{
			return this->dst > edge.dst;
		}
	}

	bool operator != (const EdgeValue& edge) const{
		return this->src != edge.src || this->dst != edge.dst;
	}

	bool operator == (const EdgeValue& edge) const{
		return this->src == edge.src && this->dst == edge.dst;
	}

	void print() const{
		cout<<src<<"  "<<dst<<"  "<<value<<endl;
	}
};