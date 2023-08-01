#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include "Edge.h"
#include <vector>
#include <set>
#include <string>
using namespace std;

class SparseMatrix{
public:
	SparseMatrix(int size);

	~SparseMatrix();

	void insert(Edge edge);
	// insert this edge into two set with two different order.
	void erase(Edge edge);

	set<Edge>::iterator find(Edge edge);
	// return the iterator which points this edge.

	set<Edge>::iterator getIterFRowBegin(int n);

	set<Edge>::iterator getIterFColBegin(int n);

	set<Edge>::iterator getIterFRowEnd(int n);

	set<Edge>::iterator getIterFColEnd(int n);	

	SparseMatrix* operator * (double c);

	SparseMatrix* operator * (SparseMatrix* matrix);

	SparseMatrix* operator + (SparseMatrix* matrix);

	void sparsifyByValue(double minValue);

	int getSize();

	int getNonZeroNum();

	void print();

	void sparsifyByRatio(double beta);

	void rowNormalize();

	void saveToLocal(string name, string address = "./", map<int, int>* realId = NULL);

	bool ifExist(Edge edge);

	int outDegree(int node_id);

	int inDegree(int node_id);

private:
	vector<set<Edge>*>* rowCol;
	vector<set<Edge>*>* colRow;
	int size;
	int nonZeroNum;
};

#endif