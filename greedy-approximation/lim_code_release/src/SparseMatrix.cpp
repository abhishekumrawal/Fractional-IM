#include "SparseMatrix.h"
#include <iostream>
#include <fstream> 
#include <sstream>
#include <iostream>
//#include <set>

SparseMatrix::SparseMatrix(int size){
	this->size = size;
	nonZeroNum = 0;
	rowCol = new vector<set<Edge>*>(size);
	colRow = new vector<set<Edge>*>(size);
    	for(int i = 0; i<size; i++){
        	rowCol->at(i) = new set<Edge>;
        	colRow->at(i) = new set<Edge>;
    	}
}

SparseMatrix::~SparseMatrix(){
	for(int i=0; i<size; i++){
		delete rowCol->at(i);
		delete colRow->at(i);
	}
	delete rowCol;
	delete colRow;
}

void SparseMatrix::insert(Edge edge){
	if(rowCol->at(edge.src)->find(edge) != rowCol->at(edge.src)->end()) return;
	rowCol->at(edge.src)->insert(edge);
	colRow->at(edge.dst)->insert(edge);
	nonZeroNum++;
}

set<Edge>::iterator SparseMatrix::find(Edge edge){
	return rowCol->at(edge.src)->find(edge);
}

void SparseMatrix::erase(Edge edge){
	rowCol->at(edge.src)->erase(edge);
	colRow->at(edge.dst)->erase(edge);
	nonZeroNum--;
}

set<Edge>::iterator SparseMatrix::getIterFRowBegin(int n){
	if(n>= size){
		cout << "Use the outIndex of matrix"<<endl;
		exit(0);
	}
	return rowCol->at(n)->begin();
}

set<Edge>::iterator SparseMatrix::getIterFColBegin(int n){
	if(n>= size){
		cout << "Use the outIndex of matrix"<<endl;
		exit(0);
	}
	return colRow->at(n)->begin();
}

set<Edge>::iterator SparseMatrix::getIterFRowEnd(int n){
	if(n>= size){
		cout << "Use the outIndex of matrix"<<endl;
		exit(0);
	}
	return rowCol->at(n)->end();
}

set<Edge>::iterator SparseMatrix::getIterFColEnd(int n){
	if(n>= size){
		cout << "Use the outIndex of matrix"<<endl;
		exit(0);
	}
	return colRow->at(n)->end();
}

SparseMatrix* SparseMatrix::operator * (double c){
	SparseMatrix* result = new SparseMatrix(size);
	set<Edge>::iterator it;
	for(int i=0; i<size; i++){
		for(it = getIterFRowBegin(i); it != getIterFRowEnd(i); it++){
			result->insert(Edge((*it).src, (*it).dst, (*it).value * c));
		}
	}
	return result;
}

SparseMatrix* SparseMatrix::operator * (SparseMatrix* matrix){
	cout << "begin to multiply a matrix"<<endl;
	cout << "Two matrices' sizes are "<< this->getNonZeroNum() << " and " << matrix->getNonZeroNum() <<endl;
	SparseMatrix* result = new SparseMatrix(size);
	set<Edge>::iterator it, it1, it2;
	double value;
	double progressBar = 0.0;
	cout<<"this multiplication is processing "<<progressBar;
	for(int i=0; i<size; i++){
		if(double(i)/size > progressBar+0.1){
			progressBar += 0.1;
			cout<<"                                          \r";
			cout<<"this multiplication is processing "<<progressBar;
		}
		for(it1 = getIterFColBegin(i); it1 != getIterFColEnd(i); it1++){
			for(it2 = matrix->getIterFRowBegin(i); it2 != matrix->getIterFRowEnd(i); it2++){
				it = result->find(Edge(it1->src, it2->dst, 0));
				if(it != result->getIterFRowEnd(it1->src)){
					value = it->value + it1->value * it2->value;
					result->erase(Edge(it1->src, it2->dst, 0));
					result->insert(Edge(it1->src, it2->dst, value));
				}
				else{
					value = it1->value * it2->value;
					result->insert(Edge(it1->src, it2->dst, value));
				}
			}
		}
	}
	cout<<"                                     \r";
	cout<<"finish this multiplication!"<<endl;
	return result;
}

SparseMatrix* SparseMatrix::operator + (SparseMatrix* matrix){
	cout << "begin to add a matrix"<<endl;
	SparseMatrix* result = new SparseMatrix(size);
	set<Edge>::iterator it1, it2;
	double progressBar = 0.0;
	cout << "Two matrices' sizes are "<< this->getNonZeroNum() << " and " << matrix->getNonZeroNum() <<endl;
        cout<<"this multiplication is processing "<<progressBar;

	for(int i=0; i<size; i++){
		if(double(i)/size > progressBar+0.1){
                        progressBar += 0.1;
                        cout<<"                                          \r";
                        cout<<"this multiplication is processing "<<progressBar;
                }
		it1 = this->getIterFRowBegin(i);
		it2 = matrix->getIterFRowBegin(i);
		bool ifcontinue = true;
		while(ifcontinue){
			if(it1 != this->getIterFRowEnd(i) && it2 != matrix->getIterFRowEnd(i)){
				if((*it1) < (*it2)){
					result->insert((*it1));
					it1++;
				}
				else if((*it1) > (*it2)){
					result->insert((*it2));
					it2++;
				}
				else{
					result->insert(Edge(it1->src, it1->dst, it1->value + it2->value));
					it1++;
					it2++;
				}
			}
			else if(it1 == this->getIterFRowEnd(i)){
				while(it2 != matrix->getIterFRowEnd(i)){
					result->insert((*it2));
					it2++;
				}
				break;			
			}
			else{
				while(it1 != this->getIterFRowEnd(i)){
					result->insert((*it1));
					it1++;
				}
				break;				
			}
		}
	}
	cout<<"                                          \r";
	cout<<"finish this addition!"<<endl;
	return result;
}

int SparseMatrix::getSize(){
	return size;
}

int SparseMatrix::getNonZeroNum(){
	return nonZeroNum;
}

void SparseMatrix::print(){
	set<Edge>::iterator it;
	for(int i=0; i<size; i++){
		for(it = getIterFRowBegin(i); it != getIterFRowEnd(i); it++){
			cout<<it->src<< " " <<it->dst << " "<< it->value << endl;
		}
	}
}

void SparseMatrix::saveToLocal(string name, string address, map<int, int>* realId){
	ofstream file;
    string filename= address + name;
    cout<<"begin to save "<<filename<<endl;
    file.open(filename, ios::out);
	file << "nodeNum: " << size << "\n";
    file<<"edgeNum: "<< nonZeroNum <<"\n";

    if(realId == NULL){
	    set<Edge>::iterator it, it1, it2;
	    for(int i=0; i<size; i++){
	    	it1 = rowCol->at(i)->begin();
	    	it2 = rowCol->at(i)->end();
	    	for(it = it1; it!= it2; it++){
	    		file << it->src << " " << it->dst << " " << it->value << endl;
	    	}
	    }
	}
	else{
		set<Edge>::iterator it, it1, it2;
	    for(int i=0; i<size; i++){
	    	it1 = rowCol->at(i)->begin();
	    	it2 = rowCol->at(i)->end();
	    	for(it = it1; it!= it2; it++){
	    		file << realId->at(it->src) << " " << realId->at(it->dst) << " " << it->value << endl;
	    	}
	    }
	}
    cout<<"finish to save "<<filename<<endl;
    file.close();
}


void SparseMatrix::sparsifyByRatio(double beta){
	set<Edge>::iterator it;
	for(int i=0; i<this->getSize(); i++){
		set<EdgeValue>* edgeValue = new set<EdgeValue>;

		double tempSum = 0;
		for(it = getIterFRowBegin(i); it != getIterFRowEnd(i); it++){		
			tempSum = tempSum + it->value;
			edgeValue->insert(EdgeValue(it->src, it->dst, it->value));
		}
		tempSum = tempSum * beta;
		double sum = 0;
		set<EdgeValue>::iterator itValue;
		for(itValue = edgeValue->begin(); itValue != edgeValue->end(); itValue++){
			sum = sum + itValue->value;
			if(sum > tempSum) break;
			this->erase(Edge(itValue->src, itValue->dst, itValue->value));
		}
		delete edgeValue;
	}
}

void SparseMatrix::sparsifyByValue(double minValue){
	set<Edge>::iterator it;
	for(int i=0; i<this->getSize(); i++){
		for(it = getIterFRowBegin(i); it != getIterFRowEnd(i); ){
			if(it->value < minValue){
				Edge temp(it->src, it->dst, it->value);
				it++;
				this->erase(temp);
			}
			else{
				it++;
			}
		}
	}
}

void SparseMatrix::rowNormalize(){
	set<Edge>::iterator it;
	for(int i=0; i<this->getSize(); i++){
		double sum = 0;
		for(it = getIterFRowBegin(i); it != getIterFRowEnd(i); it++){
			sum = sum + it->value;
		}
		for(it = getIterFRowBegin(i); it != getIterFRowEnd(i); ){
			Edge temp1(it->src, it->dst, it->value);
			Edge temp2(it->src, it->dst, it->value/sum);
			it++;
			this->erase(temp1);
			this->insert(temp2);
		}
		
	}
}

bool SparseMatrix::ifExist(Edge edge){
	return (rowCol->at(edge.src)->find(edge))!=(rowCol->at(edge.src)->end());
}

int SparseMatrix::outDegree(int node_id){
	return rowCol->at(node_id)->size();
}

int SparseMatrix::inDegree(int node_id){
	return colRow->at(node_id)->size();
}