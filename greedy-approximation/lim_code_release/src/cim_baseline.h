#include <set>
#include <vector>
#include "SparseMatrix.h"
#include <queue>

using namespace std;
class UD{
public:
	double function(double x);
	pair<double, vector<int>*>* run(SparseMatrix* gf, int k, double c, int theta, vector<vector<int>*>* RR_sets = NULL);
	pair<double, vector<int>*>* greedy(vector<vector<int>*>* RR_sets, int budget_size, double c, int node_num, int n);
	vector< vector<int>*>* generate_RR_sets(SparseMatrix* gf, int theta);
	void _generate_RR_set(SparseMatrix* gf, int start_node, vector<bool>* active, vector<int>* RR_set);
};

class CD{
public:
	double function(double x);
	//vector<double>* cubic_root(double a, double b, double c, double d);
	vector<double>* run(SparseMatrix* gf, int k, double c, int theta);
	void _exchange(int id1, int id2, vector<double>* value, vector<double>* s, vector<vector<int>*>* RR_index);
};

class HeuristicsDegree{
public:
	vector<double>* run(SparseMatrix* gf, int k, int M);
};

class OriginalGreedy{
public:
	OriginalGreedy(SparseMatrix* gf, vector<vector<pair<int, double>*>*>* prob_activate=NULL, int type=0, int sample_num=100000){
		this->sample_num = sample_num;
		this->prob_activate = prob_activate;
		this->type = type;
		this->gf = gf;
	};
	SparseMatrix* gf;
	vector<vector<pair<int, double>*>*>* prob_activate;
	int type;
	double sample_num;
	double run(int k, double delta, int cur_k, double spread_point, vector<double>* res, vector<int>* iter_number, priority_queue<pair<double, int>, vector<pair<double, int>>, less<pair<double, int>>>* q);
	double spread(vector<double>* prob);
	vector<double>* prob_1(vector<double>* res);
	vector<double>* prob_2(vector<double>* res);
	double _spread(int start_node, vector<bool>* active);
};