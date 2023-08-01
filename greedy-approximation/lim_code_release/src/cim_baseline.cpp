#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <algorithm>
#include <cassert>
#include <complex> 
#include "cim_baseline.h"


using namespace std;

double UD::function(double x)
{
	return x;
	//return 2 * x - x * x;
}

void UD::_generate_RR_set(SparseMatrix* gf, int start_node, vector<bool>* active, vector<int>* RR_set){
	active->at(start_node) = true;
	RR_set->push_back(start_node);
	for(set<Edge>::iterator iter = gf->getIterFColBegin(start_node); iter != gf->getIterFColEnd(start_node); iter++){
		int target = iter->src;
		if(active->at(target)){
			continue;
		}
		double prob = iter->value;
		bool if_activated = ((rand() / (RAND_MAX+1.0)) < prob);
		if(if_activated){
			_generate_RR_set(gf, target, active, RR_set);
		}
	}
}

vector<vector<int>*>* UD::generate_RR_sets(SparseMatrix* gf, int theta){
	vector<vector<int>*>* res = new vector<vector<int>*>;

	int n = gf->getSize();

	for (int it = 0; it < theta; it++){
        vector<bool>* active = new vector<bool>;
        active->resize(n, false);
		vector<int>* RR_set = new vector<int>;
		int start_node = rand()%n;
		_generate_RR_set(gf, start_node, active, RR_set);
		delete active;
		res->push_back(RR_set);
	}
	return res;
}

pair<double, vector<int>*>* UD::greedy(vector<vector<int>*>* RR_sets, int budget_size, double c, int node_num, int n)
{
	pair<double, vector<int>*>* res = new pair<double, vector<int>*>;
	res->first = 0;
	res->second = new vector<int>;

	vector<bool>* if_selected = new vector<bool>;
	if_selected->resize(node_num, false);

	int theta = RR_sets->size();
	vector<double>* s = new vector<double>;
	s->resize(theta, 1.0);
	vector<double>* delta = new vector<double>;
	delta->resize(n, 0.0);

	for(int i = 0; i < budget_size; i++){
		for(vector<double>::iterator it = delta->begin();  it < delta->end(); it++)
		{
			*it = 0.0;
		}

		for(int j = 0; j < theta; j++){
			vector<int>* RR_set = RR_sets->at(j);
			for(vector<int>::iterator it_node = RR_set->begin(); it_node != RR_set->end(); it_node++)
			{
				int id = *it_node;
				if(!if_selected->at(id)){
					delta->at(id) += s->at(j) * (function(c) - function(0)) / (1 - function(0));
				}
			}
		}

		int id = 0;
		int opt_id;
		double max = 0;
		for(vector<double>::iterator it = delta->begin(); it != delta->end(); it++){
			if(*it > max){
				opt_id = id;
				max = *it;
			}
			id++;
		}

		if_selected->at(opt_id) = true;

		for(int j = 0; j < theta; j++){
			vector<int>* RR_set = RR_sets->at(j);
			for(vector<int>::iterator it_node = RR_set->begin(); it_node < RR_set->end(); it_node++)
			{
				id = *it_node;
				if(id == opt_id){
					s->at(j) = s->at(j) * (1 - function(c)) / (1 - function(0));
					break;
				}
			}
		}
	}

	double spread = 0;
	for(int i = 0; i < theta; i++){
		spread = spread + 1 - s->at(i);
	}

	spread *= (double (n) / theta);
	res->first = spread;

	for(int i = 0; i < node_num; i++){
		if(if_selected->at(i)){
			res->second->push_back(i);
		}
	}

	delete if_selected;
	delete delta;
	delete s;
	return res;
}

pair<double, vector<int>*>* UD::run(SparseMatrix* gf, int k, double c, int theta, vector<vector<int>*>* RR_sets){
	pair<double, vector<int>*>* res = new pair<double, vector<int>*>;
	res->first = 0; //optimal c under exhaustive search
	res->second = new vector<int>; //selected ndoes

	if(RR_sets == NULL){
		RR_sets = generate_RR_sets(gf, theta);
	}

	int n = gf->getSize();

	int max_round = int(1.0 / c);

	double opt_c = 0;

	for(int i = 0; i < max_round; i++){
		int budget_size = int(k / (c * (i + 1)));
		pair<double, vector<int>*>* temp = greedy(RR_sets, budget_size, c * (i + 1), n, n);
		if(temp->first > res->first){
			delete res->second;
			delete res;
			res = temp;
			opt_c = c * (i + 1);
		}
	}

	res->first = opt_c;

	return res;
}

double CD::function(double x)
{
	return x;
	//return 2 * x - x * x;
}

complex<double> neg_cubic_pow(complex<double> x){
	if(real(x) < 0){
		return -pow(-x, 1.0 / 3);
	}
	else{
		return pow(x, 1.0 / 3);
	}
}

vector<double>* cubic_root(double a, double b, double c, double d){
	vector<double>* res = new vector<double>;
	double eps = 1e-8;
	double f = ((3 * c / a) - (b * b / (a * a))) / 3;
	double g = ((2 * b * b * b / (a * a * a)) - (9 * b * c / (a * a)) + (27 * d / a)) / 27;
	double h = g * g / 4 + f * f * f / 27;

	if (abs (f) <= eps && abs (g) <= eps && abs (h) <= eps) { // all 3 roots are real and equal
		res->push_back(-1.0 * pow (d / a, 1.0 / 3.0));
	} else if (h > 0) { //only 1 root is real
		double R = -(g / 2) + pow (h, 0.5);
		double S = pow (R, 1.0 / 3.0);
		double T = -(g / 2) - pow (h, 0.5);
		double U;
		if (T < 0)
			U = -1 * pow (-T, 1.0 / 3.0);
		else
			U = pow (T, 1.0 / 3.0);
//				Console.WriteLine ("{0}, {1}, {2}, {3}", R, S, T, U);
		res->push_back(S + U - b / (3 * a));
	} else if (h <= 0) { // all 3 roots are real
		double i = pow(g*g/4-h, 0.5);
		double j = pow (i, 1.0 / 3.0);
		double K = acos (-g / (2 * i));
		double L = -j;
		double M = cos (K / 3);
		double N = pow (3.0, 0.5) * sin (K / 3);
		double P = -1 * (b / (3 * a));
//				Console.WriteLine ("{0}, {1}, {2}, {3}, {4}, {5}, {6}", i, j, K, L, M, N, P);
		res->push_back(2 * j * cos (K / 3) - b / (3 * a));
		res->push_back(L * (M + N) + P);
		res->push_back(L * (M - N) + P);
	}

	return res;
}


vector<double>* quadratic_root(double a, double b, double c){
	vector<double>* res = new vector<double>;
	if(b * b - 4 * a * c > 0){
		res->push_back((-b + pow(b * b - 4 * a * c, 0.5)) / (2 * a));
		res->push_back((-b - pow(b * b - 4 * a * c, 0.5)) / (2 * a));
	}
	else if(b * b - 4 * a * c == 0){
		res->push_back(-b / (2 * a));
	}
	return res;
}

vector<double>* linear_root(double a, double b){
	vector<double>* res = new vector<double>;
	res->push_back(-b / a);
	return res;
}


void CD::_exchange(int id1, int id2, vector<double>* value, vector<double>* s, vector<vector<int>*>* RR_index){
	double sum = value->at(id1) + value->at(id2);

	vector<double>* x1;
	vector<double>* x2 = new vector<double>;
	vector<double>* new_value = new vector<double>;

	//case 1:
	double A = 0;
	double B = 0;
	double C = 0;

	set<int>* S_id1 = new set<int>;
	set<int>* S_id2 = new set<int>;
	set<int>* S_id1_id2 = new set<int>;

	for(int i = 0; i < RR_index->at(id1)->size(); i++){
		S_id1->insert(RR_index->at(id1)->at(i));
	}
	for(int i = 0; i < RR_index->at(id2)->size(); i++){
		if(S_id1->find(RR_index->at(id2)->at(i)) == S_id1->end()){
			S_id2->insert(RR_index->at(id2)->at(i));
		}
		else{
			S_id1_id2->insert(RR_index->at(id2)->at(i));
			S_id1->erase(RR_index->at(id2)->at(i));
		}

	}

	set<int>::iterator it;
	for(it = S_id1->begin(); it != S_id1->end(); it++){
		A = A + s->at(*it) / (1 - function(value->at(id1)));
	}
	for(it = S_id2->begin(); it != S_id2->end(); it++){
		C = C + s->at(*it) / (1 - function(value->at(id2)));
	}
	for(it = S_id1_id2->begin(); it != S_id1_id2->end(); it++){
		B = B + s->at(*it) / ((1 - function(value->at(id1))) * (1 - function(value->at(id2))));
	}

	double a = 2 * B;
	double b = -3 * sum * B;
	double c = A + C + (sum * sum + 2 * sum - 2) * B;
	double d = -(sum * sum - sum) * B - A - (sum - 1) * C;

	if(a != 0){
		x1 = cubic_root(a, b, c, d);
	}
	else if(b != 0){
		x1 = quadratic_root(b, c, d);
	}
	else{
		x1 = linear_root(c, d);
	}
	
	for(int i = 0; i < x1->size(); i++){
		if(x1->at(i) < max(0.0, sum - 1) || x1->at(i) > min(1.0, sum)){
			x1->at(i) = max(0.0, sum - 1);
		}
		x2->push_back(sum - x1->at(i));
	}

	//case 2:
	x1->push_back(max(0.0, sum - 1));
	x2->push_back(sum - x1->at(x1->size() - 1));

	//case 3:
	x1->push_back(min(1.0, sum));
	x2->push_back(sum - x1->at(x1->size() - 1));

	//case 4:
	x1->push_back(value->at(id1));
	x2->push_back(value->at(id2));
	
	new_value->resize(x1->size(), 0.0);
	for(int j = 0; j < new_value->size(); j++){
		new_value->at(j) += A * (1 - function(x1->at(j)));
		new_value->at(j) += B * (1 - function(x1->at(j))) * (1 - function(x2->at(j)));
		new_value->at(j) += C * (1 - function(x2->at(j)));
	}

	

	double min_value = new_value->at(0) + 1;
	int min_id;
	for(int i = 0; i < new_value->size(); i++){
		if(new_value->at(i) < min_value){
			min_value = new_value->at(i);
			min_id = i;
		}
	}

	for(it = S_id1->begin(); it != S_id1->end(); it++){
		s->at(*it) /= (1 - function(value->at(id1)));
		s->at(*it) *= (1 - function(x1->at(min_id)));
	}
	for(it = S_id2->begin(); it != S_id2->end(); it++){
		s->at(*it) /= (1 - function(value->at(id2)));
		s->at(*it) *= (1 - function(x2->at(min_id)));
	}
	for(it = S_id1_id2->begin(); it != S_id1_id2->end(); it++){
		s->at(*it) /= (1 - function(value->at(id1)));
		s->at(*it) *= (1 - function(x1->at(min_id)));
		s->at(*it) /= (1 - function(value->at(id2)));
		s->at(*it) *= (1 - function(x2->at(min_id)));
	}

	value->at(id1) = x1->at(min_id);
	value->at(id2) = x2->at(min_id);

	delete x1;
	delete x2;
	delete S_id1;
	delete S_id2;
	delete S_id1_id2;
}

vector<double>* CD::run(SparseMatrix* gf, int k, double c, int theta){
	int n = gf->getSize();

	UD* ud = new UD();
	vector< vector<int>*>* RR_sets = ud->generate_RR_sets(gf, theta);
	pair<double, vector<int>*>* res_ud = ud->run(gf, k, c, theta, RR_sets);

	vector<double>* value = new vector<double>;
	value->resize(n, 0);

	for(int i = 0; i < res_ud->second->size(); i++){
		value->at(res_ud->second->at(i)) = res_ud->first;
	}

	//preparation
	vector<double>* s = new vector<double>;
	s->resize(RR_sets->size(), 1.0);

	vector<vector<int>*>* RR_index = new vector<vector<int>*>;
	for(int i = 0; i < n; i++){
		RR_index->push_back(new vector<int>);
	}

	for(int i = 0; i < RR_sets->size(); i++){
		for(int j = 0; j < RR_sets->at(i)->size(); j++){
			s->at(i) *= (1 - function(value->at(RR_sets->at(i)->at(j))));
			RR_index->at(RR_sets->at(i)->at(j))->push_back(i);
		}
	}
	
	for(int i = 0; i < 100; i++){
		for(int j = 0; j < res_ud->second->size(); j++){
			for(int l = 0; l < res_ud->second->size(); l++){
				if(j == l) continue;
				if(abs(1 - value->at(res_ud->second->at(j))) < 1e-8 || abs(1 - value->at(res_ud->second->at(l))) < 1e-8) continue;
				_exchange(res_ud->second->at(j), res_ud->second->at(l), value, s, RR_index);
			}
		}
	}

	return value;
}

vector<double>* HeuristicsDegree::run(SparseMatrix* gf, int k, int M){
	vector<double>* res = new vector<double>;
	res->resize(gf->getSize(), 0.0);
	set<pair<int, int>>* down_degree = new set<pair<int, int>>;
	for(int i = 0; i < gf->getSize(); i++){
		down_degree->insert(pair<int, int>(gf->inDegree(i) + gf->outDegree(i), i));
	}

	set<pair<int, int>>::reverse_iterator it;
	int m = 0;
	double sum = 0;
	for(it = down_degree->rbegin(); it != down_degree->rend(); it++){
		sum += it->first;
		m++;
		if(m > M) break;
	}

	m = 0;
	for(it = down_degree->rbegin(); it != down_degree->rend(); it++){
		res->at(it->second) = double(k) * it->first / sum;
		m++;
		if(m > M) break;
	}

	return res;
}

double OriginalGreedy::_spread(int start_node, vector<bool>* active){
	double res = 1;
	active->at(start_node) = true;
	for(set<Edge>::iterator iter = gf->getIterFRowBegin(start_node); iter != gf->getIterFRowEnd(start_node); iter++){
		int target = iter->dst;
		if(active->at(target)) continue;
		double prob = iter->value;
		bool if_activated = ((rand() / (RAND_MAX+1.0)) < prob);
		if(if_activated){
			res += _spread(target, active);
		}
	}
	return res;
}

double OriginalGreedy::spread(vector<double>* prob){
	double res = 0;
	for(int i = 0; i < sample_num; i++){
		vector<int>* seed = new vector<int>;
		for(int i = 0; i < prob->size(); i++){
			if((rand() / (RAND_MAX+1.0)) < prob->at(i)){
				seed->push_back(i);
			}
		}
		vector<bool>* active = new vector<bool>;
		active->resize(gf->getSize(), false);
		for(int i = 0; i < seed->size(); i++){
			if(active->at(seed->at(i))) continue;
			res += _spread(seed->at(i), active);
		}
		delete active;
		delete seed;
	}
	return res / sample_num;
}

vector<double>* OriginalGreedy::prob_1(vector<double>* value){
	vector<double>* prob = new vector<double>;
	prob->resize(value->size(), 0.0);
	for(int i = 0; i < prob->size(); i++){
		prob->at(i) = 2 * value->at(i) - value->at(i) *  value->at(i);
	}
	return prob;
}

vector<double>* OriginalGreedy::prob_2(vector<double>* value){
	vector<double>* prob = new vector<double>;
	prob->resize(gf->getSize(), 1.0);
	for(int i = 0; i < prob_activate->size(); i++){
		vector<pair<int, double>*>* tmp_vec_pair = prob_activate->at(i);
		for(int j = 0; j < tmp_vec_pair->size(); j++){
			prob->at(tmp_vec_pair->at(j)->first) *= pow(tmp_vec_pair->at(j)->second, value->at(i));
		}
	}
	for(int i = 0; i < gf->getSize(); i++){
		prob->at(i) = 1 - prob->at(i);
	}
	return prob;
}

double OriginalGreedy::run(int k, double delta, int cur_k, double spread_point, vector<double>* res, vector<int>* iter_number, priority_queue<pair<double, int>, vector<pair<double, int>>, less<pair<double, int>>>* q){
	int d;
	if(type == 0){
		d = gf->getSize();
	}
	else if(type == 1){
		d = prob_activate->size();
	}
	if(res->size() == 0){
		res->resize(d, 0.0);
	}
	if(iter_number->size() == 0){
		iter_number->resize(d, 0);
	}

	vector<double>* prob = NULL;

	if(q->size() == 0){
		for(int i = 0; i < d; i++){
			res->at(i) = res->at(i) + delta;
			if(type == 0){
				prob = prob_1(res);
			}
			else if(type == 1){
				prob = prob_2(res);
			}
			double temp_spread = spread(prob);
			delete prob;
			res->at(i) = res->at(i) - delta;
			pair<double, int> temp_pair(temp_spread - spread_point, i);
			q->push(temp_pair);
		}
	}

	for(int i = cur_k; i < cur_k + k / delta; i++){
		cout << "iter_number = " << i << endl;
		bool loop = true;
		while(loop){
			pair<double, int> temp_pair = q->top();
			int id = temp_pair.second;
			if(iter_number->at(id) < i){
				res->at(id) += delta;
				if(type == 0){
					prob = prob_1(res);
				}
				else if(type == 1){
					prob = prob_2(res);
				}
				double temp_spread = spread(prob);
				res->at(id) -= delta;
				delete prob;
				q->pop();
				double spread_delta = temp_spread - spread_point;
				pair<double, int> insert_pair(spread_delta, id);
				q->push(insert_pair);
				iter_number->at(id) = i;
			}
			else{
				spread_point += temp_pair.first;
				res->at(id) += delta;

				res->at(id) += delta;
				if(type == 0){
					prob = prob_1(res);
				}
				else if(type == 1){
					prob = prob_2(res);
				}
				double temp_spread = spread(prob);
				res->at(id) -= delta;
				delete prob;
				q->pop();
				double spread_delta = temp_spread - spread_point;
				pair<double, int> insert_pair(spread_delta, id);
				q->push(insert_pair);
				iter_number->at(id) = i + 1;
				loop = false;
			}
		}
	}
	return spread_point;
}