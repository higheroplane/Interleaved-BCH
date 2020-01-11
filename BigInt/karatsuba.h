#include <stdio.h>
#include <vector>

using namespace std;

void print_vector (vector <int> vec);
void print_ivector (vector <int> vec);
void finalize(vector <int> num);

vector<int> naive_mul(const vector<int>& x, const vector<int>& y);
vector<int> karatsuba_mul(const vector<int>& x, const vector<int>& y);
