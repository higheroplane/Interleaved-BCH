#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime> 
#include <cstdlib>
#include <assert.h>
#include <vector>
#include <boost/math/special_functions/binomial.hpp>

void make_hilbert_errors (std::vector<int*>& errv, int pgb, int pbg, int pe, int N);

std::vector<int> numerical_encoder (int k, int len, int weight);


