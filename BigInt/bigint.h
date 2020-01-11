#include "karatsuba.h"
#include <cmath>

#define _BI(i) BigInt(vector<int>({i}))

class BigInt
{
    public:
    vector <int> num;
    bool sign = 0;
    static bool resize_mode;
    //int base;
    vector<int> naive_mul(const vector<int>& x, const vector<int>& y);
    vector<int> karatsuba_mul(const vector<int>& x, const vector<int>& y);
    void finalize();
    int len() const {return num.size();}
    int size() const
    {
        int size = 0; 
        for (size = len() - 1; num[size] == 0 && size > 0; size --);
        return size + 1;
    }
    void print()const;

    
    //BigInt(char* str);
    BigInt (){}
    BigInt (int cap){num.reserve(cap);}
    BigInt (vector <int> vec): num(vec) {finalize();}
    BigInt (bool _sign, vector <int> vec): num(vec), sign(_sign) {finalize();}
    BigInt (const BigInt& big): num (big.num), sign(big.sign) {finalize();}
    int operator[] (int i)const  {return num[i];}
    /*const BigInt& operator = (const BigInt& val)
    {
	    return BigInt (val.sign, val.num);
    }*/
};

const BigInt plus (const BigInt& _left, const BigInt& _right); 
const BigInt mminus(const BigInt& _left, const BigInt& _right);

const BigInt operator-(const BigInt&  val); 

const BigInt operator-(const BigInt& left, const BigInt& right);
const BigInt operator+(const BigInt& left, const BigInt& right);
const BigInt operator*(const BigInt& left, const BigInt& right);
const BigInt operator/(const BigInt& left, const BigInt& right);
const BigInt operator%(const BigInt& left, const BigInt& right);
const BigInt abs      (const BigInt& val);

const bool operator> (const BigInt& left, const BigInt& right);
const bool operator< (const BigInt& left, const BigInt& right);
const bool operator>=(const BigInt& left, const BigInt& right);
const bool operator<=(const BigInt& left, const BigInt& right);
const bool operator==(const BigInt& left, const BigInt& right);
const bool operator!=(const BigInt& left, const BigInt& right);


