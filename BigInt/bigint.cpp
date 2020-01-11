#include "bigint.h"

const BigInt plus(const BigInt& _left, const BigInt& _right) 
{
    BigInt left (_left), right (_right);        
    if (right > left) 
    {        
        right.num.push_back (0);
        for (int i = 0; i < left.num.size(); i++) right.num[i] += left.num[i];
        right.finalize();
        return right;
    }
    else
    {    
        left.num.push_back (0);    
        for (int i = 0; i < right.num.size(); i++) left.num[i] += right.num[i];
        left.finalize();
        return left;
    }
}

const BigInt mminus(const BigInt& _left, const BigInt& _right) 
{
    BigInt left (_left), right (_right);        
    if (right > left) 
    {        
        right.num.push_back (0);
        for (int i = 0; i < left.num.size(); i++) right.num[i] -= left.num[i];
        right.finalize();
        return right;
    }
    else
    {    
        left.num.push_back (0);    
        for (int i = 0; i < right.num.size(); i++) left.num[i] -= right.num[i];
        left.finalize();
        return left;
    }
}

const BigInt operator+(const BigInt&  left, const BigInt&  right)
{
    if ( left.sign &&  right.sign) return -(plus(left, right));
    else if (!left.sign && !right.sign) return (plus(left, right));
    else if (abs(left) >= abs(right))
    {
        if ( left.sign && !right.sign) return -(mminus (left, right));
        if (!left.sign &&  right.sign) return  (mminus (left, right));
    }
    else
    {
        if ( left.sign && !right.sign) return  (mminus (right, left));
        if (!left.sign &&  right.sign) return -(mminus (right, left));
    }
}

const BigInt operator-(const BigInt&  left, const BigInt&  right)
{
    return left + (-right);
}

const BigInt operator-(const BigInt&  val)
{
    return BigInt (!val.sign, val.num);
}

const BigInt operator*(const BigInt& left, const BigInt& right) 
{
    BigInt ret (karatsuba_mul(left.num, right.num));
    ret.finalize();
    return ret;
}

const BigInt operator/(const BigInt& _left, const BigInt& _right) 
{
    //BigInt::resize_mode = false;
    BigInt q (_left.size()), r(_right.size()), left(_left), right(_right);

    int llen = _left.size() + 1;
    int rlen = _right.size() + 1;
    int debug = 0;

    for (int j = 0; j <= _left.size() - _right.size(); j ++)
    {
        //if (debug) printf ("rlen = %d, llen = %d\n", rlen, llen);
        left.num.insert(left.num.end(), llen - rlen, 0);
        BigInt temp (vector <int> {left.num.begin() + _left.size() - _right.size() - j, left.num.begin() + _left.size()});
        //if (debug) printf ("rlen = %d, llen = %d\n", rlen, llen);
        //temp.num.insert(temp.num.end(), llen - rlen, 0);
        left.finalize();
        
        int i = 0;
        for (i = 1; BigInt(vector<int>({i})) * right <= temp && i <= 11; i ++){}
        if (i == 11) i = 0;
        else i --;

        if (debug)
        {
            printf ("left :");
            left.print();
            printf ("right:");
            right.print();
            printf ("temp :");
            temp.print();
        }

        
        q.num.insert (q.num.begin(), 1, i);
        if (debug) printf ("q:");
        if (debug) q.print();
        
        if (debug) printf ("i = %d len= %d\n", i, left.size() - right.size());
        if (debug) printf ("r*ye :");

        BigInt eye (vector<int> (llen - rlen - j + 1, 0));
        eye.num[llen - rlen - j] = i;
        if (debug) (right*eye).print();
        //if (debug) print_ivector((right*eye).num);
        
        left = left - (right*eye);
        left.finalize();
        //(right*BigInt(vector<int>({i}))*BigInt(vector<int>({(int) pow(10, _left.size() - _right.size() - j)}))).print();
        if (debug) left.print();
        if (debug) print_ivector (left.num);
        if (debug) printf("\n\n");
        if (debug)getchar();
    }   
    //left.print();
  //  printf("\n\n");
    
    if (q.num.size() != 0) return q;
    else return BigInt (0, vector<int>({0}));
    //BigInt::resize_mode = true;
}

const BigInt operator%(const BigInt& _left, const BigInt& _right) 
{
    BigInt q(_left.size()), r(_right.size()), left(_left), right(_right);
    for (int j = 0; left >= right; j ++)
    {
        int llen = left.size() + 1;
        int rlen = right.size() + 1;
        //printf ("rlen = %d\n", rlen);
        BigInt temp (vector <int> {left.num.begin() + _left.size() - _right.size() - j, left.num.begin() + _left.size()});
        //temp.num.insert(temp.num.begin(), llen - rlen, 0);
        
        int i = 0;
        for (i = 1; BigInt(vector<int>({i})) * right <= temp && i <= 10; i ++){}
        if (i == 10) i = 0;
        else i --;
    //  printf ("left:");
    //  left.print();
    //   printf ("right:");
    //  right.print();
    // printf ("temp:");
    // temp.print();

        
        //q.num.insert (q.num.begin(), 1, i);
    //  printf ("q:");
    // q.print();
        
    // printf ("i = %d len= %d\n", i, llen - rlen);
    // printf ("ne left:");
        
        left = left - right*BigInt(vector<int>({i}))*BigInt(vector<int>({(int) pow(10, _left.size() - _right.size() - j)}));
    // left.print();
    // printf("\n\n");
    }
    //left.print();
  //  printf("\n\n");
    
    return left;

}

const bool operator>(const BigInt& left, const BigInt& right) 
{
    if      (!left.sign &&  right.sign) return 1;  
    else if ( left.sign && !right.sign) return 0;   
    if      (left.size() > right.size()) return 1 + left.sign;
    else if (right.size() > left.size()) return 0 + left.sign;
    else
    {
        for (int i = left.size() - 1; i >= 0; i--)
        {
            if (left [i] > right[i]) return 1 + left.sign;
            if (right[i] > left [i]) return 0 + left.sign;
        }
        return 0 + left.sign;
    }
}

const bool operator<(const BigInt& left, const BigInt& right) 
{
    if      (!left.sign &&  right.sign) return 0;  
    else if ( left.sign && !right.sign) return 1; 
    if      (left.size() < right.size()) return 1 + left.sign;
    else if (right.size() < left.size()) return 0 + left.sign;
    else
    {
        for (int i = left.size() - 1; i >= 0; i--)
        {
            if (left [i] < right[i]) return 1 + left.sign;
            if (right[i] < left [i]) return 0 + left.sign;
        }
        return 0 + left.sign;
    }
}

const bool operator>=(const BigInt& left, const BigInt& right) 
{
    return !(left < right);
}

const bool operator<=(const BigInt& left, const BigInt& right) 
{
    return !(left > right);
}

const bool operator==(const BigInt& left, const BigInt& right) 
{
    if      (left.sign != right.sign) return 0;
    else if (right.size() < left.size()) return 0;
    else if (right.size() > left.size()) return 0;
    else
    {
        for (int i = left.size() - 1; i >= 0; i--)
            if (left[i] != right[i]) return 0;
        return 1;
    }
}

const bool operator!=(const BigInt& left, const BigInt& right) 
{
    return !(left == right);
}


const BigInt abs (const BigInt& val)
{
    return BigInt(0, val.num);
}

void BigInt::print() const
{
    if (sign == 1) printf ("-");   
    int if_zeros = 0; 
    for (int i = num.size() - 1; i > 0; i --)
    {
        if (if_zeros == 0) if(num[i] != 0) if_zeros = 1;
        if (if_zeros == 1) printf ("%d", num[i]);
    }
    printf ("%d", num[0]);
    printf ("\n");
}

void BigInt::finalize()
{
    num.insert(num.end(), 7, 0);
    for (int i = 0; i < num.size(); i ++)
    {

        num[i + 1] += num[i] / 10;
        if (num[i] >= 0) num[i] %= 10;
        else {num[i] = 10 + num[i]%10; num[i+1] --;}
    }
    num.resize (size());
}
