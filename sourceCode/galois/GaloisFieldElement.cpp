#include "GaloisFieldElement.h"
#include <assert.h>

namespace galois
{

   GaloisField * GaloisFieldElement::last_gf = nullptr;
   GaloisFieldElement::GaloisFieldElement(GaloisField* _gf, GFSymbol v)
   {
      if (_gf != NULL)
      {
         gf         = _gf;
         last_gf = gf;
         poly_value = v & gf->size();
      }
      else
        poly_value = v;
   }

   GaloisFieldElement::GaloisFieldElement(const int num)
   {
      assert (num == 1 || num == 0);
      gf = last_gf; 
      if (num == 0) poly_value = 0;
      else if (num == 1) poly_value = 1;
      else (exit(1));
   }



   GaloisFieldElement::GaloisFieldElement(const GaloisFieldElement& gfe)
   {
      gf          = gfe.gf;
      last_gf = gf;
      poly_value  = gfe.poly_value;
   }


   std::ostream& operator << (std::ostream& os, const GaloisFieldElement& gfe)
   {
      os << gfe.poly_value;
      return os;
   }


   GaloisFieldElement operator+(const GaloisFieldElement& a, const GaloisFieldElement& b)
   {
      GaloisFieldElement result  = a;
      result += b;
      return result;
   }


   GaloisFieldElement operator-(const GaloisFieldElement& a, const GaloisFieldElement& b)
   {
      GaloisFieldElement result  = a;
      result -= b;
      return result;
   }


   GaloisFieldElement operator*(const GaloisFieldElement& a, const GaloisFieldElement& b)
   {
      GaloisFieldElement result  = a;
      result *= b;
      return result;
   }


   GaloisFieldElement operator*(const GaloisFieldElement& a, const GFSymbol& b)
   {
      GaloisFieldElement result  = a;
      result *= b;
      return result;
   }


   GaloisFieldElement operator*(const GFSymbol& a, const GaloisFieldElement& b)
   {
      GaloisFieldElement result  = b;
      result *= a;
      return result;
   }


   GaloisFieldElement operator/(const GaloisFieldElement& a, const GaloisFieldElement& b)
   {
      GaloisFieldElement result  = a;
      result /= b;
      return result;
   }


   GaloisFieldElement operator^(const GaloisFieldElement& a, const int& b)
   {
      GaloisFieldElement result  = a;
      result ^= b;
      return result;
   }

}
