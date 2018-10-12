// The file test_utils.tpp contains the templated definitions of the 
// functions with prototypes in the file test_utils.h.

template <class ComplexType, class RealType>
ComplexType random_complex ()
{
   double angle = 3.14*((double) rand())/RAND_MAX;
   double re_part = cos(angle);
   double im_part = sin(angle);

   RealType real_part = RealType(re_part);
   RealType imag_part = RealType(im_part);

   ComplexType result = ComplexType(real_part, imag_part);

   return result;
}

template <class ComplexType, class RealType>
void random_point ( int dim, ComplexType* coordinates )
{
   for(int idx=0; idx<dim; idx++)
      coordinates[idx] = random_complex<ComplexType,RealType>();
}

void random_exponents ( int dim, int expmax, int* exponents )
{
   for(int idx=0; idx<dim; idx++)
      exponents[idx] = rand() % (expmax+1);
}

template <class ComplexType, class RealType>
PolyMon<ComplexType,RealType> random_monomial ( int dim, int expmax )
{
   ComplexType ran = random_complex<ComplexType,RealType>();
   RealType coefficient[2];
   coefficient[0] = ran.real;
   coefficient[1] = ran.imag;

   int exponents[dim];
   random_exponents(dim,expmax,exponents);

   PolyMon<ComplexType,RealType> result
      = PolyMon<ComplexType,RealType>(dim, exponents, coefficient);

   return result;
}

template <class ComplexType, class RealType>
PolyMon<ComplexType,RealType> prompted_monomial ( int dim )
{
   RealType coefficient[2];

   cout << "-> give the real part of the coefficient : ";
   cin >> coefficient[0];
   cout << "-> give the imaginary part of the coefficient : ";
   cin >> coefficient[1];

   int exponents[dim];

   for(int idx=0; idx<dim; idx++)
   {
      cout << "-> give the exponent for variable " << idx+1 << " : ";
      cin >> exponents[idx];
   }

   PolyMon<ComplexType,RealType> result
      = PolyMon<ComplexType,RealType>(dim, exponents, coefficient);

   return result;
}

template <class ComplexType, class RealType>
ComplexType plain_eval ( PolyMon<ComplexType,RealType>& m, ComplexType* x )
{
   ComplexType result = m.coef;

   for(int varidx=0; varidx<m.n_var; varidx++)
   {
      int idx = m.pos[varidx];

      for(int expidx=0; expidx<m.exp[varidx]; expidx++)
         result = result*x[idx];
   }
   return result;
}

template <class ComplexType, class RealType>
void plain_diff
 ( PolyMon<ComplexType,RealType>& m, ComplexType *x, ComplexType *deri )
{
   for(int varidx=0; varidx<m.n_var; varidx++)
   {
      int idx = m.pos[varidx];

      deri[varidx] = m.coef;
      if(m.exp[varidx] > 1)
      {
         ComplexType factor = ComplexType(m.exp[varidx],0.0);
         deri[varidx] = factor*deri[varidx];
         for(int expidx=0; expidx<m.exp[varidx]-1; expidx++)
            deri[varidx] = deri[varidx]*x[idx];
      }
      for(int otheridx=0; otheridx<m.n_var; otheridx++)
         if(otheridx != varidx)
         {
            idx = m.pos[otheridx];
            for(int expidx=0; expidx<m.exp[otheridx]; expidx++)
               deri[varidx] = deri[varidx]*x[idx];
         }
   }
}

template <class ComplexType, class RealType>
void print_data ( PolyMon<ComplexType,RealType>& monomial )
{
   cout << endl << "The coefficient : " << monomial.coef << endl;

   cout << "dim : " << monomial.dim;
   cout << "  n_var : " << monomial.n_var;
   cout << "  n_base : " << monomial.n_base << endl;

   cout << endl << "The position array of the monomial :";
   for(int idx=0; idx<monomial.n_var; idx++)
      cout << " " << monomial.pos[idx];
   cout << endl;

   cout << "The exponent array of the monomial :";
   for(int idx=0; idx<monomial.n_var; idx++)
      cout << " " << monomial.exp[idx];
   cout << endl;

   cout << endl << "The base position array of the monomial :";
   for(int idx=0; idx<monomial.n_base; idx++)
      cout << " " << monomial.pos_base[idx];
   cout << endl;

   cout << "The base exponent array of the monomial :";
   for(int idx=0; idx<monomial.n_base; idx++)
      cout << " " << monomial.exp_base[idx];
   cout << endl << endl;
}

template <class ComplexType, class RealType>
PolyEq<ComplexType,RealType> random_polynomial
 ( int dim, int nbterms, int expmax )
{
   PolyEq<ComplexType,RealType> polynomial(dim);
   polynomial.constant = random_complex<ComplexType,RealType>();

   polynomial.n_mon = nbterms-1;

   for(int idx=0; idx<nbterms-1; idx++)
   {
      PolyMon<ComplexType,RealType>* term
         = new PolyMon<ComplexType,RealType>
                  (random_monomial<ComplexType,RealType>(dim,expmax));
      polynomial.mon.push_back(term);
   }
   return polynomial;
}

template <class ComplexType, class RealType>
void write_monomial ( int dim, PolyMon<ComplexType,RealType>& m )
{
   m.print_tableau(dim);
}

template <class ComplexType, class RealType>
void write_polynomial ( PolyEq<ComplexType,RealType>& p )
{
   cout << scientific << setprecision(4);

   cout << p.constant.real << "  " << p.constant.imag;
   for(int idx=0; idx<p.dim; idx++) cout << " 0";
   cout << endl;

   for(int idx=0; idx<p.n_mon; idx++)
   {
      write_monomial<ComplexType,RealType>(p.dim, *p.mon[idx]);

      cout << endl;
   }
}

template <class ComplexType, class RealType>
ComplexType plain_eval
 ( PolyEq<ComplexType,RealType>& p, ComplexType* x )
{
   ComplexType result = p.constant;

   for(int idx=0; idx<p.n_mon; idx++)
      result += plain_eval<ComplexType,RealType>(*p.mon[idx],x);

   return result;
}

template <class ComplexType, class RealType>
void plain_diff
 ( PolyEq<ComplexType,RealType>& p, ComplexType* x, ComplexType *deri )
{
   ComplexType monderi[p.dim];

   for(int idx=0; idx<p.dim; idx++) deri[idx].init(0.0,0.0);

   for(int idx=0; idx<p.n_mon; idx++)
   {
      PolyMon<ComplexType,RealType> *m = p.mon[idx];

      plain_diff<ComplexType,RealType>(*m, x, monderi);

      for(int jdx=0; jdx<m->n_var; jdx++)
         deri[m->pos[jdx]] += monderi[jdx];
   }
}

template <class ComplexType, class RealType>
void powertable
 ( int dim, int* maxdeg, ComplexType* point, ComplexType** powers )
{
   for(int idx=0; idx<dim; idx++)
   {
      int size = maxdeg[idx]+1;
      powers[idx][0] = point[idx];
      for(int powidx=1; powidx<size; powidx++)
         powers[idx][powidx] = powers[idx][powidx-1]*point[idx];
   }
}
