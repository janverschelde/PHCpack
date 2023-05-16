// The file complexinc_job.h defines the class ComplexIncrementJob
// to store the data of one job to increment the first operand,
// of either the real or imaginary part of complex data.

#include "complexinc_job.h"

ComplexIncrementJob::ComplexIncrementJob
 ( int monomix, int kind, int ix, bool isreal )
{
   monidx = monomix; // monomial index
   inkind = kind;    // kind of increment
   incidx = ix;      // index of the increment
   rornot = isreal;  // real or not?
}

int ComplexIncrementJob::get_monomial_index ( void ) const
{
   return monidx;
}

int ComplexIncrementJob::get_increment_kind ( void ) const
{
   return inkind;
}

int ComplexIncrementJob::get_increment_index ( void ) const
{
   return incidx;
}

bool ComplexIncrementJob::get_isreal ( void ) const
{
   return rornot;
}

std::ostream& operator<<
 ( std::ostream& os, const ComplexIncrementJob& job )
{
   os << "monomial " << job.monidx << " : ";

   if(job.inkind ==  1) os << "f[" << job.incidx << "]";
   if(job.inkind ==  2) os << "b[" << job.incidx << "]";
   if(job.inkind ==  3) os << "c[" << job.incidx << "]";
   if(job.rornot)
      os << "re^a += ";
   else
      os << "im^a += ";

   if(job.inkind ==  1) os << "f[" << job.incidx << "]";
   if(job.inkind ==  2) os << "b[" << job.incidx << "]";
   if(job.inkind ==  3) os << "c[" << job.incidx << "]";
   if(job.rornot)
      os << "re^b";
   else
      os << "im^b";

   return os;
}
