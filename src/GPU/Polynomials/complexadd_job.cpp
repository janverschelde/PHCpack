// The file complexadd_job.cpp defines the methods of the class
// ComplexAdditionJob, specified in "complexadd_job.h".

#include "complexadd_job.h"

ComplexAdditionJob::ComplexAdditionJob
 ( int atp, int itp, int monix1, int monix2, int ix1, int ix2 )
{
   adtype = atp;
   intype = itp;
   updmon = monix1;
   incmon = monix2;
   updidx = ix1;
   incidx = ix2;
}

int ComplexAdditionJob::get_addition_type ( void ) const
{
   return adtype;
}

int ComplexAdditionJob::get_increment_type ( void ) const
{
   return intype;
}

int ComplexAdditionJob::get_update_monomial ( void ) const
{
   return updmon;
}

int ComplexAdditionJob::get_update_index ( void ) const
{
   return updidx;
}

int ComplexAdditionJob::get_increment_monomial ( void ) const
{
   return incmon;
}

int ComplexAdditionJob::get_increment_index ( void ) const
{
   return incidx;
}

std::ostream& operator<< ( std::ostream& os, const ComplexAdditionJob& job )
{
   if(job.adtype == 1)
   {
      // os << "f[" << job.updmon << "," << job.updidx << "]re += ";
      if(job.incmon < 0) // no distinction between -1 and -2
      {
         if(job.incidx < 0)
            os << "f[" << job.updmon << "," << job.updidx << "]re += "
               << "cstre";
         else
            os << "f[" << job.updmon << "," << job.updidx << "]re += "
               << "cff[" << job.incidx << "]re";
      }
      else
      {
         if(job.intype == 1)
            os << "f[" << job.updmon << "," << job.updidx << "]re += "
               << "f[" << job.incmon << "," << job.incidx << "]re^a";
         else if(job.intype == 2)
            os << "f[" << job.updmon << "," << job.updidx << "]re += "
               << "b[" << job.incmon << "," << job.incidx << "]re^a";
         else if(job.intype == 3)
            os << "f[" << job.updmon << "," << job.updidx << "]re += "
               << "c[" << job.incmon << "," << job.incidx << "]re^a";
         else if(job.intype == 4)
            os << "f[" << job.updmon << "," << job.updidx << "]re -= "
               << "f[" << job.incmon << "," << job.incidx << "]re^b";
         else if(job.intype == 5)
            os << "f[" << job.updmon << "," << job.updidx << "]re -= "
               << "b[" << job.incmon << "," << job.incidx << "]re^b";
         else if(job.intype == 6)
            os << "f[" << job.updmon << "," << job.updidx << "]re -= "
               << "c[" << job.incmon << "," << job.incidx << "]re^b";
         else
            os << "invalid operation";
      }
   }
   else if(job.adtype == 2)
   {
      // os << "b[" << job.updmon << "," << job.updidx << "]re += ";
      if(job.incmon < 0) // no distinction between -1 and -2
         os << "b[" << job.updmon << "," << job.updidx << "]re += "
            << "cff[" << job.incidx << "]re";
      else
      {
         if(job.intype == 1)
            os << "b[" << job.updmon << "," << job.updidx << "]re += "
               << "f[" << job.incmon << "," << job.incidx << "]re^a";
         else if(job.intype == 2)
            os << "b[" << job.updmon << "," << job.updidx << "]re += "
               << "b[" << job.incmon << "," << job.incidx << "]re^a";
         else if(job.intype == 3)
            os << "b[" << job.updmon << "," << job.updidx << "]re += "
               << "c[" << job.incmon << "," << job.incidx << "]re^a";
         else if(job.intype == 4)
            os << "b[" << job.updmon << "," << job.updidx << "]re -= "
               << "f[" << job.incmon << "," << job.incidx << "]re^b";
         else if(job.intype == 5)
            os << "b[" << job.updmon << "," << job.updidx << "]re -= "
               << "b[" << job.incmon << "," << job.incidx << "]re^b";
         else if(job.intype == 6)
            os << "b[" << job.updmon << "," << job.updidx << "]re -= "
               << "c[" << job.incmon << "," << job.incidx << "]re^b";
         else
            os << "invalid operation";
      }
   }
   else if(job.adtype == 3)
   {
      // os << "c[" << job.updmon << "," << job.updidx << "]re += ";
      if(job.incmon < 0) // no distinction between -1 and -2
         os << "c[" << job.updmon << "," << job.updidx << "]re += "
            << "cff[" << job.incidx << "]re";
      else
      {
         if(job.intype == 1)
            os << "c[" << job.updmon << "," << job.updidx << "]re += "
               << "f[" << job.incmon << "," << job.incidx << "]re^a";
         else if(job.intype == 2)
            os << "c[" << job.updmon << "," << job.updidx << "]re += "
               << "b[" << job.incmon << "," << job.incidx << "]re^a";
         else if(job.intype == 3)
            os << "c[" << job.updmon << "," << job.updidx << "]re += "
               << "c[" << job.incmon << "," << job.incidx << "]re^a";
         else if(job.intype == 4)
            os << "c[" << job.updmon << "," << job.updidx << "]re -= "
               << "f[" << job.incmon << "," << job.incidx << "]re^b";
         else if(job.intype == 5)
            os << "c[" << job.updmon << "," << job.updidx << "]re -= "
               << "b[" << job.incmon << "," << job.incidx << "]re^b";
         else if(job.intype == 6)
            os << "c[" << job.updmon << "," << job.updidx << "]re -= "
               << "c[" << job.incmon << "," << job.incidx << "]re^b";
         else
            os << "invalid operation";
      }
   }
   else if(job.adtype == 4)
   {
      os << "f[" << job.updmon << "," << job.updidx << "]im += ";
      if(job.incmon < 0) // no distinction between -1 and -2 ...
      {
         if(job.incidx < 0)
            os << "cstim";
         else
            os << "cff[" << job.incidx << "]im";
      }
      else
      {
         if(job.intype == 7)
            os << "f[" << job.incmon << "," << job.incidx << "]im^a";
         else if(job.intype == 8)
            os << "b[" << job.incmon << "," << job.incidx << "]im^a";
         else if(job.intype == 9)
            os << "c[" << job.incmon << "," << job.incidx << "]im^a";
         else if(job.intype == 10)
            os << "f[" << job.incmon << "," << job.incidx << "]im^b";
         else if(job.intype == 11)
            os << "b[" << job.incmon << "," << job.incidx << "]im^b";
         else if(job.intype == 12)
            os << "c[" << job.incmon << "," << job.incidx << "]im^b";
         else
            os << "invalid operation";
      }
   }
   else if(job.adtype == 5)
   {
      os << "b[" << job.updmon << "," << job.updidx << "]im += ";
      if(job.incmon < 0) // no distinction between -1 and -2 ...
         os << "cff[" << job.incidx << "]im";
      else
      {
         if(job.intype == 7)
            os << "f[" << job.incmon << "," << job.incidx << "]im^a";
         else if(job.intype == 8)
            os << "b[" << job.incmon << "," << job.incidx << "]im^a";
         else if(job.intype == 9)
            os << "c[" << job.incmon << "," << job.incidx << "]im^a";
         else if(job.intype == 10)
            os << "f[" << job.incmon << "," << job.incidx << "]im^b";
         else if(job.intype == 11)
            os << "b[" << job.incmon << "," << job.incidx << "]im^b";
         else if(job.intype == 12)
            os << "c[" << job.incmon << "," << job.incidx << "]im^b";
         else
            os << "invalid operation";
      }
   }
   else if(job.adtype == 6)
   {
      os << "c[" << job.updmon << "," << job.updidx << "]im += ";
      if(job.incmon < 0) // no distinction between -1 and -2 
         os << "cff[" << job.incidx << "]im";
      else
      {
         if(job.intype == 7)
            os << "f[" << job.incmon << "," << job.incidx << "]im^a";
         else if(job.intype == 8)
            os << "b[" << job.incmon << "," << job.incidx << "]im^a";
         else if(job.intype == 9)
            os << "c[" << job.incmon << "," << job.incidx << "]im^a";
         else if(job.intype == 10)
            os << "f[" << job.incmon << "," << job.incidx << "]im^b";
         else if(job.intype == 11)
            os << "b[" << job.incmon << "," << job.incidx << "]im^b";
         else if(job.intype == 12)
            os << "c[" << job.incmon << "," << job.incidx << "]im^b";
         else
            os << "invalid operation";
      }
   }
   else
      os << "invalid operation";

   return os;
}
