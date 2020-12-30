// The file addition_job.cpp defines the methods of the class AdditionJob,
// specified in "addition_job.h".

#include "addition_job.h"

AdditionJob::AdditionJob
 ( int atp, int itp, int monix1, int monix2, int ix1, int ix2 )
{
   adtype = atp;
   intype = itp;
   updmon = monix1;
   incmon = monix2;
   updidx = ix1;
   incidx = ix2;
}

int AdditionJob::get_addition_type ( void ) const
{
   return adtype;
}

int AdditionJob::get_increment_type ( void ) const
{
   return intype;
}

int AdditionJob::get_update_monomial ( void ) const
{
   return updmon;
}

int AdditionJob::get_update_index ( void ) const
{
   return updidx;
}

int AdditionJob::get_increment_monomial ( void ) const
{
   return incmon;
}

int AdditionJob::get_increment_index ( void ) const
{
   return incidx;
}

std::ostream& operator<< ( std::ostream& os, const AdditionJob& job )
{
   if(job.adtype == 1)
   {
      os << "f[" << job.updmon << "," << job.updidx << "] += ";
      if(job.incmon < 0)
      {
         if(job.incidx < 0)
            os << "cst";
         else
            os << "cff[" << job.incidx << "]";
      }
      else
      {
         if(job.intype == 1)
            os << "f[" << job.incmon << "," << job.incidx << "]";
         else if(job.intype == 2)
            os << "b[" << job.incmon << "," << job.incidx << "]";
         else
            os << "c[" << job.incmon << "," << job.incidx << "]";
      }
   }
   else if(job.adtype == 2)
   {
      os << "b[" << job.updmon << "," << job.updidx << "] += ";
      if(job.incmon < 0)
         os << "cff[" << job.incidx << "]";
      else
      {
         if(job.intype == 1)
            os << "f[" << job.incmon << "," << job.incidx << "]";
         else if(job.intype == 2)
            os << "b[" << job.incmon << "," << job.incidx << "]";
         else
            os << "c[" << job.incmon << "," << job.incidx << "]";
      }
   }
   else
   {
      os << "c[" << job.updmon << "," << job.updidx << "] += ";
      if(job.incmon < 0)
         os << "cff[" << job.incidx << "]";
      else
      {
         if(job.intype == 1)
            os << "f[" << job.incmon << "," << job.incidx << "]";
         else if(job.intype == 2)
            os << "b[" << job.incmon << "," << job.incidx << "]";
         else
            os << "c[" << job.incmon << "," << job.incidx << "]";
      }
   }
   return os;
}
