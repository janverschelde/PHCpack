#ifndef DICT_H
#define DICT_H

#include <iostream>
#include <map>
#include <string>
#include <algorithm>

using namespace std;

// Dictionary class from variables to positions
/*
   Dictionary maps from variables to positions.
   The positions of new variable are their order, 
   when inserted into this dictionary.
*/

template <class T>
class Dict
{
   map<T, int> job; // Map from string to integer

   public:

      int n_job; // Number of jobs

      Dict(){ n_job = 0; }

      // Constructor builds an empty dictionary.
      ~Dict() // Destructor
      {
         // cout << "Dictionary destructed" << endl;
      }

      int at ( T job_key );
      // Get the value string
      /*
         Map jobkey to its value.
         @param job_key a string of job key
         @return integer of job value
       */

      bool has ( T job_key );
      // Check whether job is in the dictionary
      /*
         use map.find()
         @param job_key
         @return boolean varible of the result to search job
       */

      int insert ( T new_job );
      // Insert an string into dictionary
      /*
         New job value is n_job. Increase n_job by 1.
         @param new_job a string of job key
         @return void
       */

      void insert(T new_job, int val);

     int get ( T job_key );

     void print(); // Prints the jobs with their values.

     T* reverse();
};

class JobKey
{
   public:

      int n_job;
      int* job_begin;
      int last_level;
      JobKey ( int* job_begin, int n_job, int last_level )
      {
         this->job_begin = job_begin;
         this->n_job = n_job;
         this->last_level = last_level;
      };
      bool operator< ( const JobKey& job_key ) const
      {
         if(n_job<job_key.n_job)
         {
            return true;
         }
         else if(n_job>job_key.n_job)
         {
            return false;
         }
         if(last_level<job_key.last_level)
         {
            return true;
         }
         else if(last_level>job_key.last_level)
         {
            return false;
         }
         for(int i=0; i<n_job; i++)
         {
            if(job_begin[i] < job_key.job_begin[i])
            {
               return true;
            }
            else if(job_begin[i] > job_key.job_begin[i])
            {
               return false;
            }
        }
        return false;
    }

    friend ostream& operator<<(ostream& output, const JobKey& job_key);
};

inline ostream& operator<< ( ostream& output, const JobKey& job_key )
{
   output << job_key.n_job << " "
          << job_key.last_level << " ";
   for(int i=0; i<job_key.n_job; i++)
   {
      output << job_key.job_begin[i] << " ";
   }
   return output;  // for multiple << operators.
}

template <class T>
int Dict<T>::insert(T new_job)
{
   job[new_job]=n_job;
   n_job++;
   return n_job-1;
}

template <class T>
void Dict<T>::insert(T new_job, int val)
{
   job[new_job]=val;
   n_job++;
}

template <class T>
int Dict<T>::at(T job_key)
{
   return job[job_key];
}

template <class T>
bool Dict<T>::has(T job_key)
{
   if(job.find(job_key)== job.end())
      return false;
   else
      return true;
}

template <class T>
void Dict<T>::print()
{
   int i = 0;
   for(typename map<T, int>::const_iterator it = job.begin();
       it!=job.end(); ++it)
   {
      cout << i++ << " ";
      cout << it->first << " "<< it->second << endl;
   }
}

template <class T>
int Dict<T>::get ( T job_key )
{
   typename map<T, int>::const_iterator result = job.find(job_key);
   if(result== job.end())
   {
      job[job_key]=n_job;
      n_job++;
      return n_job-1;
   }
   else
   {
      return result->second;
   }
}

template <class T>
T* Dict<T>::reverse()
{
   T* pos_var = new T[n_job];
   for(typename map<T, int>::const_iterator it = job.begin();
       it!=job.end(); ++it)
   {
      pos_var[it->second] = it->first;
   }
   return pos_var;
}

typedef Dict<JobKey> JobDict;
typedef Dict<string> VarDict;

#endif
