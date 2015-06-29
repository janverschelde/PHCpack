template <class T>
complexH<T>::complexH( const complexH<T>& source)
{
   real = source.real;
   imag = source.imag;
}

template <class T>
std::ostream& operator<< ( std::ostream& os, const complexH<T>& number )
{
   os << number.real << " + i*" << number.imag << endl;

   return os;
}

template <class T>
std:: ifstream& operator >> ( std::ifstream& is, complexH<T>& number)
{
   string a,b;
   is >> a >> b;
   char c[100];
   c[0] = b[0];
   is.ignore(10,' ');
   is.get(c+1, 99, 'i');
   std::cout << c << std::endl;
   number = complexH<T>(a.c_str(),c);
   is.get();

   return is;
}

template <class T>
complexH<T> complexH<T>::operator=(const complexH<T>& a)
{
   real = a.real;
   imag = a.imag;

   return *this;
}

template <class T>
complexH<T> complexH<T>::operator=(const int& a)
{
   real = a;
   imag = 0.0;

   return *this;
}

template <class T>
complexH<T> complexH<T>::operator=(const double& a)
{
   real = a;
   imag = 0.0;

   return *this;
}

template <class T>
complexH<T> complexH<T>::operator=(const dd_real& a)
{
   real = a;
   imag = 0.0;

   return *this;
}

template <class T>
complexH<T> complexH<T>::operator=(const qd_real& a)
{
   real = a;
   imag = 0.0;

   return *this;
}

template<class T>
inline complexH<T>::complexH(const double& a, const double& b )
{
   real = a;
   imag = b;
}

template<class T>
inline complexH<T>::complexH(const dd_real& a, const dd_real& b )
{
   real = a;
   imag = b;
}

template<class T>
inline complexH<T>::complexH(const qd_real& a, const qd_real& b )
{
   real = a;
   imag = b;
}

template<>
inline complexH<qd_real>::complexH(const char* a)
{
   real = qd_real(a);
   imag = 0.0;
}

template<>
inline complexH<dd_real>::complexH(const char* a)
{
   real = dd_real(a);
   imag = 0.0;
}

template<>
inline complexH<double>::complexH(const char* a)
{
   real = atof(a);
   imag = 0.0;
}

template<>
inline complexH<qd_real>::complexH(const char* a, const char* b)
{
   real = qd_real(a);
   imag = qd_real(b);
}

template<>
inline complexH<dd_real>::complexH(const char* a, const char* b)
{
   real = dd_real(a);
   imag = dd_real(b);
}

template<>
inline complexH<double>::complexH(const char* a, const char* b)
{
   real = atof(a);
   imag = atof(b);
}

template<class T>
inline complexH<T>::complexH(const int& a)
{
   real = a;
   imag = 0.0;
}

template<class T>
inline complexH<T>::complexH(const double& a)
{
   real = a;
   imag = 0.0;
}

template <class T>
void complexH<T>::init(const double& a, const double& b)
{
   complexH<T> temp(a,b);

   real = temp.real;
   imag = temp.imag;
}

template <class T>
void complexH<T>::init(const dd_real& a, const dd_real& b)
{
   real = a;
   imag = b;
}

template <class T>
void complexH<T>::init(const qd_real& a, const qd_real& b)
{
   real = a;
   imag = b;
}


/*template <class T>
__device__ complexH<T> complexH<T>::operator+(complexH<T> a)
{
   return complexH(real+a.real,imag+a.imag,1);
}*/

template <class T>
complexH<T> complexH<T>::operator+(const complexH<T>& a)
{
   return complexH(real+a.real,imag+a.imag);
}

template <class T>
void complexH<T>::operator+=(const complexH<T>& a)
{
   real += a.real;
   imag += a.imag;
}

template <class T>
void complexH<T>::operator-=(const complexH<T>& a)
{
   real -= a.real;
   imag -= a.imag;
}

template <class T>
complexH<T> complexH<T>::operator-(const complexH<T>& a)
{
   return complexH(real-a.real,imag-a.imag);
}

template <class T>
complexH<T> complexH<T>::operator*(const complexH<T>& a) const
{
   return complexH(real*a.real-imag*a.imag,imag*a.real+real*a.imag);
}

template <class T>
void complexH<T>::operator*=(const complexH<T>& a)
{
   T real_tmp = real;
   real = real*a.real-imag*a.imag;
   imag = imag*a.real+real_tmp*a.imag;
}

template <class T>
complexH<T> complexH<T>::operator*(const T& a)
{
   return complexH(real*a,imag*a);
}

template <class T>
complexH<T> complexH<T>::operator*(const int& a)
{
   return complexH(real*a,imag*a);
}

template <class T>
complexH<T> complexH<T>::operator/(const T& a)
{
   return complexH(real/a,imag/a);
}

template <class T>
complexH<T> complexH<T>::operator/(const complexH<T>& a)
{
   return complexH((real*a.real+imag*a.imag)/(a.real*a.real+a.imag*a.imag),
                   (imag*a.real-real*a.imag)/(a.real*a.real+a.imag*a.imag));
}
