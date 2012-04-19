/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _vec_H
#define	_vec_H

#include <iostream>
#include <math.h>
#include <stdexcept>
#include <string>
//#include "tokenizer.h"

namespace votca { namespace tools {
using namespace std;
/**
    \brief Vector class for a 3 component vector

    This class represents a 3 component vector to store e.g. postitions, velocities, forces, ...
    Operators for basic vector-vector and vector-scalar operations are defined.
    you can access the elements with the functions x(), y(), z(), both reading and writing is possible;
    x + v.x();
    v.x() = 5.;
*/

class vec {
public:
    
    vec();
    vec(const vec &v);
    vec(const double r[3]);
    vec(const double &x, const double &y, const double &z);
    
    vec &operator=(const vec &v);
    vec &operator+=(const vec &v);
    vec &operator-=(const vec &v);
    vec &operator*=(const double &d);
    vec &operator/=(const double &d);
    
    /**
     * \brief get full access to x element
     * @return reference to x
     */
    double &x() { return _x; }
    /**
     * \brief get full access to y element
     * @return reference to y
     */
    double &y() { return _y; }
    /**
     * \brief get full access to z element
     * @return reference to z
     */
    double &z() { return _z; }
    
    void setX(const double &x) { _x = x; }
    void setY(const double &y) { _y = y; }
    void setZ(const double &z) { _z = z; }
    
    /**
     * \brief read only access to x element
     * @return x const reference to x
     *
     * This function can be usefule when const is used to allow for better
     * optimization. Always use getX() instead of x() if possible.
     */
    const double &getX() const { return _x; }
    /**
     * \brief read only access to y element
     * @return x const reference to y
     *
     * This function can be usefule when const is used to allow for better
     * optimization. Always use getY() instead of y() if possible.
     */
    const double &getY() const { return _y; }
    /**
     * \brief read only access to z element
     * @return x const reference to z
     *
     * This function can be usefule when const is used to allow for better
     * optimization. Always use getZ() instead of Z() if possible.
     */
    const double &getZ() const { return _z; }
    
    /**
     * \brief normalize the vector
     * @return normalized vector
     * This function normalizes the vector and returns itself after normalization.
     * After this call, the vector stores the normalized value.
     */
    vec &normalize();
    
    private:
        double _x, _y, _z;
};

inline vec::vec() {}

inline vec::vec(const vec &v)
    : _x(v._x), _y(v._y), _z(v._z) {}
        
inline vec::vec(const double r[3])
    : _x(r[0]), _y(r[1]), _z(r[2]) {}
    
inline vec::vec(const double &x, const double &y, const double &z)
        : _x(x), _y(y), _z(z) {}
    
inline bool operator==(const vec &v1, const vec &v2)
{
    return ((v1.getX()==v2.getX()) && (v1.getY()==v2.getY()) && (v1.getZ()==v2.getZ()));
}

inline bool operator!=(const vec &v1, const vec &v2)
{
    return ((v1.getX()!=v2.getX()) || (v1.getY()!=v2.getY()) || (v1.getZ()==v2.getZ()));
}

inline vec &vec::operator=(const vec &v)
{ 
        _x=v._x; _y=v._y; _z=v._z;
        return *this;
}    

inline vec &vec::operator+=(const vec &v)
{ 
        _x+=v._x; _y+=v._y; _z+=v._z;
        return *this;
}    
        
inline vec &vec::operator-=(const vec &v)
{ 
        _x-=v._x; _y-=v._y; _z-=v._z;
        return *this;
}    

inline vec &vec::operator*=(const double &d)
{ 
        _x*=d; _y*=d; _z*=d;
        return *this;
}    

inline vec &vec::operator/=(const double &d)
{ 
        _x/=d; _y/=d; _z/=d;
        return *this;
}    

inline vec operator+(const vec &v1, const vec &v2)
{
    return (vec(v1)+=v2);
}

inline vec operator-(const vec &v1, const vec &v2)
{
    return (vec(v1)-=v2);
}

inline vec operator-(const vec &v1){
    return vec (-v1.getX(), -v1.getY(), -v1.getZ());
}

inline vec operator*(const vec &v1, const double &d)
{
    return (vec(v1)*=d);
}

inline vec operator*(const double &d, const vec &v1)
{
    return (vec(v1)*=d);
}

inline vec operator/(const vec &v1, const double &d)
{
    return (vec(v1)/=d);
}

inline std::ostream &operator<<(std::ostream &out, const vec& v)
{
      out << '[' << v.getX() << ", " << v.getY() << ", " << v.getZ() << ']';
      return out;
}

/*inline std::istream &operator>>(std::istream &in, vec& v)
{
    char c;
    in.get(c);
    if(c != '[') {
        throw std::runtime_error("error, invalid character in vector string");
    }
    
    string str;
    while (in.good()) {
        in.get(c);
        if(c==']') { // found end of vector
            Tokenizer tok(str, ",");
            vector<double> d;
            tok.ConvertToVector(d);
            if(d.size() != 3)
                throw std::runtime_error("error, invalid number of entries in vector");
            v.setX(d[0]);
            v.setY(d[1]);
            v.setZ(d[2]);
	    return in;
        }
        str += c;
    }
    throw std::runtime_error("did not find closing bracket in string to vec conversion");

    return in;
}*/
    
/// dot product
inline double operator*(const vec &v1, const vec &v2)
{
    return v1.getX()*v2.getX() + v1.getY()*v2.getY() + v1.getZ()*v2.getZ();
}

/// cross product
inline vec operator^(const vec &v1, const vec &v2)
{
    return vec(
        v1.getY()*v2.getZ() - v1.getZ()*v2.getY(),
        v1.getZ()*v2.getX() - v1.getX()*v2.getZ(),
        v1.getX()*v2.getY() - v1.getY()*v2.getX()
    );
}

inline double abs(const vec &v)
{
    return sqrt(v*v);
}

inline vec &vec::normalize()
{ 
    return ((*this)*=1./abs(*this));
}

}}

#endif	/* _vec_H */

