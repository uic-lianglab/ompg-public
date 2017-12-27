#ifndef POINT_H
#define POINT_H

#include <math.h>

/**
 * 3-D coordinate
 */
class Point {
public:
    /**
     * 3-D (x,y,z) coordinates
     */
    double x,y,z;

    /**
     * Default constructor
     */
    inline Point() : x(0), y(0), z(0) {}

    /**
     * Constructor initializes 3-D coordinates
     */
    inline Point(const double X,
                 const double Y,
                 const double Z): x(X), y(Y), z(Z) {}

    /**
     * Add two points
     */
    inline Point operator+(const Point &p) const {
        return Point(x+p.x,y+p.y,z+p.z);
    }

    /**
     * In-place add two points
     */
    inline Point& operator+=(const Point &p) {
        x += p.x; y += p.y; z += p.z;
        return *this;
    }

    /**
     * Subtract two points
     */
    inline Point operator-(const Point &p) const {
        return Point(x-p.x,y-p.y,z-p.z);
    }

    /**
     * In-place subtract two points
     */
    inline Point& operator-=(const Point &p) {
        x -= p.x; y -= p.y; z -= p.z;
        return *this;
    }
   
    /**
     * Divide by scalar
     */
    inline Point operator/(const double d) const {
        return Point(x / d, y / d, z / d);
    }

    /**
     * In-place divide by scalar
     */
    inline Point& operator/=(const double d) {
        x /= d; y /= d; z /= d;
        return *this;
    }

    /**
     * Multiply by scalar
     */
    inline Point operator*(const double d) const {
        return Point(x*d, y*d, z*d);
    }

    /**
     * In-place multiply by scalar
     */
    inline Point& operator*=(const double d) {
        x *= d; y *= d; z *= -d;
        return *this;
    }

    /**
     * @return sum of squared coordinates
     */
    inline double square() const {
        return x*x+y*y+z*z;
    }

    /**
     * @return cross-product
     */
    inline Point cross(const Point &v) const {
        return Point(y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x);
    }

    /**
     * @return dot-product
     */
    inline double dot(const Point& v) const {
        return x*v.x + y*v.y + z*v.z;
    }
    
    /**
     * @return distance to parameter point
     */
    inline double dis(const Point& v) const {
        return sqrt(disquare(v));
    }

    /**
     * @return squared distance to parameter point
     */
    inline double disquare(const Point& v) const {
        return ((x - v.x)*(x - v.x) + (y - v.y)*(y - v.y) + (z - v.z)*(z - v.z));
    }
};

#endif // POINT_H
