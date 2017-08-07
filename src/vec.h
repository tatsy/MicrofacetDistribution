#ifndef _VEC_H_
#define _VEC_H_

#include <stdexcept>

//! Vector class
class Vec {
public:
    Vec()
        : x{0.0}, y{0.0}, z{0.0} {
    }

    explicit Vec(double x_)
        : x{x_}, y{x_}, z{x_} {
    }

    Vec(double x_, double y_, double z_) 
        : x{x_}, y{y_}, z{z_} {
    }

    Vec(const Vec &v)
        : x{v.x}, y{v.y}, z{v.z} {
    }

    Vec &operator=(const Vec &v) {
        x = v.x, y = v.y, z = v.z;
        return *this;
    }

    Vec &operator+=(const Vec &b) { 
        x += b.x, y += b.y, z += b.z;
        return *this;
    }

    Vec &operator-=(const Vec &b) {
        return this->operator+=(-b);
    }

    Vec operator-() const {
        return Vec(-x, -y, -z);
    }

    Vec &operator*=(double s) {
        x *= s, y *= s, z *= s;
        return *this;
    }

    Vec &operator*=(const Vec &b) {
        x *= b.x, y *= b.y, z *= b.z;
        return *this;
    }

    Vec &operator/=(double s) {
        if (s == 0.0) throw std::runtime_error("Zero division!!");
        return this->operator*=(1.0 / s);
    }

    double length() const { return std::sqrt(this->dot(*this)); }
    Vec& normalize() { return this->operator/=(this->length()); }
    double dot(const Vec &b) const { return x*b.x + y*b.y + z*b.z; }
    Vec cross(const Vec &b) {return Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }

    double x, y, z;
};

Vec operator+(const Vec &a, const Vec &b) {
    Vec ret = a;
    ret += b;
    return ret;
}

Vec operator-(const Vec &a, const Vec &b) {
    Vec ret = a;
    ret -= b;
    return ret;
}

Vec operator*(const Vec &a, const Vec &b) {
    Vec ret = a;
    ret *= b;
    return ret;
}

Vec operator*(const Vec &a, double s) {
    Vec ret = a;
    ret *= s;
    return ret;
}

Vec operator*(double s, const Vec &a) {
    Vec ret = a;
    ret *= s;
    return ret;
}

Vec operator/(const Vec &a, double s) {
    Vec ret = a;
    ret /= s;
    return ret;
}

#endif  // _VEC_H_
