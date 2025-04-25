//
// Created by William Zhao on 4/4/25.
//

#ifndef MATH_H
#define MATH_H

#include <cmath>

struct Vector3
{
    float x, y, z;
    Vector3() : x(0), y(0), z(0) {}
    Vector3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
    Vector3 operator+(const Vector3& other) const
    {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }
    bool operator==(const Vector3& other) const
    {
        return (x == other.x && y == other.y && z == other.z);
    }
    Vector3 operator-(const Vector3& other) const
    {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }
    Vector3 operator*(float scalar) const
    {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }
    Vector3 operator/(float scalar) const
    {
        return Vector3(x / scalar, y / scalar, z / scalar);
    }

    static float Dot(const Vector3& a, const Vector3& b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
    static Vector3 Cross(const Vector3& a, const Vector3& b)
    {
        return Vector3(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
        );
    }
    float Length() const
    {
        return std::sqrt(x * x + y * y + z * z);
    }
    float LengthSquared() const
    {
        return x * x + y * y + z * z;
    }
    void Normalize()
    {
        const float length = Length();
        if (length > 0)
        {
            x /= length;
            y /= length;
            z /= length;
        }
    }

    static void Normalize(Vector3& vect);

    static const Vector3 UnitX;
    static const Vector3 UnitY;
    static const Vector3 UnitZ;

    // Method to return a pointer to the first element
    const float* data() const {
        return &x;
    }
};

struct Vector4
{
    float x, y, z, w;
    Vector4() : x(0), y(0), z(0), w(0) {}
    Vector4(float x_, float y_, float z_, float w_) : x(x_), y(y_), z(z_), w(w_) {}
    Vector4 operator+(const Vector4& other) const
    {
        return Vector4(x + other.x, y + other.y, z + other.z, w + other.w);
    }
    Vector4 operator-(const Vector4& other) const
    {
        return Vector4(x - other.x, y - other.y, z - other.z, w - other.w);
    }
    Vector4 operator*(float scalar) const
    {
        return Vector4(x * scalar, y * scalar, z * scalar, w * scalar);
    }
    Vector4 operator/(float scalar) const
    {
        return Vector4(x / scalar, y / scalar, z / scalar, w / scalar);
    }

    static float Dot(const Vector4& a, const Vector4& b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
    }
    static Vector4 Cross(const Vector4& a, const Vector4& b, const Vector4& c)
    {
        return Vector4(
            a.y * (b.z * c.w - b.w * c.z) - a.z * (b.y * c.w - b.w * c.y) + a.w * (b.y * c.z - b.z * c.y),
            -(a.x * (b.z * c.w - b.w * c.z) - a.z * (b.x * c.w - b.w * c.x) + a.w * (b.x * c.z - b.z * c.x)),
            a.x * (b.y * c.w - b.w * c.y) - a.y * (b.x * c.w - b.w * c.x) + a.w * (b.x * c.y - b.y * c.x),
            -(a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) + a.z * (b.x * c.y - b.y * c.x))
        );
    }
    float Length() const
    {
        return std::sqrt(x * x + y * y + z * z + w * w);
    }
    float LengthSquared() const
    {
        return x * x + y * y + z * z + w * w;
    }
    void Normalize()
    {
        const float length = Length();
        if (length > 0)
        {
            x /= length;
            y /= length;
            z /= length;
            w /= length;
        }
    }

    static const Vector4 UnitX;
    static const Vector4 UnitY;
    static const Vector4 UnitZ;
    static const Vector4 UnitW;

    // Method to return a pointer to the first element
    const float* data() const {
        return &x;
    }
};

inline Vector3 operator*(float scalar, const Vector3& vec)
{
    return Vector3(vec.x * scalar, vec.y * scalar, vec.z * scalar);
}

inline Vector4 operator*(float scalar, const Vector4& vec)
{
    return Vector4(vec.x * scalar, vec.y * scalar, vec.z * scalar, vec.w * scalar);
}

namespace Math {
    const static float Pi = 3.14159265358979323846f;
    const static float TwoPi = 6.28318530717958647692f;
}

struct Ray
{
    Vector3 p,d;
    Ray(Vector3 pVect, Vector3 dVect)
    {
        pVect.Normalize();
        dVect.Normalize();
        p = pVect;
        d = dVect;
    }

    Ray(Vector3 dVect)
    : p()
    {
        dVect.Normalize();
        d = dVect;
    }

    Ray() 
    : p(), d() 
    {
    }

    Vector3 at(float t) const {
        return p + t * d;
    }
};

#endif //MATH_H
