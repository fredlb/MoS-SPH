#pragma once

struct vec2
{
    float   x, y;

    vec2(): x(0.0), y(0.0){}

    vec2(float _x, float _y): x(_x), y(_y){}
    vec2(const vec2 &v): x(v.x), y(v.y){}
	
    inline vec2 operator - () const
    {
        return vec2(-x, -y);
    }

    inline vec2 & operator = (const vec2 &v)
    {
        x = v.x;
        y = v.y;
        return *this;
    }

    inline vec2 & operator += (const vec2 &v)
    {
        x += v.x;
        y += v.y;
        return *this;
    }

    inline vec2 & operator -= (const vec2 &v)
    {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    inline vec2 & operator *= (float scalar)
    {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    inline vec2 & operator /= (float scalar)
    {
        x /= scalar;
        y /= scalar;
        return *this;
    }
	
};

inline vec2 operator + (const vec2 &v1, const vec2 &v2)
{
    return vec2(v1.x + v2.x, v1.y + v2.y);
}

inline vec2 operator - (const vec2 &v1, const vec2 &v2)
{
    return vec2(v1.x - v2.x, v1.y - v2.y);
}

inline vec2  operator * (const vec2 &v, float scalar)
{
    return vec2(v.x * scalar, v.y * scalar);
}

inline vec2 operator * (float scalar, const vec2 &v)
{
    return vec2(scalar * v.x, scalar * v.y);
}

inline vec2 operator * (const vec2 &v1, const vec2 &v2)
{
    return vec2(v1.x * v2.x, v1.y * v2.y);
}

inline vec2 operator / (const vec2 &v, float scalar)
{
    return vec2(v.x / scalar, v.y / scalar);
}

inline vec2 operator / (float scalar, const vec2 &v)
{
    return vec2(scalar / v.x, scalar / v.y);
}

inline vec2 operator / (const vec2 &v1, const vec2 &v2)
{
    return vec2(v1.x / v2.x, v1.y / v2.y);
}

inline float dot(const vec2 &v1, const vec2 &v2)
{
	return v1.x*v2.x + v1.y*v2.y;
}
