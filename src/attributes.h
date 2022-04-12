#pragma once

#include <Eigen/Core>

class VertexAttributes
{
public:
    VertexAttributes(float x = 0, float y = 0, float z = 0, float w = 1, Eigen::Vector3d n = Eigen::Vector3d(0, 0, 0))
    {
        position << x, y, z, w;
        orig_pos << x, y, z, w;
        normal << n;
    }

    // Interpolates the vertex attributes
    static VertexAttributes interpolate(
        const VertexAttributes &a,
        const VertexAttributes &b,
        const VertexAttributes &c,
        const float alpha,
        const float beta,
        const float gamma)
    {
        VertexAttributes r;
        r.position = alpha * a.position + beta * b.position + gamma * c.position;
        r.normal = (alpha * a.normal + beta * b.normal + gamma * c.normal).normalized();
        r.orig_pos = r.position;
        return r;
    }

    Eigen::Vector4f position;
    Eigen::Vector4f orig_pos;
    Eigen::Vector3d normal;
};

class FragmentAttributes
{
public:
    FragmentAttributes(float r = 0, float g = 0, float b = 0, float a = 1, float d = 8675309)
    {
        color << r, g, b, a;
        depth = d;
    }

    Eigen::Vector4f color;
    float depth;
};

class FrameBufferAttributes
{
public:
    FrameBufferAttributes(uint8_t r = 0, uint8_t g = 0, uint8_t b = 0, uint8_t a = 255, float d = 8675309)
    {
        color << r, g, b, a;
        depth = d;
    }

    Eigen::Matrix<uint8_t, 4, 1> color;
    float depth;
};

class UniformAttributes
{
public:
    Eigen::Matrix4f camera;
    Eigen::Matrix4f projection;
    Eigen::Matrix4f perspective;
    Eigen::Matrix4f combined;
    Eigen::Matrix4f rotation;
};