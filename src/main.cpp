// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "raster.h"

#include <gif.h>
#include <fstream>

#include <Eigen/Geometry>
// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

using namespace std;
using namespace Eigen;

//Image height
const int H = 480;

//Camera settings
const double near_plane = 1.5; //AKA focal length
const double far_plane = near_plane * 100;
const double field_of_view = 0.7854; //45 degrees
const double aspect_ratio = 1.5;
const bool is_perspective = false;
const Vector3d camera_position(0, 0, 3);
const Vector3d camera_gaze(0, 0, -1);
const Vector3d camera_top(0, 1, 0);

//Object
const std::string data_dir = DATA_DIR;
const std::string mesh_filename(data_dir + "bunny.off");
MatrixXd vertices; // n x 3 matrix (n points)
MatrixXi facets;   // m x 3 matrix (m triangles)

//Material for the object
const Vector3d obj_diffuse_color(0.5, 0.5, 0.5);
const Vector3d obj_specular_color(0.2, 0.2, 0.2);
const double obj_specular_exponent = 256.0;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector3d> light_colors;
//Ambient light
const Vector3d ambient_light(0.3, 0.3, 0.3);

//Fills the different arrays
void setup_scene()
{
    //Loads file
    std::ifstream in(mesh_filename);
    if (!in.good())
    {
        std::cerr << "Invalid file " << mesh_filename << std::endl;
        exit(1);
    }
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i)
    {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i)
    {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16);
}

void build_uniform(UniformAttributes &uniform)
{
    //NOTE: setup uniform

    //NOTE: setup camera, compute w, u, v
    Matrix4f camera(4, 4);
    camera.setZero();
    Vector3f e = camera_position.cast <float> ();
    Vector3f g = camera_gaze.cast <float> ();
    Vector3f top = camera_top.cast <float> ();

    Vector3f w = -g.normalized();
    Vector3f u = (top.cross(w)).normalized();
    Vector3f v = w.cross(u);

    //NOTE: compute the camera transformation
    camera.block(0, 0, 3, 1) = u;
    camera.block(0, 1, 3, 1) = v;
    camera.block(0, 2, 3, 1) = w;
    camera.block(0, 3, 3, 1) = e;
    camera(3, 3) = 1.0;
    uniform.camera = camera.inverse();

    // NOTE: setup projection matrix
    float image_y = near_plane * tan(field_of_view / 2);
    float image_x = aspect_ratio * image_y;

    float l = -image_x;
    float b = -image_y;
    float n = -near_plane;

    float r = image_x;
    float t = image_y;
    float f = -far_plane;

    Matrix4f projection;
    projection.setZero();
    projection(0, 0) = 2 / (r - l);
    projection(1, 1) = 2 / (t - b);
    projection(2, 2) = 2 / (n - f);
    projection(3, 3) = 1;
    projection(0, 3) = -(r + l) / (r - l);
    projection(1, 3) = -(t + b) / (t - b);
    projection(2, 3) = -(n + f) / (n - f);
    uniform.projection = projection;

    Matrix4f P;
    if (is_perspective)
    {
        //TODO setup perspective camera


        uniform.combined = uniform.projection * uniform.perspective * uniform.camera;
    }
    else
    {
        uniform.combined = uniform.projection * uniform.camera;
    }
}

void simple_render(Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //NOTE: fill the shader
        VertexAttributes out;
        out.position = uniform.combined * va.position;
        return out;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //NOTE: fill the shader
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //NOTE: fill the shader
        return FrameBufferAttributes(fa.color[0] * 255, fa.color[1] * 255, fa.color[2] * 255, fa.color[3] * 255);
    };

    std::vector<VertexAttributes> vertex_attributes;
    //NOTE: build the vertex attributes from vertices and facets
    for (int i = 0; i < facets.rows(); i++)
    {
        Vector3d a = vertices.row(facets(i, 0));
        Vector3d b = vertices.row(facets(i, 1));
        Vector3d c = vertices.row(facets(i, 2));

        vertex_attributes.push_back(VertexAttributes(a[0], a[1], a[2]));
        vertex_attributes.push_back(VertexAttributes(b[0], b[1], b[2]));
        vertex_attributes.push_back(VertexAttributes(c[0], c[1], c[2]));
    }

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

Matrix4f compute_rotation(const double alpha)
{
    //NOTE: Compute the rotation matrix of angle alpha on the y axis around the object barycenter
    Matrix4f res;

    // find barycenter by taking average of all triangle barycenters
    Vector3f center(0, 0, 0);
    for (int i = 0; i < facets.rows(); i++)
    {
        Vector3f a = vertices.row(facets(i, 0)).cast <float> ();
        Vector3f b = vertices.row(facets(i, 1)).cast <float> ();
        Vector3f c = vertices.row(facets(i, 2)).cast <float> ();
        center += (a + b + c) / 3;
    }
    center /= facets.rows();

    Matrix4f translate_to_origin;
    translate_to_origin.setZero();
    translate_to_origin(0, 0) = 1;
    translate_to_origin(1, 1) = 1;
    translate_to_origin(2, 2) = 1;
    translate_to_origin(3, 3) = 1;
    translate_to_origin.block(0, 3, 3, 1) = -center;

    Matrix4f rotation;
    double angle = alpha * 2 * 3.1415;
    rotation.setZero();
    rotation(0, 0) = cos(angle);
    rotation(0, 2) = -sin(angle);
    rotation(2, 0) = sin(angle);
    rotation(2, 2) = cos(angle);
    rotation(1, 1) = 1;
    rotation(3, 3) = 1;

    Matrix4f translate_from_origin;
    translate_from_origin.setZero();
    translate_from_origin(0, 0) = 1;
    translate_from_origin(1, 1) = 1;
    translate_from_origin(2, 2) = 1;
    translate_from_origin(3, 3) = 1;
    translate_from_origin.block(0, 3, 3, 1) = center;

    res = translate_from_origin * rotation * translate_to_origin;

    return res;
}

void wireframe_render(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    uniform.rotation = compute_rotation(alpha);

    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        // NOTE: fill the shader
        VertexAttributes out;
        out.position = uniform.combined * uniform.rotation * va.position;
        return out;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //NOTE: fill the shader
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //NOTE: fill the shader
        return FrameBufferAttributes(fa.color[0] * 255, fa.color[1] * 255, fa.color[2] * 255, fa.color[3] * 255);
    };

    std::vector<VertexAttributes> vertex_attributes;

    //NOTE: generate the vertex attributes for the edges and rasterize the lines
    //NOTE: use the transformation matrix
    for (int i = 0; i < facets.rows(); i++)
    {
        Vector3d a = vertices.row(facets(i, 0));
        Vector3d b = vertices.row(facets(i, 1));
        Vector3d c = vertices.row(facets(i, 2));

        vertex_attributes.push_back(VertexAttributes(a[0], a[1], a[2]));
        vertex_attributes.push_back(VertexAttributes(b[0], b[1], b[2]));
        vertex_attributes.push_back(VertexAttributes(b[0], b[1], b[2]));
        vertex_attributes.push_back(VertexAttributes(c[0], c[1], c[2]));
        vertex_attributes.push_back(VertexAttributes(c[0], c[1], c[2]));
        vertex_attributes.push_back(VertexAttributes(a[0], a[1], a[2]));
    }

    rasterize_lines(program, uniform, vertex_attributes, 0.5, frameBuffer);
}

void get_shading_program(Program &program)
{
    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //HACK: transform the position and the normal
        VertexAttributes out;
        out.position = uniform.combined * uniform.rotation * va.position;
        out.normal = va.normal;
        return out;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //NOTE: compute the correct lighting
        const Vector3d pos = Vector3d(va.position[0], va.position[1], va.position[2]);
        Vector3d v = (camera_position - pos).normalized();

        Vector3d lights_color(0, 0, 0);
        for (int i = 0; i < light_positions.size(); i++)
        {
            const Vector3d &light_position = light_positions[i];
            const Vector3d &light_color = light_colors[i];

            const Vector3d pos = Vector3d(va.position[0], va.position[1], va.position[2]);
            const Vector3d Li = (light_position - pos).normalized();
            
            // Diffuse contribution
            const Vector3d diffuse = obj_diffuse_color * std::max(Li.dot(va.normal), 0.0);

            // Specular contribution, use obj_specular_color
            Vector3d h = (v + Li).normalized();
            const Vector3d specular = obj_specular_color * std::pow(std::max(h.dot(va.normal), 0.0), obj_specular_exponent);

            // Attenuate lights according to the squared distance to the lights
            const Vector3d D = light_position - pos;
            lights_color += (diffuse + specular).cwiseProduct(light_color) / D.squaredNorm();
        }

        lights_color += ambient_light;
        return FragmentAttributes(lights_color[0], lights_color[1], lights_color[2], 1, va.position[2]);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //NOTE: implement the depth check
        if (fa.depth > previous.depth)
        {
            return FrameBufferAttributes(fa.color[0] * 255, fa.color[1] * 255, fa.color[2] * 255, fa.color[3] * 255, fa.depth);
        }
        else
        {
            return previous;
        }
    };
}

void flat_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    uniform.rotation = compute_rotation(alpha);
    Program program;
    get_shading_program(program);

    std::vector<VertexAttributes> vertex_attributes;
    //NOTE: compute the normals
    for (int i = 0; i < facets.rows(); i++)
    {
        Vector3d a = vertices.row(facets(i, 0));
        Vector3d b = vertices.row(facets(i, 1));
        Vector3d c = vertices.row(facets(i, 2));

        Vector4f a4d(a[0], a[1], a[2], 1);
        Vector4f b4d(b[0], b[1], b[2], 1);
        Vector4f c4d(c[0], c[1], c[2], 1);

        a4d = uniform.rotation * a4d;
        b4d = uniform.rotation * b4d;
        c4d = uniform.rotation * c4d;

        Vector3d a_rot(a4d[0], a4d[1], a4d[2]);
        Vector3d b_rot(b4d[0], b4d[1], b4d[2]);
        Vector3d c_rot(c4d[0], c4d[1], c4d[2]);

        Vector3d normal = ((b_rot - a_rot).cross(c_rot - a_rot)).normalized();

        vertex_attributes.push_back(VertexAttributes(a[0], a[1], a[2], 1, normal));
        vertex_attributes.push_back(VertexAttributes(b[0], b[1], b[2], 1, normal));
        vertex_attributes.push_back(VertexAttributes(c[0], c[1], c[2], 1, normal));
    }

    //NOTE: set material colors

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

void pv_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    uniform.rotation = compute_rotation(alpha);
    Program program;
    get_shading_program(program);


    //NOTE: compute the vertex normals as vertex normal average
    std::vector<Vector3d> face_normals;
    std::vector<Vector3d> vertex_normals(vertices.size());
    for (int i = 0; i < facets.rows(); i++)
    {
        Vector3d a = vertices.row(facets(i, 0));
        Vector3d b = vertices.row(facets(i, 1));
        Vector3d c = vertices.row(facets(i, 2));

        Vector4f a4d(a[0], a[1], a[2], 1);
        Vector4f b4d(b[0], b[1], b[2], 1);
        Vector4f c4d(c[0], c[1], c[2], 1);

        a4d = uniform.rotation * a4d;
        b4d = uniform.rotation * b4d;
        c4d = uniform.rotation * c4d;

        Vector3d a_rot(a4d[0], a4d[1], a4d[2]);
        Vector3d b_rot(b4d[0], b4d[1], b4d[2]);
        Vector3d c_rot(c4d[0], c4d[1], c4d[2]);

        Vector3d normal = ((b_rot - a_rot).cross(c_rot - a_rot)).normalized();

        face_normals.push_back(normal);

        vertex_normals[facets(i, 0)] += normal;
        vertex_normals[facets(i, 1)] += normal;
        vertex_normals[facets(i, 2)] += normal;
    }

    for (int i = 0; i < vertex_normals.size(); i++)
    {
        vertex_normals[i] = vertex_normals[i].normalized();
    }

    std::vector<VertexAttributes> vertex_attributes;
    //NOTE: create vertex attributes
    //NOTE: set material colors

    for (int i = 0; i < facets.rows(); i++)
    {
        Vector3d a = vertices.row(facets(i, 0));
        Vector3d b = vertices.row(facets(i, 1));
        Vector3d c = vertices.row(facets(i, 2));

        vertex_attributes.push_back(VertexAttributes(a[0], a[1], a[2], 1, vertex_normals[facets(i, 0)]));
        vertex_attributes.push_back(VertexAttributes(b[0], b[1], b[2], 1, vertex_normals[facets(i, 1)]));
        vertex_attributes.push_back(VertexAttributes(c[0], c[1], c[2], 1, vertex_normals[facets(i, 2)]));
    }

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

int main(int argc, char *argv[])
{
    setup_scene();

    int W = H * aspect_ratio;
    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> frameBuffer(W, H);
    vector<uint8_t> image;

    simple_render(frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("simple.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    frameBuffer = Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> (W, H);
    wireframe_render(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("wireframe.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    frameBuffer = Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> (W, H);
    flat_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("flat_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    frameBuffer = Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> (W, H);
    pv_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("pv_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    //TODO: add the animation
    int delay = 25;
    GifWriter g;

    frameBuffer = Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> (W, H);
    const char* wireframe_file = "wireframe.gif";
    GifBegin(&g, wireframe_file, frameBuffer.rows(), frameBuffer.cols(), delay);
    for (float i = 1; i > 0; i -= 0.05)
    {
        frameBuffer.setConstant(FrameBufferAttributes());
        wireframe_render(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
    }
    GifEnd(&g);
    
    frameBuffer = Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> (W, H);
    const char* flat_file = "flat_shading.gif";
    GifBegin(&g, flat_file, frameBuffer.rows(), frameBuffer.cols(), delay);
    for (float i = 1; i > 0; i -= 0.05)
    {
        frameBuffer.setConstant(FrameBufferAttributes());
        flat_shading(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
    }
    GifEnd(&g);

    frameBuffer = Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> (W, H);
    const char* pv_file = "pv_shading.gif";
    GifBegin(&g, pv_file, frameBuffer.rows(), frameBuffer.cols(), delay);
    for (float i = 1; i > 0; i -= 0.05)
    {
        frameBuffer.setConstant(FrameBufferAttributes());
        pv_shading(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
    }
    GifEnd(&g);

    return 0;
}
