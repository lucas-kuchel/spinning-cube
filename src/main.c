#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.14159265358979323846

typedef struct vec4 {
    union {
        struct {
            double x;
            double y;
            double z;
            double w;
        };

        double values[4];
    };
} vec4;

typedef struct vec2 {
    union {
        struct {
            double x;
            double y;
        };

        double values[2];
    };
} vec2;

typedef struct uvec3 {
    union {
        struct {
            size_t x;
            size_t y;
            size_t z;
        };

        size_t values[3];
    };
} uvec3;

typedef struct triangle {
    union {
        struct {
            vec2 a;
            vec2 b;
            vec2 c;
        };

        vec2 values[3];
    };
} triangle;

typedef struct mat4x4 {
    vec4 values[4];
} mat4x4;

double vec4_dot(const vec4* a, const vec4* b) {
    return (a->x * b->x) + (a->y * b->y) + (a->z * b->z) + (a->w * b->w);
}

vec2 vec2_sub(const vec2* a, const vec2* b) {
    vec2 result;

    result.x = a->x - b->x;
    result.y = a->y - b->y;

    return result;
}

double vec2_dot(const vec2* a, const vec2* b) {
    return (a->x * b->x) + (a->y * b->y);
}

vec2 vec2_rdiag(const vec2* a) {
    vec2 result;

    result.x = a->y;
    result.y = -a->x;

    return result;
}

vec2 vec2_ldiag(const vec2* a) {
    vec2 result;

    result.x = -a->y;
    result.y = a->x;

    return result;
}

bool triangle_test(const triangle* tri, const vec2* point) {
    vec2 ab_edge = vec2_sub(&tri->b, &tri->a);
    vec2 ca_edge = vec2_sub(&tri->a, &tri->c);
    vec2 bc_edge = vec2_sub(&tri->c, &tri->b);

    vec2 ab_diag = vec2_rdiag(&ab_edge);
    vec2 ca_diag = vec2_rdiag(&ca_edge);
    vec2 bc_diag = vec2_rdiag(&bc_edge);

    vec2 a_pdist = vec2_sub(point, &tri->a);
    vec2 c_pdist = vec2_sub(point, &tri->c);
    vec2 b_pdist = vec2_sub(point, &tri->b);

    double a_dist = vec2_dot(&ab_diag, &a_pdist);
    double b_dist = vec2_dot(&ca_diag, &c_pdist);
    double c_dist = vec2_dot(&bc_diag, &b_pdist);

    return (a_dist <= 0.0 && b_dist <= 0.0 && c_dist <= 0.0);
}

vec4 row(const mat4x4* mat, size_t n) {
    vec4 result;

    for (size_t i = 0; i < 4; i++) {
        result.values[i] = mat->values[i].values[n];
    }

    return result;
}

mat4x4 transpose(const mat4x4* mat) {
    mat4x4 result;

    for (size_t i = 0; i < 4; i++) {
        result.values[i] = row(mat, i);
    }

    return result;
}

mat4x4 multiply(const mat4x4* a, const mat4x4* b) {
    mat4x4 ax = transpose(a);
    mat4x4 result;

    for (size_t i = 0; i < 4; i++) {
        const vec4* col = &b->values[i];

        for (size_t j = 0; j < 4; j++) {
            const vec4* row = &ax.values[j];

            result.values[i].values[j] = vec4_dot(row, col);
        }
    }

    return result;
}

vec4 transform(const mat4x4* a, const vec4* b) {
    mat4x4 ax = transpose(a);
    vec4 result;

    for (size_t i = 0; i < 4; i++) {
        const vec4* row = &ax.values[i];

        result.values[i] = vec4_dot(row, b);
    }

    return result;
}

mat4x4 identity() {
    mat4x4 result;

    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            result.values[i].values[j] = (i == j);
        }
    }

    return result;
}

mat4x4 empty() {
    mat4x4 result;

    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            result.values[i].values[j] = 0;
        }
    }

    return result;
}

mat4x4 rotation_x(double rad) {
    mat4x4 result = identity();

    result.values[1].y = cos(rad);
    result.values[1].z = -sin(rad);

    result.values[2].y = sin(rad);
    result.values[2].z = cos(rad);

    return result;
}

mat4x4 rotation_y(double rad) {
    mat4x4 result = identity();

    result.values[0].x = cos(rad);
    result.values[0].z = -sin(rad);

    result.values[2].x = sin(rad);
    result.values[2].z = cos(rad);

    return result;
}

mat4x4 rotation_z(double rad) {
    mat4x4 result = identity();

    result.values[0].x = cos(rad);
    result.values[0].y = -sin(rad);

    result.values[1].x = sin(rad);
    result.values[1].y = cos(rad);

    return result;
}

void clear(char* screen, size_t width, size_t height) {
    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            screen[y * width + x] = ' ';
        }
    }

    printf("\x1b[H");
}

void present(char* screen, size_t width, size_t height) {
    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            putchar(screen[y * width + x]);
        }

        putchar('\n');
    }
}

size_t ndc_axis_to_pixel(double axis, size_t res) {
    if (axis < -1.0) {
        axis = -1.0;
    }

    if (axis > 1.0) {
        axis = 1.0;
    }

    double normalized = (axis + 1.0) * 0.5;
    double pixel = normalized * res;

    return (size_t)round(pixel);
}

double pixel_axis_to_ndc(size_t axis, size_t res) {
    double daxis = (double)axis;
    double dres = (double)res;

    return (2 * (daxis / dres)) - 1;
}

void set(vec2* point, char value, char* screen, size_t width, size_t height) {
    double x = point->x;
    double y = point->y;

    if (x > 1.0 || x < -1.0 || y > 1.0 || y < -1.0) {
        return;
    }

    size_t ix = ndc_axis_to_pixel(x, width);
    size_t iy = ndc_axis_to_pixel(y, height);

    screen[iy * width + ix] = value;
}

mat4x4 perspective(double fov_rad, double aspect, double near, double far) {
    mat4x4 result = empty();

    double f = 1.0 / tan(fov_rad / 2.0);

    result.values[0].x = f / aspect;
    result.values[1].y = f;
    result.values[2].z = -(far + near) / (far - near);
    result.values[2].w = -(2.0 * far * near) / (far - near);
    result.values[3].z = -1.0;
    result.values[3].w = 0.0;

    return result;
}

mat4x4 translate(double x, double y, double z) {
    mat4x4 result = identity();

    result.values[3].x = x;
    result.values[3].y = y;
    result.values[3].z = z;

    return result;
}

mat4x4 scale(double x, double y, double z) {
    mat4x4 result = identity();

    result.values[0].x = x;
    result.values[1].y = y;
    result.values[2].z = z;

    return result;
}

int main(void) {
    const size_t WIDTH = 80;
    const size_t HEIGHT = 40;

    char screen[WIDTH][HEIGHT];

    double rx = 0.0;
    double ry = 0.0;
    double rz = 0.0;

    const double TARGET_FPS = 60.0;
    const double TARGET_FRAME_TIME = 1.0 / TARGET_FPS;

    vec4 model[12] = {
        {-0.525731, 0.000000, 0.850651, 1.0},
        {0.525731, 0.000000, 0.850651, 1.0},
        {-0.525731, 0.000000, -0.850651, 1.0},
        {0.525731, 0.000000, -0.850651, 1.0},
        {0.000000, 0.850651, 0.525731, 1.0},
        {0.000000, 0.850651, -0.525731, 1.0},
        {0.000000, -0.850651, 0.525731, 1.0},
        {0.000000, -0.850651, -0.525731, 1.0},
        {0.850651, 0.525731, 0.000000, 1.0},
        {-0.850651, 0.525731, 0.000000, 1.0},
        {0.850651, -0.525731, 0.000000, 1.0},
        {-0.850651, -0.525731, 0.000000, 1.0},
    };

    uvec3 tris[20] = {
        {0, 4, 1},
        {0, 9, 4},
        {9, 5, 4},
        {4, 5, 8},
        {4, 8, 1},
        {8, 10, 1},
        {8, 3, 10},
        {5, 3, 8},
        {5, 2, 3},
        {2, 7, 3},
        {7, 10, 3},
        {7, 6, 10},
        {7, 11, 6},
        {11, 0, 6},
        {0, 1, 6},
        {6, 1, 10},
        {9, 0, 11},
        {9, 11, 2},
        {9, 2, 5},
        {7, 2, 11},
    };

    triangle triangles[20];

    clock_t last_time = clock();

    mat4x4 proj = perspective(PI * 0.5, (double)WIDTH / (double)HEIGHT, 0.1, 100.0);
    mat4x4 trans = translate(0.0, 0.0, -4.0);
    mat4x4 scl = scale(1.0, 1.0, 1.0);

    while (true) {
        clock_t now = clock();
        double delta_time = (double)(now - last_time) / CLOCKS_PER_SEC;

        if (delta_time >= TARGET_FRAME_TIME) {
            last_time = now;

            clear((char*)screen, WIDTH, HEIGHT);

            mat4x4 final = identity();

            mat4x4 rot_x = rotation_x(rx);
            mat4x4 rot_y = rotation_y(ry);
            mat4x4 rot_z = rotation_z(rz);

            final = multiply(&final, &proj);
            final = multiply(&final, &trans);
            final = multiply(&final, &rot_z);
            final = multiply(&final, &rot_y);
            final = multiply(&final, &rot_x);
            final = multiply(&final, &scl);

            vec2 points[8];

            for (size_t i = 0; i < 12; i++) {
                vec4 point = transform(&final, &model[i]);

                point.x /= point.w;
                point.y /= point.w;
                point.z /= point.w;
                point.y *= 0.65;

                if (point.x > 1.0 && point.x < -1.0) {
                    continue;
                }

                if (point.y > 1.0 && point.y < -1.0) {
                    continue;
                }

                if (point.z > 1.0 && point.z < -1.0) {
                    continue;
                }

                vec2 trunc = {
                    .x = point.x,
                    .y = point.y,
                };

                points[i] = trunc;
            }

            for (size_t i = 0; i < 20; i++) {
                uvec3* indices = &tris[i];

                triangle* tri = &triangles[i];

                tri->a = points[indices->x];
                tri->b = points[indices->y];
                tri->c = points[indices->z];
            }

            for (size_t y = 0; y < HEIGHT; y++) {
                for (size_t x = 0; x < WIDTH; x++) {
                    vec2 p = {
                        .x = pixel_axis_to_ndc(x, WIDTH),
                        .y = pixel_axis_to_ndc(y, HEIGHT),
                    };

                    for (size_t i = 0; i < 20; i++) {
                        if (triangle_test(&triangles[i], &p)) {
                            set(&p, i + 80, (char*)screen, WIDTH, HEIGHT);
                        }
                    }
                }
            }

            rx += 0.25 * PI * delta_time;
            ry += 0.25 * PI * delta_time;
            rz += 0.25 * PI * delta_time;

            present((char*)screen, WIDTH, HEIGHT);
        }
    }

    return 0;
}