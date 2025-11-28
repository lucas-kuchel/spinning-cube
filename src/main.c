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

typedef struct mat4x4 {
    vec4 values[4];
} mat4x4;

double dot(const vec4* a, const vec4* b) {
    return (a->x * b->x) + (a->y * b->y) + (a->z * b->z) + (a->w * b->w);
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

            result.values[i].values[j] = dot(row, col);
        }
    }

    return result;
}

vec4 transform(const mat4x4* a, const vec4* b) {
    mat4x4 ax = transpose(a);
    vec4 result;

    for (size_t i = 0; i < 4; i++) {
        const vec4* row = &ax.values[i];

        result.values[i] = dot(row, b);
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

    // printf("\x1b[H");
}

void present(char* screen, size_t width, size_t height) {
    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            putchar(screen[y * width + x]);
        }

        putchar('\n');
    }
}

size_t ndc_axis_to_pixel(float axis, size_t min, size_t max) {
    if (axis < -1.0) {
        axis = -1.0;
    }

    if (axis > 1.0) {
        axis = 1.0;
    }

    double normalized = (axis + 1.0) * 0.5;
    double pixel = min + normalized * (max - min);

    return (size_t)round(pixel);
}

void set(vec4* point, char value, char* screen, size_t width, size_t height) {
    double x = point->x;
    double y = point->y;

    if (x > 1.0 || x < -1.0 || y > 1.0 || y < -1.0) {
        return;
    }

    size_t ix = ndc_axis_to_pixel(x, 0, width);
    size_t iy = ndc_axis_to_pixel(y, 0, height);

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

int main(void) {
    const size_t WIDTH = 80;
    const size_t HEIGHT = 40;

    char screen[WIDTH][HEIGHT];

    double rx = 0.0;
    double ry = 0.0;
    double rz = 0.0;

    const double TARGET_FPS = 60.0;
    const double TARGET_FRAME_TIME = 1.0 / TARGET_FPS;

    vec4 model[8] = {
        {.x = 0.5, .y = 0.5, .z = 0.5, .w = 1.0},
        {.x = -0.5, .y = 0.5, .z = 0.5, .w = 1.0},
        {.x = 0.5, .y = -0.5, .z = 0.5, .w = 1.0},
        {.x = -0.5, .y = -0.5, .z = 0.5, .w = 1.0},
        {.x = 0.5, .y = 0.5, .z = -0.5, .w = 1.0},
        {.x = -0.5, .y = 0.5, .z = -0.5, .w = 1.0},
        {.x = 0.5, .y = -0.5, .z = -0.5, .w = 1.0},
        {.x = -0.5, .y = -0.5, .z = -0.5, .w = 1.0},
    };

    clock_t last_time = clock();

    mat4x4 proj = perspective(PI * 0.5, (double)WIDTH / (double)HEIGHT, 0.1, 100.0);

    while (true) {
        clock_t now = clock();
        double delta_time = (double)(now - last_time) / CLOCKS_PER_SEC;
        mat4x4 trans = translate(0.0, 0.0, -3.0);

        if (delta_time >= TARGET_FRAME_TIME) {
            last_time = now;

            clear((char*)screen, WIDTH, HEIGHT);

            mat4x4 rot_x = rotation_x(rx);
            mat4x4 rot_y = rotation_y(ry);
            mat4x4 rot_z = rotation_z(rz);

            mat4x4 finalTransform = multiply(&proj, &trans);

            finalTransform = multiply(&finalTransform, &rot_z);
            finalTransform = multiply(&finalTransform, &rot_y);
            finalTransform = multiply(&finalTransform, &rot_x);

            for (size_t i = 0; i < 8; i++) {
                vec4 point = transform(&finalTransform, &model[i]);

                point.x /= point.w;
                point.y /= point.w;
                point.z /= point.w;
                point.y *= 0.5;

                set(&point, '@', (char*)screen, WIDTH, HEIGHT);
            }

            rx += 0.25 * PI * delta_time;
            ry += 0.25 * PI * delta_time;
            rz += 0.25 * PI * delta_time;

            present((char*)screen, WIDTH, HEIGHT);
        }
    }

    return 0;
}