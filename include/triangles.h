#include <iostream>
#include <cmath>
template <typename T>
struct Vector3D {
    T x;
    T y;
    T z;

    Vector3D operator+(const Vector3D<T>& r) const { return {x+r.x, y+r.y, z+r.z}; } //перегрузка операторов
    Vector3D operator-(const Vector3D<T>& r) const { return {x-r.x, y-r.y, z-r.z}; }
    Vector3D operator*(float s) const          { return {x*s,   y*s,   z*s  }; }
    Vector3D operator/(float s) const          { return {x/s,   y/s,   z/s  }; }
};

template <typename T>
struct Point3D {
    T x;
    T y;
    T z;

     Point3D operator+(const Vector3D<T>& v) const {
        return {x + v.x, y + v.y, z + v.z};
    }
    Point3D operator-(const Vector3D<T>& v) const {
        return {x - v.x, y - v.y, z - v.z};
    }

    void display() const {
        std::cout << "Value: (" << x << " " << y << " " << z << ")" << std::endl;
    }
};

template <typename T>
struct Line3D {
    Point3D<T> origin;
    Vector3D<T> direction;
};

template <typename T>
struct Intersection {
    bool hit;      // true если есть пересечение
    T t;
    T u, v;
    Point3D<T> intersection_point;
};


template <typename T>
struct Triangle3D {
        Point3D<T> v0;
        Point3D<T> v1;
        Point3D<T> v2;
};


template <typename T>
Vector3D<T> crossProduct(const Vector3D<T>& a, const Vector3D<T>& b);

template <typename T>
T scalarProduct(const Vector3D<T>& a, const Vector3D<T>& b);

template <typename T>
Intersection<T> MollerTrumbore(const Triangle3D<T>& tri,
                               const Line3D<T>& ray,
                               T epsilon = 1e-6);

template <typename T>
Intersection<T> TriangleIntersection(const Triangle3D<T>& tri_1,
                                    const Triangle3D<T>& tri_2,
                                    T epsilon = 1e-6);

template <typename T>
Vector3D<T> normal(const Triangle3D<T>& tri);


template <typename T>
bool pointOnPlane(const Triangle3D<T>& tri, const Point3D<T>& p, T eps);


template <typename T>
bool pointInTriangle(const Triangle3D<T>& tri, const Point3D<T>& p, T eps);




template <typename T>
Vector3D<T> crossProduct(const Vector3D<T>& a, const Vector3D<T>& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

template <typename T>
T scalarProduct(const Vector3D<T>& a, const Vector3D<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;

}

template <typename T>
Vector3D<T> normal(const Triangle3D<T>& tri) {
    Vector3D<T> e1{tri.v1.x - tri.v0.x, tri.v1.y - tri.v0.y, tri.v1.z - tri.v0.z};
    Vector3D<T> e2{tri.v2.x - tri.v0.x, tri.v2.y - tri.v0.y, tri.v2.z - tri.v0.z};
    return crossProduct(e1, e2);
}

template <typename T>
bool pointOnPlane(const Triangle3D<T>& tri, const Point3D<T>& p, T eps) {
    Vector3D<T> n = normal(tri);
    Vector3D<T> v{p.x - tri.v0.x, p.y - tri.v0.y, p.z - tri.v0.z};
    T d = scalarProduct(n, v);
    return std::fabs(d) <= eps * std::sqrt(scalarProduct(n,n));
}

template <typename T>
bool pointInTriangle(const Triangle3D<T>& tri, const Point3D<T>& p, T eps) {

    Vector3D<T> v0{tri.v2.x - tri.v0.x, tri.v2.y - tri.v0.y, tri.v2.z - tri.v0.z};
    Vector3D<T> v1{tri.v1.x - tri.v0.x, tri.v1.y - tri.v0.y, tri.v1.z - tri.v0.z};
    Vector3D<T> v2{p.x - tri.v0.x, p.y - tri.v0.y, p.z - tri.v0.z};

    T dot00 = scalarProduct(v0,v0);
    T dot01 = scalarProduct(v0,v1);
    T dot02 = scalarProduct(v0,v2);
    T dot11 = scalarProduct(v1,v1);
    T dot12 = scalarProduct(v1,v2);

    T denom = dot00*dot11 - dot01*dot01;
    if (std::fabs(denom) < eps) return false;

    T u = (dot11*dot02 - dot01*dot12) / denom;
    T v = (dot00*dot12 - dot01*dot02) / denom;
    return u >= -eps && v >= -eps && (u+v) <= 1.0 + eps;
}



template <typename T>
Intersection<T> MollerTrumbore(const Triangle3D<T>& tri,
                               const Line3D<T>& ray,
                               T epsilon)
{
    Vector3D<T> e1 = {tri.v1.x - tri.v0.x, tri.v1.y - tri.v0.y, tri.v1.z - tri.v0.z};
    Vector3D<T> e2 = {tri.v2.x - tri.v0.x, tri.v2.y - tri.v0.y, tri.v2.z - tri.v0.z};

    Vector3D<T> p = crossProduct(ray.direction, e2);
    T det = scalarProduct(e1, p);

    if (det > -epsilon && det < epsilon)
        return {false, 0, 0, 0, {0,0,0}}; // параллельность

    T invDet = 1 / det;
    Vector3D<T> s = {ray.origin.x - tri.v0.x, ray.origin.y - tri.v0.y, ray.origin.z - tri.v0.z};
    T u = scalarProduct(s, p) * invDet;
    if (u < 0 || u > 1) return {false, 0, 0, 0, {0,0,0}};

    Vector3D<T> q = crossProduct(s, e1);
    T v = scalarProduct(ray.direction, q) * invDet;
    if (v < 0 || u + v > 1) return {false, 0, 0, 0, {0,0,0}};

    T t = scalarProduct(e2, q) * invDet;
    Point3D<T> intersectPoint = ray.origin + ray.direction * t;

    return {true, t, u, v, intersectPoint};
}


template <typename T>
Intersection<T> TriangleIntersection(const Triangle3D<T>& tri_1,
                                     const Triangle3D<T>& tri_2,
                                     T epsilon)
{

    for (int i = 0; i < 3; ++i) { // // Проверяем все рёбра tri_2
        const Point3D<T>& p0 = (&tri_2.v0)[i];
        const Point3D<T>& p1 = (&tri_2.v0)[(i+1)%3];
        Line3D<T> edge;
        edge.origin = p0;
        edge.direction = {p1.x - p0.x, p1.y - p0.y, p1.z - p0.z};
        auto hit = MollerTrumbore(tri_1, edge, epsilon);
        if (hit.hit && hit.t >= -epsilon && hit.t <= 1.0 + epsilon) return hit;
    }


    if (pointOnPlane(tri_1, tri_2.v0, epsilon) &&  //проверка того, лежит ли tri_2 целиком в плоскости tri_1 и хотя бы одна вершина внутри
        pointOnPlane(tri_1, tri_2.v1, epsilon) &&
        pointOnPlane(tri_1, tri_2.v2, epsilon))
    {
        for (int i = 0; i < 3; ++i) {
            const Point3D<T>& p = (&tri_2.v0)[i];
            if (pointInTriangle(tri_1, p, epsilon)) {
                return {true, 0, 0, 0, p}; // точка пересечения – любая вершина
            }
        }
    }

    if (pointOnPlane(tri_2, tri_1.v0, epsilon) &&  //проверка того, лежит ли tri_1 внутри tri_2
        pointOnPlane(tri_2, tri_1.v1, epsilon) &&
        pointOnPlane(tri_2, tri_1.v2, epsilon))
    {
        for (int i = 0; i < 3; ++i) {
            const Point3D<T>& p = (&tri_1.v0)[i];
            if (pointInTriangle(tri_2, p, epsilon)) {
                return {true, 0, 0, 0, p};
            }
        }
    }
    return {false, 0, 0, 0, {0,0,0}};
}





