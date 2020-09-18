#include <iostream>
#include <cmath>
#include <cstring>
#include <array>

template<class T, std::size_t n = 3>
struct Point {
    using value_type = T;
    std::array<T, n> mat;

    constexpr std::size_t size() const {
        return n;
    }


    Point(const Point<T, n> &cp) : mat(cp.mat) {
    }


    Point(Point<T, n> &&cp) : mat(std::move(cp.mat)) {
    }

    // default constructor
    Point() {
        mat.fill(0);
    }

    // a variadic constructor with sfinae ...
    template<typename... E, typename std::enable_if<sizeof...(E) == n, bool>::type = true>
    Point(E &&...e) :mat{{std::forward<E>(e)...}} {
    }

    Point<T, n> &operator=(const Point<T, n> &o) {
        mat = o.mat;
        return *this;
    }


    const T &operator[](std::size_t i) const {
//        assert(i < n);
        return mat[i];
    }

    T &operator[](std::size_t i) {
//        assert(i < n);
        return mat[i];
    }


    //sum of two Points
    Point<T, n> &operator+=(const Point<T, n> &pt) {
        for (std::size_t i = 0; i < n; i++)
            mat[i] += pt[i];

        return (*this);
    }

    //Point + some value
    Point<T, n> &operator+=(const T &val) {
        for (std::size_t i = 0; i < n; i++)
            mat[i] += val;

        return (*this);
    }


    Point<T, n> operator+(const Point<T, n> &val) const {
        Point<T, n> ret;
        for (std::size_t i = 0; i < n; i++)
            ret[i] = mat[i] + val[i];
        return ret;
    }

    Point<T, n> operator+(const T &val) const {
        Point<T, n> ret;
        for (std::size_t i = 0; i < n; i++)
            ret[i] = mat[i] + val;
        return ret;
    }


    Point<T, n> operator-(const Point<T, n> &val) const {
        Point<T, n> ret;
        for (std::size_t i = 0; i < n; i++)
            ret[i] = mat[i] - val[i];
        return ret;
    }

    Point<T, n> operator-(const T &val) const {
        Point<T, n> ret;
        for (std::size_t i = 0; i < n; i++)
            ret[i] = mat[i] - val;
        return ret;
    }

    Point<T, n> operator*(const Point<T, n> &val) const {
        Point<T, n> ret;
        for (std::size_t i = 0; i < n; i++)
            ret[i] = mat[i] * val[i];
        return ret;
    }

    T dot(const Point<T, n> &val) const { // dot-product
        T ret = 0;
        for (std::size_t i = 0; i < n; i++)
            ret += mat[i] * val[i];
        return ret;
    }

    Point<T, n> operator*(const T &val) const { // dot-product
        Point<T, n> ret;
        for (std::size_t i = 0; i < n; i++)
            ret[i] = mat[i] * val;
        return ret;
    }

    Point<T, n> operator/(const T &val) const { // dot-product
        Point<T, n> ret;
        for (std::size_t i = 0; i < n; i++)
            ret[i] = mat[i] / val;
        return ret;
    }

    Point<T, n> &operator-=(const Point<T, n> &val) {
        for (std::size_t i = 0; i < n; i++)
            mat[i] -= val[i];
        return *this;
    }

    Point<T, n> &operator/=(const Point<T, n> &val) {
        for (std::size_t i = 0; i < n; i++)
            mat[i] /= val[i];
        return *this;
    }

    Point<T, n> &operator/=(const T val) {
        for (std::size_t i = 0; i < n; i++)
            mat[i] /= val;
        return *this;
    }

    T norm() const {
        T ret = 0;
        for (std::size_t i = 0; i < n; i++)
            ret += mat[i] * mat[i];
        return sqrt(ret);
    }

/*    void sanitycheck(double lower, double upper) {
        for (int i = 0; i < n; i++) {
            assert(mat[i] < upper && mat[i] > lower);
        }
    }*/
};


// https://stackoverflow.com/questions/7054272/how-to-draw-smooth-curve-through-n-points-using-javascript-html5-canvas

/*
 * @params:
 * - vertices: the points that need to be interpolated. The actual type can be anything as long as some basic operators are defined.
 * - num_vertices: the amount of points
 * - granularity: the desired amount of segments between consecutive points
 * - tension: magic value between 0 and 1. see https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Cardinal_spline
 *         0.5 gives reasonable results.
 * - close: wether the sequence of vertices should be considered a closed loop or not.
 * - cf: callbackfunction, it will be called for each generated vertex.
 */

template<typename point_type, typename callbackfct>
void cardinalpoints(point_type const *const vertices,
                    const std::size_t num_vertices,
                    const std::size_t granularity,
                    const double tension,
                    const bool close,
                    const callbackfct &cf) {

    const std::size_t numOfSeg = granularity;

    struct Cardinals {
        double c1, c2, c3, c4;
    };
    Cardinals *cardinals = (Cardinals *) alloca(sizeof(Cardinals) * numOfSeg);
    memset(cardinals, sizeof(Cardinals) * numOfSeg, 0);

    for (std::size_t i = 0; i < numOfSeg; i++) {
        //c1 =   2 * Math.pow(st, 3)  - 3 * Math.pow(st, 2) + 1;
        //c2 = -(2 * Math.pow(st, 3)) + 3 * Math.pow(st, 2);
        //c3 =       Math.pow(st, 3)  - 2 * Math.pow(st, 2) + st;
        //c4 =       Math.pow(st, 3)  -     Math.pow(st, 2);

        const double st = double(i) / double(numOfSeg - 1);
        const double stpow2 = st * st;
        const double stpow3 = stpow2 * st;
        const double st2pow3 = stpow3 * 2.0;
        const double st3pow2 = stpow2 * 3.0;

        cardinals[i] = {st2pow3 - st3pow2 + 1.0,
                        st3pow2 - st2pow3,
                        stpow3 - 2.0 * stpow2 + st,
                        stpow3 - stpow2,
        };
    }

    const auto dealWithSection = [&](const point_type &vm1,
                                     const point_type &v0,
                                     const point_type &v1,
                                     const point_type &v2) {

        const auto t1 = (v1 - vm1) * tension;
        const auto t2 = (v2 - v0) * tension;
//        cf(v0);
        for (std::size_t t = 0; t < numOfSeg; t++) {
            const auto &c = cardinals[t];
            cf(v0 * c.c1 + v1 * c.c2 + t1 * c.c3 + t2 * c.c4);
        }
    };

    if (close) {
        dealWithSection(vertices[num_vertices - 1],
                        vertices[0],
                        vertices[1],
                        vertices[2]);
    } else {
        dealWithSection(vertices[0],
                        vertices[0],
                        vertices[1],
                        vertices[2]);
    }

    for (std::size_t i = 1; i + 2 < num_vertices; ++i) {
        dealWithSection(vertices[i - 1],
                        vertices[i],
                        vertices[i + 1],
                        vertices[i + 2]);
    }

    if (close) {
        dealWithSection(vertices[num_vertices - 3],
                        vertices[num_vertices - 2],
                        vertices[num_vertices - 1],
                        vertices[0]);
        dealWithSection(vertices[num_vertices - 2],
                        vertices[num_vertices - 1],
                        vertices[0],
                        vertices[1]);
    } else {
        dealWithSection(vertices[num_vertices - 3],
                        vertices[num_vertices - 2],
                        vertices[num_vertices - 1],
                        vertices[num_vertices - 1]);
    }

}

template<typename points_container>
void printCardinalSpline(std::size_t granularity, double tension, double close, const points_container &input) {


//    for (const auto &i : input) {
//        std::cerr << "1 " << i[0] << ' ' << i[1] << std::endl;
//    }
//    std::cerr << std::endl;

    cardinalpoints(input.data(),
                   input.size(),
                   granularity,
                   tension,
                   close,
                   [tension](const Point<double, 2> &p) {
                       std::cerr << "0 " << p[0] << " " << p[1] << '\n';
                   });
    std::cerr << std::endl;
}

int main(int argc, char *argv[]) {
    const std::size_t granularity = argc > 1 ? atoi(argv[1]) : 20;

    //C
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{1., 2.},
            Point<double, 2>{0., 2.},
            Point<double, 2>{0., 0.},
            Point<double, 2>{1., 0.},
    });


    //A
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 6>{
            Point<double, 2>{1.5, 0.},
            Point<double, 2>{2., 2.},
            Point<double, 2>{2.25, 1.},
            Point<double, 2>{1.75, 1.},
            Point<double, 2>{2.25, 1.},
            Point<double, 2>{2.5, 0.},
    });
    //R
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 5>{
            Point<double, 2>{3., 0.},
            Point<double, 2>{3., 2.},
            Point<double, 2>{4., 1.5},
            Point<double, 2>{3., 1.},
            Point<double, 2>{4., 0.},
    });
    //D
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{4.5, 0.},
            Point<double, 2>{4.5, 2.},
            Point<double, 2>{5.5, 1.},
            Point<double, 2>{4.5, 0.},
    });

    // i dot
    printCardinalSpline(granularity, 0.5, true, std::array<Point<double, 2>, 4>{
            Point<double, 2>{6., 1.75},
            Point<double, 2>{6. - 0.125, 1.75 + 0.125},
            Point<double, 2>{6., 2.},
            Point<double, 2>{6. + 0.125, 1.75 + 0.125},
    });

    // i base
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 3>{
            Point<double, 2>{6., 0.},
            Point<double, 2>{6., 1.5},
            Point<double, 2>{6., 0.},
    });

    // n
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{6.5, 0.},
            Point<double, 2>{6.5, 2.},
            Point<double, 2>{7.5, 0.},
            Point<double, 2>{7.5, 2.},
    });

    //A
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 6>{
            Point<double, 2>{8., 0.},
            Point<double, 2>{8.5, 2.},
            Point<double, 2>{8.75, 1.},
            Point<double, 2>{8.25, 1.},
            Point<double, 2>{8.75, 1.},
            Point<double, 2>{9., 0.},
    });

    //L
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 3>{
            Point<double, 2>{9.5, 2.},
            Point<double, 2>{9.5, 0.},
            Point<double, 2>{10.5, 0.},
    });


    //S
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{1., -1.},
            Point<double, 2>{0., -1.5},
            Point<double, 2>{1., -2.5},
            Point<double, 2>{0., -3.},
    });

    //P
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{1.5, -2.},
            Point<double, 2>{2.5, -1.5},
            Point<double, 2>{1.5, -1.},
            Point<double, 2>{1.5, -3.},
    });


    //L
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 3>{
            Point<double, 2>{3., -1.},
            Point<double, 2>{3., -3.},
            Point<double, 2>{4., -3.},
    });

// i dot
    printCardinalSpline(granularity, 0.5, true, std::array<Point<double, 2>, 4>{
            Point<double, 2>{4.5, 1.75 - 3.0},
            Point<double, 2>{4.5 - 0.125, 1.75 + 0.125 - 3.0},
            Point<double, 2>{4.5, 2 - 3.0},
            Point<double, 2>{4.5 + 0.125, 1.75 + 0.125 - 3.0},
    });

    // i base
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 3>{
            Point<double, 2>{4.5, 0 - 3.0},
            Point<double, 2>{4.5, 1.5 - 3.0},
            Point<double, 2>{4.5, 0 - 3.0},
    });

    // n
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{5., -3.0},
            Point<double, 2>{5., -1.0},
            Point<double, 2>{6., -3.0},
            Point<double, 2>{6., -1.0},
    });

    // E
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 7>{
            Point<double, 2>{7.5, 0 - 1.0},
            Point<double, 2>{6.5, 0 - 1.0},
            Point<double, 2>{6.5, 0 - 2.0},
            Point<double, 2>{7.5, 0 - 2.0},
            Point<double, 2>{6.5, 0 - 2.0},
            Point<double, 2>{6.5, 0 - 3.0},
            Point<double, 2>{7.5, 0 - 3.0},
    });

    //S
    printCardinalSpline(granularity, 0.5, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{1 + 8., -1.},
            Point<double, 2>{0 + 8., -1.5},
            Point<double, 2>{1 + 8., -2.5},
            Point<double, 2>{0 + 8., -3.},
    });


    return 0;
}