#include <iostream>
#include <cmath>
#include <cstring>
#include <array>

#include "cardinalsplines.hpp"

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

    double tension = 0;

    //C
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{1., 2.},
            Point<double, 2>{0., 2.},
            Point<double, 2>{0., 0.},
            Point<double, 2>{1., 0.},
    });


    //A
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 6>{
            Point<double, 2>{1.5, 0.},
            Point<double, 2>{2., 2.},
            Point<double, 2>{2.25, 1.},
            Point<double, 2>{1.75, 1.},
            Point<double, 2>{2.25, 1.},
            Point<double, 2>{2.5, 0.},
    });
    //R
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 5>{
            Point<double, 2>{3., 0.},
            Point<double, 2>{3., 2.},
            Point<double, 2>{4., 1.5},
            Point<double, 2>{3., 1.},
            Point<double, 2>{4., 0.},
    });
    //D
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{4.5, 0.},
            Point<double, 2>{4.5, 2.},
            Point<double, 2>{5.5, 1.},
            Point<double, 2>{4.5, 0.},
    });

    // i dot
    printCardinalSpline(granularity, tension, true, std::array<Point<double, 2>, 4>{
            Point<double, 2>{6., 1.75},
            Point<double, 2>{6. - 0.125, 1.75 + 0.125},
            Point<double, 2>{6., 2.},
            Point<double, 2>{6. + 0.125, 1.75 + 0.125},
    });

    // i base
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 3>{
            Point<double, 2>{6., 0.},
            Point<double, 2>{6., 1.5},
            Point<double, 2>{6., 0.},
    });

    // n
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{6.5, 0.},
            Point<double, 2>{6.5, 2.},
            Point<double, 2>{7.5, 0.},
            Point<double, 2>{7.5, 2.},
    });

    //A
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 6>{
            Point<double, 2>{8., 0.},
            Point<double, 2>{8.5, 2.},
            Point<double, 2>{8.75, 1.},
            Point<double, 2>{8.25, 1.},
            Point<double, 2>{8.75, 1.},
            Point<double, 2>{9., 0.},
    });

    //L
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 3>{
            Point<double, 2>{9.5, 2.},
            Point<double, 2>{9.5, 0.},
            Point<double, 2>{10.5, 0.},
    });


    //S
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{1., -1.},
            Point<double, 2>{0., -1.5},
            Point<double, 2>{1., -2.5},
            Point<double, 2>{0., -3.},
    });

    //P
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{1.5, -2.},
            Point<double, 2>{2.5, -1.5},
            Point<double, 2>{1.5, -1.},
            Point<double, 2>{1.5, -3.},
    });


    //L
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 3>{
            Point<double, 2>{3., -1.},
            Point<double, 2>{3., -3.},
            Point<double, 2>{4., -3.},
    });

// i dot
    printCardinalSpline(granularity, tension, true, std::array<Point<double, 2>, 4>{
            Point<double, 2>{4.5, 1.75 - 3.0},
            Point<double, 2>{4.5 - 0.125, 1.75 + 0.125 - 3.0},
            Point<double, 2>{4.5, 2 - 3.0},
            Point<double, 2>{4.5 + 0.125, 1.75 + 0.125 - 3.0},
    });

    // i base
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 3>{
            Point<double, 2>{4.5, 0 - 3.0},
            Point<double, 2>{4.5, 1.5 - 3.0},
            Point<double, 2>{4.5, 0 - 3.0},
    });

    // n
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{5., -3.0},
            Point<double, 2>{5., -1.0},
            Point<double, 2>{6., -3.0},
            Point<double, 2>{6., -1.0},
    });

    // E
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 7>{
            Point<double, 2>{7.5, 0 - 1.0},
            Point<double, 2>{6.5, 0 - 1.0},
            Point<double, 2>{6.5, 0 - 2.0},
            Point<double, 2>{7.5, 0 - 2.0},
            Point<double, 2>{6.5, 0 - 2.0},
            Point<double, 2>{6.5, 0 - 3.0},
            Point<double, 2>{7.5, 0 - 3.0},
    });

    //S
    printCardinalSpline(granularity, tension, false, std::array<Point<double, 2>, 4>{
            Point<double, 2>{1 + 8., -1.},
            Point<double, 2>{0 + 8., -1.5},
            Point<double, 2>{1 + 8., -2.5},
            Point<double, 2>{0 + 8., -3.},
    });


    return 0;
}