#ifndef SPLINES_CARDINALSPLINES_HPP
#define SPLINES_CARDINALSPLINES_HPP


// https://stackoverflow.com/questions/7054272/how-to-draw-smooth-curve-through-n-points-using-javascript-html5-canvas

/*
 * @params:
 * - vertices: the points that need to be interpolated. The actual type can be anything as long as some basic operators are defined.
 * - num_vertices: the amount of points
 * - granularity: the desired amount of segments between consecutive points
 * - tension: magic value between 0 and 1. see https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Cardinal_spline
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

        cardinals[i] = {st2pow3 - st3pow2 + 1.0,    // coefficient of p0
                        st3pow2 - st2pow3,          // coefficient of p1
                        stpow3 - 2.0 * stpow2 + st, // coefficient of m0
                        stpow3 - stpow2,            // coefficient of m1
        };
    }

    const auto dealWithSection = [&](const point_type &vm1,
                                     const point_type &v0,
                                     const point_type &v1,
                                     const point_type &v2) {

        const auto m0 = (v1 - vm1) * (1.0 - tension) / 2.0; // starting tangent
        const auto m1 = (v2 - v0) * (1.0 - tension) / 2.0; // ending tangent

        for (std::size_t t = 0; t < numOfSeg; t++) {
            const auto &c = cardinals[t];
            cf(v0 * c.c1 + v1 * c.c2 + m0 * c.c3 + m1 * c.c4);
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

#endif //SPLINES_CARDINALSPLINES_HPP
