#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <algorithm>

#include "common.h"
#include "timer.h"
#include "vec.h"
#include "random.h"
#include "argparse.h"
#include "microfacet.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

//! Microfacet distribution
static std::unique_ptr<MicrofacetDistribution> microfacet = nullptr;

//! Ray class
struct Ray { Vec o, d; Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };

//! BSDF types
enum Refl_t {
    DIFFUSE,          // Diffuse
    CONDUCTOR,        // Smooth conductor
    DIELECTRIC,       // Smooth dielectric
    ROUGH_CONDUCTOR,  // Rough conductor
    ROUGH_DIELECTRIC  // Rough dielectric
};

//! Sphere class
struct Sphere {
    double rad;
    Vec pos, emit, color;
    Refl_t refl;

    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_)
        : rad(rad_), pos(p_), emit(e_), color(c_), refl(refl_) {}

    double intersect(const Ray &r) const {
        Vec op = pos - r.o;
        double t, eps = 1e-4, b = op.dot(r.d), det = b*b - op.dot(op) + rad*rad;
        if (det < 0) return 0; else det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

Sphere spheres[] = {  //Scene: radius, position, emission, color, material
  Sphere(600, Vec(50,681.6 - .27,81.6), Vec(12,12,12),  Vec(), DIFFUSE),          // Light
  Sphere(1e5, Vec(1e5 + 1,40.8,81.6),   Vec(),Vec(.95,.85,.15),DIFFUSE),          // Left
  Sphere(1e5, Vec(-1e5 + 99,40.8,81.6), Vec(),Vec(.15,.85,.95),DIFFUSE),          // Rght
  Sphere(1e5, Vec(50,40.8, 1e5),        Vec(),Vec(.75,.75,.75),DIFFUSE),          // Back
  Sphere(1e5, Vec(50,40.8,-1e5 + 170),  Vec(),Vec(),           DIFFUSE),          // Front
  Sphere(1e5, Vec(50, 1e5, 81.6),       Vec(),Vec(.75,.75,.75),DIFFUSE),          // Bottom
  Sphere(1e5, Vec(50,-1e5 + 81.6,81.6), Vec(),Vec(.75,.75,.75),DIFFUSE),          // Top
  Sphere(16.5,Vec(27,16.5,47),          Vec(),Vec(1,1,1)*.999, ROUGH_CONDUCTOR),  // Mirror
  Sphere(16.5,Vec(73,16.5,78),          Vec(),Vec(1,1,1)*.999, ROUGH_DIELECTRIC)  // Glass
};

const Sphere &light = spheres[0];

inline bool intersect(const Ray &r, double &t, int &id) {
    double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
    for (int i = int(n); i--;) if ((d = spheres[i].intersect(r)) && d < t) { t = d; id = i; }
    return t < inf;
}

//! Compute radiance
Vec radiance(const Ray &r) {
    Vec beta(1.0, 1.0, 1.0);
    Vec L(0.0);
    double pdf = 1.0;
    Ray ray = r;

    for (int depth = 0; ; depth++) {
        // Intersection test
        double t;
        int id = 0;
        if (!intersect(ray, t, id) || depth >= 16) {
            return Vec(0.0);
        }

        const Sphere &obj = spheres[id];
        Vec x = ray.o + ray.d * t;
        Vec n = (x - obj.pos).normalize();
        Vec nl = n.dot(ray.d) < 0 ? n : -n;
        Vec f = obj.color;

        L = L + beta * obj.emit / pdf;

        // Russian roulette
        double p = std::max(f.x, std::max(f.y, f.z));
        if (++depth > 5) {
            if (mt.next() < p) {
                f = f * (1 / p);
            } else {
                break;
            }
        }

        // Sample next direction
        const Vec wi = -ray.d;
        if (obj.refl == DIFFUSE) {
            // Diffuse reflection
            double r1 = 2 * Pi * mt.next();
            double r2 = mt.next();
            double r2s = sqrt(r2);
            Vec w = nl;
            Vec u = std::abs(w.x) > 0.1 ? Vec(0.0, 1.0, 0.0) : Vec(1.0, 0.0, 0.0).cross(w).normalize();
            Vec v = w.cross(u);
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).normalize();

            beta = f * beta;
            ray = Ray(x, d);
        } else if (obj.refl == CONDUCTOR) {
            // Ideal specular reflection
            beta = f * beta;
            ray = Ray(x, n * 2.0 * n.dot(wi) - wi);
        } else if (obj.refl == ROUGH_CONDUCTOR) {
            // Convert world space normal to local space normal
            Vec w = n;
            Vec u = std::abs(w.y) < 1.0-6 ? Vec(1.0, 0.0, 0.0) : Vec(0.0, 1.0, 0.0).cross(w).normalize();
            Vec v = w.cross(u);
            Vec wiLocal = Vec(u.dot(wi), v.dot(wi), w.dot(wi)).normalize();

            // Rough specular reflection
            Vec wmLocal = microfacet->sampleWm(wiLocal);
            Vec wm = u * wmLocal.x + v * wmLocal.y + w * wmLocal.z;
            Vec woLocal = wmLocal * 2.0 * wmLocal.dot(wiLocal) - wiLocal;
            Vec wo = u * woLocal.x + v * woLocal.y + w * woLocal.z;

            // Check validity of reflection
            if (wiLocal.z * woLocal.z < 0.0) {
                break;
            }

            beta = f * beta * microfacet->weight(woLocal, wiLocal, wmLocal);
            ray = Ray(x, wo);
        } else if (obj.refl == DIELECTRIC) {
            // Ideal dielectirc
            Ray reflRay(x, n * 2.0 * n.dot(wi) - wi);
            bool into = n.dot(nl) > 0;

            // Snell's rule
            double nc = 1.0;
            double nt = 1.5;
            double nnt = into ? nc / nt : nt / nc;
            double ddn = -std::abs(n.dot(wi));
            double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
            if (cos2t < 0.0) {
                // Total reflection
                beta = f * beta;
                ray = reflRay;
                continue;
            }

            // Refraction (with Schlick approximation)
            Vec tdir = (-wi * nnt - n * ((into ? 1 : -1) * (ddn*nnt + sqrt(cos2t)))).normalize();
            double a = nt - nc;
            double b = nt + nc;
            double R0 = (a * a) / (b * b);
            double c = 1.0 - (into ? -ddn : tdir.dot(n));
            double Re = R0 + (1 - R0) * std::pow(c, 5.0);
            double Tr = 1 - Re;
            double P = 0.25 + 0.5 * Re;
            double RP = Re / P;
            double TP = Tr / (1.0 - P);
            if (mt.next() < P) {
                beta = f * beta;
                ray = reflRay;
                pdf /= RP;
            } else {
                beta = f * beta;
                ray = Ray(x, tdir);
                pdf /= TP;
            }
        } else if (obj.refl == ROUGH_DIELECTRIC) {
            // Convert world space normal to local space normal
            Vec w = n;
            Vec u = std::abs(w.y) < 1.0-6 ? Vec(1.0, 0.0, 0.0) : Vec(0.0, 1.0, 0.0).cross(w).normalize();
            Vec v = w.cross(u);
            Vec wiLocal = Vec(u.dot(wi), v.dot(wi), w.dot(wi)).normalize();

            // Rough dielectric
            Vec wmLocal = microfacet->sampleWm(wiLocal);
            Vec wm = u * wmLocal.x + v * wmLocal.y + w * wmLocal.z;

            // Reflection
            Vec woLocalRe = wmLocal * 2.0 * wmLocal.dot(wiLocal) - wiLocal;
            bool into = n.dot(nl) > 0.0;

            // Refraction
            double nc = 1.0;
            double nt = 1.5;
            double nnt = into ? nc / nt : nt / nc;
            double ddn = -std::abs(wi.dot(wm));
            double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

            Vec woLocalTr = (-wiLocal * nnt - wmLocal * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalize();
            double a = nt - nc;
            double b = nt + nc;
            double R0 = (a * a) / (b * b);
            double c = 1.0 - (into ? -ddn : woLocalTr.dot(wmLocal));
            double Re = R0 + (1 - R0) * std::pow(c, 5.0);
            double Tr = 1 - Re;
            double P = 0.25 + 0.5 * Re;
            double RP = Re / P;
            double TP = Tr / (1.0 - P);

            Vec woLocal;
            bool isRefl = false;
            if (cos2t < 0.0) {
                // Total reflection
                woLocal = woLocalRe;
                isRefl = true;
            } else {
                if (mt.next() < P) {
                    // Reflection
                    woLocal = woLocalRe;
                    pdf /= RP;
                    isRefl = true;
                } else {
                    // Transmission
                    woLocal = woLocalTr;
                    pdf /= TP;
                    isRefl = false;
                }
            }

            // Check validity of reflection / transmission
            if ((wiLocal.z * woLocal.z < 0.0 && isRefl) ||
                (wiLocal.z * woLocal.z > 0.0 && !isRefl)) {
                break;
            }

            Vec wo = u * woLocal.x + v * woLocal.y + w * woLocal.z;
            beta = f * beta * microfacet->weight(woLocal, wiLocal, wmLocal);
            ray = Ray(x, wo);
        }
    }

    return L;
}

int main(int argc, char **argv) {
    auto &parser = ArgumentParser::getInstance();
    try {
        parser.addArgument("-w", "--width", 960, false, "Image width");
        parser.addArgument("-h", "--height", 720, false, "Image height");
        parser.addArgument("-s", "--samples", 32, false, "Samples per pixel");
        parser.addArgument("", "--distrib", "beckmann", false,
                           "Name of distribution (ggx or beckmann)");
        parser.addArgument("", "--samplevis", "yes", false, "Sample visible normals or not");
        parser.addArgument("-o", "--output", "image.png", false, "Output file name");
        parser.addArgument("", "--alphax", 0.1, false, "Roughness for x-direction");
        parser.addArgument("", "--alphay", 0.1, false, "Roughness for y-direction");
        parser.parse(argc, argv);
    } catch (std::runtime_error &e) {
        std::cerr << parser.helpText() << std::endl;
        std::cerr << e.what() << std::endl;
        return -1;
    }

    std::cerr << parser.helpText() << std::endl;
    parser.print();

    const int w = parser.getInt("width");
    const int h = parser.getInt("height");
    const int samples = parser.getInt("samples");
    const std::string distribution = parser.getString("distrib");
    const bool sampleVis = parser.getBool("samplevis");
    const std::string outfile = parser.getString("output");
    const double alphax = parser.getDouble("alphax");
    const double alphay = parser.getDouble("alphay");

    if (distribution == "beckmann") {
        microfacet.reset(new BeckmannDistribution(alphax, alphay, sampleVis));
    } else if (distribution == "ggx") {
        microfacet.reset(new GGXDistribution(alphax, alphay, sampleVis));
    } else {
        fprintf(stderr, "Unknown distribution type: %s\n", distribution.c_str());
        return -1;
    }

    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize());
    Vec cx = Vec(w * 0.5135 / h, 0.0, 0.0);
    Vec cy = cx.cross(cam.d).normalize() * 0.5135;
    Vec r;
    auto c = std::unique_ptr<Vec[]>(new Vec[w * h]);

    Timer timer;
    timer.start();

    #pragma omp parallel for schedule(dynamic, 1) private(r) num_threads(7)
    for (int y = 0; y < h; y++) {
        fprintf(stderr, "\rRendering (%d spp) %6.2f %%", samples * 4, 100.*y / (h - 1));
        for (int x = 0; x < w; x++) {
            for (int sy = 0, i = (h - y - 1)*w + x; sy < 2; sy++) {
                for (int sx = 0; sx < 2; sx++, r = Vec(0.0)) {
                    for (int s = 0; s < samples; s++) {
                        double r1 = 2 * mt.next(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * mt.next(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx*(((sx + .5 + dx) / 2 + x) / w - .5) +
                            cy*(((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        r = r + radiance(Ray(cam.o + d * 140, d.normalize())) * (1. / samples);
                    }
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
                }
            }
        }
    }

    printf("Time: %f sec\n", timer.stop());

    // Save image with std_image_write
    auto bytes = std::unique_ptr<uint8_t[]>(new uint8_t[w * h * 3]);
    for (int i = 0; i < w * h; i++) {
        bytes[i * 3 + 0] = toInt(c[i].x);
        bytes[i * 3 + 1] = toInt(c[i].y);
        bytes[i * 3 + 2] = toInt(c[i].z);
    }
    stbi_write_png(outfile.c_str(), w, h, 3, bytes.get(), w * 3);
    fprintf(stdout, "\nFinish!!\n");
}
