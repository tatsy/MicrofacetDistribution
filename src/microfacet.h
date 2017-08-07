#ifndef _MICROFACET_H_
#define _MICROFACET_H_

#include "common.h"
#include "vec.h"
#include "random.h"

namespace {

// ----------------------------------------------------------------------------
// Utility functions to compute trigonometric funtions for a vector
// ----------------------------------------------------------------------------

static inline double cosTheta(const Vec &w) {
    return w.z;
}

static inline double cos2Theta(const Vec &w) {
    double c = cosTheta(w);
    return c * c;
}

static inline double sinTheta(const Vec &w) {
    return std::sqrt(std::max(0.0, 1.0 - cos2Theta(w)));
}

static inline double tanTheta(const Vec &w) {
    return sinTheta(w) / cosTheta(w);
}

static inline double cosPhi(const Vec &w) {
    if (w.z == 1.0) return 0.0;

    return w.x / sinTheta(w);
}

static inline double cos2Phi(const Vec &w) {
    double c = cosPhi(w);
    return c * c;
}

static inline double sinPhi(const Vec &w) {
    if (w.z == 1.0) return 0.0;

    return w.y / sinTheta(w);
}

static inline double sin2Phi(const Vec &w) {
    double s = sinPhi(w);
    return s * s;
}

// Rational approximation of inverse error function
// This method is borrowed from PBRT v3 
static inline double erfinv(double x) {
    double w, p;
    x = clamp(x, -0.99999f, 0.99999f);
    w = -std::log((1 - x) * (1 + x));
    if (w < 5) {
        w = w - 2.5f;
        p = 2.81022636e-08f;
        p = 3.43273939e-07f + p * w;
        p = -3.5233877e-06f + p * w;
        p = -4.39150654e-06f + p * w;
        p = 0.00021858087f + p * w;
        p = -0.00125372503f + p * w;
        p = -0.00417768164f + p * w;
        p = 0.246640727f + p * w;
        p = 1.50140941f + p * w;
    } else {
        w = std::sqrt(w) - 3;
        p = -0.000200214257f;
        p = 0.000100950558f + p * w;
        p = 0.00134934322f + p * w;
        p = -0.00367342844f + p * w;
        p = 0.00573950773f + p * w;
        p = -0.0076224613f + p * w;
        p = 0.00943887047f + p * w;
        p = 1.00167406f + p * w;
        p = 2.83297682f + p * w;
    }
    return p * x;
}

// Rational approximation of error function
// This method is borrowed from PBRT v3
static inline double erf(double x) {
    // constants
    double a1 = 0.254829592f;
    double a2 = -0.284496736f;
    double a3 = 1.421413741f;
    double a4 = -1.453152027f;
    double a5 = 1.061405429f;
    double p = 0.3275911f;

    // Save the sign of x
    int sign = 1;
    if (x < 0) sign = -1;
    x = std::abs(x);

    // A&S formula 7.1.26
    double t = 1 / (1 + p * x);
    double y =
        1 -
        (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

    return sign * y;
}

}  // anonymous namespace


// ----------------------------------------------------------------------------
// The interface for microfacet distribution
// ----------------------------------------------------------------------------
class MicrofacetDistribution {
public:
    MicrofacetDistribution(double alphax, double alphay, bool sampleVisible)
        : alphax_{alphax}
        , alphay_{alphay}
        , sampleVisible_{sampleVisible} {
    }

    virtual ~MicrofacetDistribution() {
    }

    // Sample microfacet normal
    virtual Vec sampleWm(const Vec &wo) const = 0;

    // The lambda function, which appears in the slope-space sampling
    virtual double lambda(const Vec &wo) const = 0;

    // Smith's masking function
    double G1(const Vec &wo) const {
        return 1.0 / (1.0 + lambda(wo));
    }

    // Smith's masking-shadowing function
    double G(const Vec &wo, const Vec &wi) const {
        return G1(wo) * G1(wi);
    }

    // Weighting function
    double weight(const Vec &wo, const Vec &wi, const Vec &wm) const {
        if (!sampleVisible_) {
            return std::abs(wi.dot(wm)) * G(wo, wi) /
                   std::max(1.0e-8, std::abs(cosTheta(wi) * cosTheta(wm)));
        } else {
            return G1(wo);
        }
    }

    // Private parameters
    double alphax_, alphay_;
    bool sampleVisible_;
};

// ----------------------------------------------------------------------------
// Beckmann distribution
// ----------------------------------------------------------------------------
class BeckmannDistribution : public MicrofacetDistribution {
public:
    BeckmannDistribution(double alphax, double alphay, bool sampleVisible)
        : MicrofacetDistribution{alphax, alphay, sampleVisible} {
    }

    Vec sampleWm(const Vec &wo) const override {
        const double U1 = random.next();
        const double U2 = random.next();

        if (!sampleVisible_) {
            // Sample half-vector without taking normal visibility into consideration
            double tan2Theta, phi;
            if (alphax_ == alphay_) {
                // Isotropic case
                // Following smapling formula can be found in [Walter et al. 2007]
                // "Microfacet Models for Refraction through Rough Surfaces"
                double alpha = alphax_;
                double logSample = std::log(1.0 - U1);
                tan2Theta = -alpha * alpha * logSample;
                phi = 2.0 * Pi * U2;
            } else {
                // Anisotropic case
                // Following sampling strategy is analytically derived from 
                // P.15 of [Heitz et al. 2014]
                // "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs"
                double logSample = std::log(1.0 - U1);
                phi = std::atan(alphay_ / alphax_ * std::tan(2.0 * Pi * U2 + 0.5 * Pi));
                if (U2 > 0.5) {
                    phi += Pi;
                }

                double sinPhi = std::sin(phi);
                double cosPhi = std::cos(phi);
                double alphax2 = alphax_ * alphax_;
                double alphay2 = alphay_ * alphay_;
                tan2Theta = -logSample / (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
            }

            // Compute normal direction from angle information sampled above
            double cosTheta = 1.0 / std::sqrt(1.0 + tan2Theta);
            double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
            Vec wm(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);

            if (wm.z < 0.0) {
                wm = -wm;
            }

            return wm;
        } else {
            // Sample microfacet normals by considering only visible normals
            bool flip = wo.z < 0.0;
            Vec wm = sampleBeckmann(flip ? -wo : wo, alphax_, alphay_);

            if (wm.z < 0.0) {
                wm = -wm;
            }

            return wm;
        }
    }

    double lambda(const Vec &wo) const override {
        using ::erf;

        double cosThetaO = cosTheta(wo);
        double sinThetaO = sinTheta(wo);
        double absTanThetaO = std::abs(tanTheta(wo)); 
        if (std::isinf(absTanThetaO)) return 0.;

        double alpha = std::sqrt(cosThetaO * cosThetaO * alphax_ * alphax_ +
                                 sinThetaO * sinThetaO * alphay_ * alphay_);
        double a = 1.0 / (alpha * absTanThetaO);
        return (erf(a) - 1.0) / 2.0 + std::exp(-a * a) / (2.0 * a * std::sqrt(Pi) + 1.0e-6);
    }

    // Sample microfacet normal through slope distribution
    static Vec sampleBeckmann(const Vec &wi, double alphax, double alphay) {
        // 1. stretch wi
        Vec wiStretched = Vec(alphax * wi.x, alphay * wi.y, wi.z).normalize();

        // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
        double slopex, slopey;
        sampleBeckman_P22_11(cosTheta(wiStretched), &slopex, &slopey);

        // 3. rotate
        double tmp = cosPhi(wiStretched) * slopex - sinPhi(wiStretched) * slopey;
        slopey = sinPhi(wiStretched) * slopex + cosPhi(wiStretched) * slopey;
        slopex = tmp;

        // 4. unstretch
        slopex = alphax * slopex;
        slopey = alphay * slopey;

        // 5. compute normal
        return Vec(-slopex, -slopey, 1.0).normalize();    
    }

    // Sample slope distribution with alphax and alphay equal to one
    static void sampleBeckman_P22_11(double cosThetaI, double *slopex, double *slopey) {
        using ::erf;

        const double U1 = random.next();
        const double U2 = random.next();

        // The special case where the ray comes from normal direction
        // The following sampling is equivarent to the sampling of
        // microfacet normals (not slopes) on isotropic rough surface
        if (cosThetaI > 0.9999) {
            double r = std::sqrt(-std::log(1.0f - U1));
            double sinPhi = std::sin(2 * Pi * U2);
            double cosPhi = std::cos(2 * Pi * U2);
            *slopex = r * cosPhi;
            *slopey = r * sinPhi;
            return;
        }

        double sinThetaI = std::sqrt(std::max((double)0, (double)1 - cosThetaI * cosThetaI));
        double tanThetaI = sinThetaI / cosThetaI;
        double cotThetaI = 1 / tanThetaI;
        const double erfCotThetaI = erf(cotThetaI);

        // Initialize lower and higher bounds for bisection method
        double lower = -1.0;
        double higher = 1.0;

        // Normalization factor for the CDF
        // This is equivarent to (G1 / 2)^{-1}
        static const double SQRT_PI_INV = 1.f / std::sqrt(Pi);
        double normalization = 1.0 / (1.0 + erfCotThetaI + SQRT_PI_INV * tanThetaI * std::exp(-cotThetaI * cotThetaI));

        // The following bisection method acquires "erf(x_m)"
        int it = 0;
        double samplex = std::max(U1, 1e-6);
        while (std::abs(higher - lower) > 1.0e-6) {
            // Bisection
            double mid = 0.5 * (lower + higher);

            // Evaluation for current estimation
            double x_m = erfinv(mid);
            double value = normalization * (1 + mid + SQRT_PI_INV * tanThetaI * std::exp(-x_m * x_m)) - samplex;

            /* Update bisection intervals */
            if (value > 0.0) {
                higher = mid;
            } else {
                lower = mid;
            }
        }

        // Compute slopex from erf(x_m) given by above bisection method
        *slopex = erfinv(lower);

        // Sample y_m
        *slopey = erfinv(2.0f * std::max(U2, 1e-6) - 1.0);
    }
};


// ----------------------------------------------------------------------------
// GGX distribution
// ----------------------------------------------------------------------------
class GGXDistribution : public MicrofacetDistribution {
public:
    GGXDistribution(double alphax, double alphay, bool sampleVisible)
        : MicrofacetDistribution{alphax, alphay, sampleVisible} {
    }

    Vec sampleWm(const Vec &wo) const override {
        const double U1 = random.next();
        const double U2 = random.next();

        if (!sampleVisible_) {
            double tan2Theta, phi;
            if (alphax_ == alphay_) {
                // Isotropic case
                double alpha = alphax_;
                tan2Theta = alpha * alpha * U1 / (1.0 - U1);
                phi = 2.0 * Pi * U2;                
            } else {
                // Anisotropic case
                phi = std::atan(alphay_ / alphax_ * std::tan(2.0 * Pi * U2 + 0.5 * Pi));
                if (U2 > 0.5) {
                    phi += Pi;
                }

                double sinPhi = std::sin(phi);
                double cosPhi = std::cos(phi);
                double alphax2 = alphax_ * alphax_;
                double alphay2 = alphay_ * alphay_;
                double alpha2 = 1.0 / (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
                tan2Theta = U1 / (1.0 - U1) * alpha2;
            }

            // Compute normal direction from angle information sampled above
            double cosTheta = 1.0 / std::sqrt(1.0 + tan2Theta);
            double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
            Vec wm = Vec(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);

            if (wm.z < 0.0) {
                wm = -wm;
            }

            return wm;
        } else {
            // Sample microfacet normals by considering only visible normals
            bool flip = wo.z < 0.0;
            Vec wm = sampleGGX(flip ? -wo : wo, alphax_, alphay_);

            if (wm.z < 0.0) {
                wm = -wm;
            }

            return wm;            
        }
    }

    static Vec sampleGGX(const Vec &wi, double alphax, double alphay) {
        // 1. stretch wi
        Vec wiStretched = Vec(alphax * wi.x, alphay * wi.y, wi.z).normalize();

        // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
        double slopex, slopey;
        sampleGGX_P22_11(cosTheta(wiStretched), &slopex, &slopey);

        // 3. rotate
        double tmp = cosPhi(wiStretched) * slopex - sinPhi(wiStretched) * slopey;
        slopey = sinPhi(wiStretched) * slopex + cosPhi(wiStretched) * slopey;
        slopex = tmp;

        // 4. unstretch
        slopex = alphax * slopex;
        slopey = alphay * slopey;

        // 5. compute normal
        return Vec(-slopex, -slopey, 1.0).normalize();       
    }

    static void sampleGGX_P22_11(double cosThetaI, double *slopex, double *slopey) {
        const double U1 = random.next();
        double U2 = random.next();

        // The special case where the ray comes from normal direction
        // The following sampling is equivarent to the sampling of
        // microfacet normals (not slopes) on isotropic rough surface
        if (cosThetaI > 0.9999) {
            double r = std::sqrt(U1 / (1.0 - U1));
            double sinPhi = std::sin(2 * Pi * U2);
            double cosPhi = std::cos(2 * Pi * U2);
            *slopex = r * cosPhi;
            *slopey = r * sinPhi;
            return;
        }

        double sinThetaI = std::sqrt(std::max(0.0, 1.0 - cosThetaI * cosThetaI));
        double tanThetaI = sinThetaI / cosThetaI;
        double a = 1.0 / tanThetaI;
        double G1 = 2.0 / (1.0 + std::sqrt(1.0 + 1.0 / (a * a)));

        // Sample slopex
        double A = 2.0 * U1 / G1 - 1.0;
        double B = tanThetaI;
        double tmp = std::min(1.0 / (A * A - 1.0), 1.0e12);

        double D = std::sqrt(B * B * tmp * tmp - (A * A - B * B) * tmp);
        double slopex1 = B * tmp - D;
        double slopex2 = B * tmp + D;
        *slopex = (A < 0.0 || slopex2 > 1.0 / tanThetaI) ? slopex1 : slopex2;

        // Sample slopey
        double S;
        if (U2 > 0.5) {
            S = 1.0;
            U2 = 2.0 * (U2 - 0.5);
        } else {
            S = -1.0;
            U2 = 2.0 * (0.5 - U2);
        }

        double z = (U2 * (U2 * (U2 * 0.27385f - 0.73369f) + 0.46341f)) /
                   (U2 * (U2 * (U2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
        *slopey = S * z * std::sqrt(1.0 + (*slopex) * (*slopex));
    }

    double lambda(const Vec &wo) const override {
        double absTanThetaO = std::abs(tanTheta(wo)); 
        if (std::isinf(absTanThetaO)) return 0.;
        
        double alpha =
            std::sqrt(cos2Phi(wo) * alphax_ * alphax_ + sin2Phi(wo) * alphay_ * alphay_);
        double alpha2Tan2Theta = (alpha * absTanThetaO) * (alpha * absTanThetaO);
        return (-1.0 + std::sqrt(1.f + alpha2Tan2Theta)) / 2;
    }
};

#endif  // _MICROFACET_H_
