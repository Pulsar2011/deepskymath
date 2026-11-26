#include <gtest/gtest.h>
#include <DSTmath/DSTmath.h>
#include <DSTmath/DSTarray.h>
#include <DSTmath/DSTstatistic.h>


#include <string>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <type_traits>
#include <cmath>
#include <vector>
#include <algorithm>

template <typename T>
class StatsTest : public ::testing::Test
{
    protected:
        using value_type = T;
        static size_t N() { return 10000; }
};

using StatTypes = ::testing::Types<uint8_t,int8_t,int16_t,uint16_t,int32_t,uint32_t,int64_t,uint64_t,float,double>;
//using StatTypes = ::testing::Types<double>;
TYPED_TEST_SUITE(StatsTest, StatTypes);

TYPED_TEST(StatsTest, positive_uniform_distribution)
{
    using T = typename TestFixture::value_type;

    const size_t N = TestFixture::N();
    const T min = 0;
    const T max = 100;

    DST::Math::MaskedArray<T> array(N);
    for(size_t i=0; i < N; ++i)
        array.data()[i] = DST::Math::sample_uniform_distribution<T>(max);

    T dmin = array.min();
    T dmax = array.max();


    EXPECT_GE(dmin, min);
    EXPECT_LE(dmax, max);

    // KS goodness-of-fit without Boost
    constexpr double alpha = 1e-6; // very low alpha to avoid flakes in CI
    std::vector<double> xs(N);
    if constexpr (std::is_integral_v<T>) {
        // Map discrete uniform [min,max] to [0,1] using bin midpoints
        const double k = static_cast<double>(static_cast<long long>(max) - static_cast<long long>(min) + 1);
        for (size_t i=0;i<N;++i)
            xs[i] = (static_cast<double>(array.data()[i] - min) + 0.5) / k;
        double p = DST::Math::ks_uniform_pvalue(xs, 0.0, 1.0);
        EXPECT_GT(p, alpha);
    } else {
        for (size_t i=0;i<N;++i)
            xs[i] = static_cast<double>(array.data()[i]);
        double p = DST::Math::ks_uniform_pvalue(xs, static_cast<double>(min), static_cast<double>(max));
        EXPECT_GT(p, alpha);
    }
}

// Build Poisson CDF table up to kmax via stable recurrence
static std::vector<double> poisson_cdf_table(double lambda, int kmax)
{
    std::vector<double> F(std::max(kmax + 1, 0));
    if (kmax < 0) return F;
    double pk = std::exp(-lambda); // P(K=0)
    double c  = pk;
    F[0] = c;
    for (int k = 1; k <= kmax; ++k) {
        pk *= lambda / static_cast<double>(k);
        c  += pk;
        F[k] = c;
    }
    if (!F.empty()) F.back() = std::min(1.0, F.back());
    return F;
}

// KS D statistic for discrete Poisson using ECDF
static double ks_poisson_D(const std::vector<int>& xs, double lambda)
{
    if (xs.empty()) return 0.0;
    const int kmax = *std::max_element(xs.begin(), xs.end());
    auto F = poisson_cdf_table(lambda, kmax);

    std::vector<size_t> counts(kmax + 1, 0);
    for (int v : xs) {
        if (v >= 0) counts[std::min(v, kmax)]++;
    }

    size_t cum = 0;
    const double n = static_cast<double>(xs.size());
    double D = 0.0;
    for (int k = 0; k <= kmax; ++k) {
        cum += counts[k];
        const double Fe = cum / n;                 // empirical CDF at k
        const double Fp = F[k];                    // Poisson CDF at k
        D = std::max(D, std::fabs(Fe - Fp));
    }
    return D;
}

TYPED_TEST(StatsTest, uniform_distribution)
{
    using T = typename TestFixture::value_type;

    const size_t N = TestFixture::N();
    const T min = 0;
    const T max = 100;

    DST::Math::MaskedArray<T> array(N);
    for(size_t i=0; i < N; ++i)
        array.data()[i] = DST::Math::sample_uniform_distribution<T>(min,max);

    T dmin = array.min();
    T dmax = array.max();


    EXPECT_GE(dmin, min);
    EXPECT_LE(dmax, max);

    // KS goodness-of-fit without Boost
    constexpr double alpha = 1e-6; // very low alpha to avoid flakes in CI
    std::vector<double> xs(N);
    if constexpr (std::is_integral_v<T>) {
        // Map discrete uniform [min,max] to [0,1] using bin midpoints
        const double k = static_cast<double>(static_cast<long long>(max) - static_cast<long long>(min) + 1);
        for (size_t i=0;i<N;++i)
            xs[i] = (static_cast<double>(array.data()[i] - min) + 0.5) / k;
        double p = DST::Math::ks_uniform_pvalue(xs, 0.0, 1.0);
        EXPECT_GT(p, alpha);
    } else {
        for (size_t i=0;i<N;++i)
            xs[i] = static_cast<double>(array.data()[i]);
        double p = DST::Math::ks_uniform_pvalue(xs, static_cast<double>(min), static_cast<double>(max));
        EXPECT_GT(p, alpha);
    }
}

static inline double normal_cdf(double z)
{
    return 0.5 * (1.0 + std::erf(z / std::sqrt(2.0)));
}


TEST(DistributionTest, normalDistribution)
{
    const size_t N = 100000;
    const double mean = 50.0;
    const double sigma = 5.0;

    DST::Math::MaskedArray<double> array(N);
    for(size_t i=0; i < N; ++i)
        array.data()[i] = DST::Math::sample_normal_distribution(mean,sigma);

    double dmean = array.mean();
    double drms = std::sqrt(array.variance());
    double drmse= array.variance_error()/(2.*drms);

    EXPECT_NEAR(dmean, mean, 4.0* std::abs(drms)/sqrt(static_cast<double>(N)));
    EXPECT_NEAR(drms, sigma, 3.0* drmse);

    // KS on transformed Uniform(0,1): u_i = Phi((x_i - μ)/σ)
    std::vector<double> u(N);
    for (size_t i = 0; i < N; ++i) {
        const double z = (array.data()[i] - mean) / sigma;
        u[i] = normal_cdf(z);
    }
    constexpr double alpha = 1e-6;
    double p = DST::Math::ks_uniform_pvalue(u, 0.0, 1.0);
    EXPECT_GT(p, alpha);

    EXPECT_ANY_THROW(DST::Math::sample_normal_distribution(0,-1*sigma));
    EXPECT_ANY_THROW(DST::Math::sample_normal_distribution(0,0));
}

TEST(DistributionTest, poissonDistribution)
{
    const size_t N = 100000;
    const double _l = 10;

    DST::Math::MaskedArray<double> array(N);
    for(size_t i=0; i < N; ++i)
        array.data()[i] = DST::Math::sample_poisson_distribution(_l);

    // Method-of-moments checks
    const double dmean = array.mean();
    const double drms  = std::sqrt(array.variance());
    const double mean_se = std::sqrt(_l / static_cast<double>(N));
    const double drmse = array.variance_error() / (2.0 * drms);

    EXPECT_NEAR(dmean, _l, 4.0 * mean_se);
    EXPECT_NEAR(drms,  std::sqrt(_l), 3.0 * drmse);

    // KS with DKW bound (no Boost)
    std::vector<int> xs(N);
    for (size_t i = 0; i < N; ++i) xs[i] = static_cast<int>(array.data()[i]+0.5);

    const double D = ks_poisson_D(xs, _l);
    constexpr double alpha = 1e-6;
    const double eps = std::sqrt(std::log(2.0 / alpha) / (2.0 * static_cast<double>(N)));
    EXPECT_LE(D, eps);

    EXPECT_ANY_THROW(DST::Math::sample_poisson_distribution(-1.0));
}

// Composite Simpson rule of f on [a, b] with n intervals (n made even)
template <typename F>
static double simpson_integral(F&& f, double a, double b, size_t n)
{
    n += n % 2;
    const double h = (b - a) / static_cast<double>(n);
    double s = f(a) + f(b);
    for (size_t k = 1; k < n; ++k)
        s += ((k % 2) ? 4.0 : 2.0) * f(a + h * static_cast<double>(k));
    return s * h / 3.0;
}

// Standard Landau density from its definition, independent of the DENLAN coefficients:
//   phi(l) = (1/pi) int_0^inf exp(-t ln t - l t) sin(pi t) dt
// Valid for l >= -3: below, the integrand reaches exp(e^{-l-1}) and cancels to a tiny result.
static double landau_density_definition(double l)
{
    const double pi = std::acos(-1.0);
    const double T  = (l <= 1.0) ? 60.0 : 1.0 + 60.0 / l;   // integrand negligible beyond T
    auto f = [l, pi](double t) { return (t > 0.0) ? std::exp(-t * std::log(t) - l * t) * std::sin(pi * t) : 0.0; };
    return simpson_integral(f, 0.0, T, 200000) / pi;
}

// Upper tail of the standard Landau from the same definition, independent of DENLAN and of the table:
//   1 - Phi(l) = (1/pi) int_0^inf exp(-t ln t - l t) sin(pi t) / t dt   (its derivative in l is -phi(l), and it -> 0 as l -> inf)
// Valid for l >= -3, as above.
static double landau_upper_tail_definition(double l)
{
    const double pi = std::acos(-1.0);
    const double T  = (l <= 1.0) ? 60.0 : 60.0 / l;   // l > 1: integrand below e^(1/e - 60) beyond T
    auto f = [l, pi](double t) { return (t > 0.0) ? std::exp(-t * std::log(t) - l * t) * std::sin(pi * t) / t : pi; };
    return simpson_integral(f, 0.0, T, 200000) / pi;
}

TEST(DistributionTest, landauDensity)
{
    using S = DST::Math::standardLandau;

    // DENLAN against the definition, one point in each fitted branch (6 significant digits expected)
    for (double l : {-3.0, -1.5, S::lambda_mode, 0.5, 2.0, 8.0, 20.0, 100.0, 500.0})
    {
        const double ref = landau_density_definition(l);
        EXPECT_NEAR(S::density(l), ref, 1e-5 * ref) << "lambda = " << l;
    }

    // Continuity where DENLAN changes formula: a coefficient typo shows up as a jump
    for (double b : {-5.5, -1.0, 1.0, 5.0, 12.0, 50.0, 300.0})
    {
        const double right = S::density(b);
        const double left  = S::density(std::nextafter(b, -std::numeric_limits<double>::infinity()));
        EXPECT_NEAR(left, right, 1e-5 * right) << "branch edge lambda = " << b;
    }

    // Mode
    EXPECT_GT(S::density(S::lambda_mode), S::density(S::lambda_mode - 1e-2));
    EXPECT_GT(S::density(S::lambda_mode), S::density(S::lambda_mode + 1e-2));

    // Normalisation: Simpson on [-6, 20], in ln(lambda) on [20, 1e4], analytical tail above 1e4
    const double core = simpson_integral(S::density, -6.0, 20.0, 26000);
    const double tail = simpson_integral([](double s) { const double l = std::exp(s); return S::density(l) * l; },
                                         std::log(20.0), std::log(1e4), 20000);
    EXPECT_NEAR(core + tail + S::upperTail(1e4), 1.0, 1e-5);
}

TEST(DistributionTest, landauCDF)
{
    using S = DST::Math::standardLandau;
    const S& standard = S::instance();

    // Table against a direct integration of the density
    for (double l : {-2.0, S::lambda_mode, 0.0, 3.0, 15.0})
        EXPECT_NEAR(standard.cdf(l), simpson_integral(S::density, -6.0, l, 20000), 1e-5) << "lambda = " << l;

    // Table against the definition: core, uniform and geometric parts, analytical tail above the last node (1e6)
    for (double l : {-3.0, -2.0, S::lambda_mode, 0.0, 1.0, 3.0, 10.0, 50.0, 300.0, 1e3, 1e4, 1e5, 2e6})
    {
        const double ref = landau_upper_tail_definition(l);
        EXPECT_NEAR(1.0 - standard.cdf(l), ref, 1e-5 * ref) << "lambda = " << l;
    }

    // The quantile inverts the CDF: left tail, core, tabulated tail and analytical tail (above 1 - 1e-6)
    for (double u : {1e-8, 1e-3, 0.1, 0.5, 0.9, 0.999, 1.0 - 1e-7})
        EXPECT_NEAR(standard.cdf(standard.quantile(u)), u, 1e-10) << "u = " << u;

    // Limits
    EXPECT_EQ(standard.cdf(-10.0), 0.0);
    EXPECT_TRUE(std::isinf(standard.quantile(1.0)));

    const DST::Math::Landau_distribution landau(0.0, 1.0);
    EXPECT_DOUBLE_EQ(landau.getIntegral(), 1.0);
}

// Laplace transform of the standard Landau (Landau 1944; Kolbig & Schorr 1984):
//   int exp(-s l) phi(l) dl = s^s,  s > 0
// and, integrating by parts, int exp(-s l) Phi(l) dl = s^(s-1).
// Exact values, independent of DENLAN and of the table. Small s weights the heavy right tail,
// large s the double-exponential left tail (s = 3: integrand peaks near l = -2.3).
TEST(DistributionTest, landauLaplaceTransform)
{
    using S = DST::Math::standardLandau;
    const S& standard = S::instance();

    for (double s : {0.1, 0.5, 1.0, 2.0, 3.0})
    {
        const double a = -8.0;       // phi(-8) ~ exp(-e^7): nothing below
        const double b = 80.0 / s;   // exp(-s b) = e^-80: nothing above
        const size_t n = 200000;

        const double ref_phi = std::pow(s, s);
        const double L_phi   = simpson_integral([s](double l) { return std::exp(-s * l) * S::density(l); }, a, b, n);
        EXPECT_NEAR(L_phi, ref_phi, 1e-5 * ref_phi) << "density, s = " << s;

        const double ref_Phi = std::pow(s, s - 1.0);
        const double L_Phi   = simpson_integral([s, &standard](double l) { return std::exp(-s * l) * standard.cdf(l); }, a, b, n);
        EXPECT_NEAR(L_Phi, ref_Phi, 1e-5 * ref_Phi) << "cumulative, s = " << s;
    }
}

TEST(DistributionTest, landauDistribution)
{
    const size_t N = 100000;
    const double mpv = 50.0;
    const double xi  = 5.0;
    constexpr double alpha = 1e-6;

    DST::Math::Landau_distribution landau(mpv, xi);
    EXPECT_DOUBLE_EQ(landau.getMPV(), mpv);

    const double F_low = landau.cdfAt(landau.lowEdge());
    const double F_up  = landau.cdfAt(landau.upEdge());

    std::vector<double> x(N);
    for (size_t i = 0; i < N; ++i)
        x[i] = landau();

    // Every sample in the range
    EXPECT_GE(*std::min_element(x.begin(), x.end()), landau.lowEdge());
    EXPECT_LE(*std::max_element(x.begin(), x.end()), landau.upEdge());

    // Median (the Landau has neither mean nor variance): standard error 1 / (2 f(m) sqrt(N))
    const double m = landau.quantile(F_low + 0.5 * (F_up - F_low));
    std::vector<double> sorted(x);
    std::nth_element(sorted.begin(), sorted.begin() + N/2, sorted.end());
    const double median_se = (F_up - F_low) / (2.0 * landau(m) * std::sqrt(static_cast<double>(N)));
    EXPECT_NEAR(sorted[N/2], m, 4.0 * median_se);

    // Heavy tail: share of samples above mpv + 20 xi, binomial
    const double xT = mpv + 20.0 * xi;
    const double pT = (F_up - landau.cdfAt(xT)) / (F_up - F_low);
    const double nT = static_cast<double>(std::count_if(x.begin(), x.end(), [xT](double v) { return v > xT; }));
    EXPECT_NEAR(nT / static_cast<double>(N), pT, 4.0 * std::sqrt(pT * (1.0 - pT) / static_cast<double>(N)));

    // KS on transformed Uniform(0,1): u_i = (F(x_i) - F(low)) / (F(up) - F(low))
    std::vector<double> u(N);
    for (size_t i = 0; i < N; ++i)
        u[i] = (landau.cdfAt(x[i]) - F_low) / (F_up - F_low);
    EXPECT_GT(DST::Math::ks_uniform_pvalue(u, 0.0, 1.0), alpha);

    // Truncated range: samples stay inside and follow the renormalised CDF
    landau.setRange(mpv - 2.0 * xi, mpv + 10.0 * xi);
    const double G_low = landau.cdfAt(landau.lowEdge());
    const double G_up  = landau.cdfAt(landau.upEdge());
    for (size_t i = 0; i < N; ++i)
    {
        x[i] = landau();
        u[i] = (landau.cdfAt(x[i]) - G_low) / (G_up - G_low);
    }
    EXPECT_GE(*std::min_element(x.begin(), x.end()), landau.lowEdge());
    EXPECT_LE(*std::max_element(x.begin(), x.end()), landau.upEdge());
    EXPECT_GT(DST::Math::ks_uniform_pvalue(u, 0.0, 1.0), alpha);

    // The copy keeps the range
    DST::Math::Landau_distribution copy(landau);
    EXPECT_DOUBLE_EQ(copy.lowEdge(), landau.lowEdge());
    EXPECT_DOUBLE_EQ(copy.upEdge(),  landau.upEdge());

    EXPECT_ANY_THROW(DST::Math::Landau_distribution(mpv, -xi));
    EXPECT_ANY_THROW(DST::Math::Landau_distribution(mpv, 0.0));
    EXPECT_ANY_THROW(landau.setRange(mpv + xi, mpv - xi));                 // a > b
    EXPECT_ANY_THROW(landau.setRange(mpv - 40.0 * xi, mpv - 30.0 * xi));   // lambda < -6: no probability
}