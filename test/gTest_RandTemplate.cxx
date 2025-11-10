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