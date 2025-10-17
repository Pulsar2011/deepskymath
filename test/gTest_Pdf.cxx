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


class test_pdf : public DST::Math::pdf
{
public:
    test_pdf() : DST::Math::pdf("test PDF") {};
    ~test_pdf() {clear();};


    double gauss_fcn(const double& x, const DST::Math::pdf_param& p) const
    {
        if(p.size() != 3)
            throw std::invalid_argument("\033[43;31m[test_pdf::operator()] Invalid number of parameter for gaussian function (3 expected)\033[0m");
        
        double A = p[0];
        double m = p[1];
        double s = p[2];
        
        return A*(1./(s*sqrt(2.*DST::Math::MathCore::Pi())))*exp(-0.5*((x-m)/s)*((x-m)/s));
    }

    // Expose protected API for testing
    using DST::Math::pdf::brent;
};

TEST(pdf,ctor)
{
    test_pdf p = test_pdf();

    for(size_t i=0; i < 100; ++i)
        EXPECT_NEAR(p.gen(),0, std::numeric_limits<double>::epsilon());
}

TEST(pdf,managment)
{
    test_pdf p = test_pdf();


    p.add_pdf(static_cast<DST::Math::pdf_function>(&test_pdf::gauss_fcn),{1.,0,1});
    p.setRange(-3,3);

    for(double x=-3; x <= 3; x+=0.1)
        EXPECT_NEAR(p(x),p.gauss_fcn(x,{1.,0,1}),1e-5);

    p.add_pdf(static_cast<DST::Math::pdf_function>(&test_pdf::gauss_fcn),{1.,0,4});
    p.setRange(-10,10);

    for(double x=-10; x <= 10; x+=0.1)
        EXPECT_NEAR(p(x),p.gauss_fcn(x,{1.,0,1})+p.gauss_fcn(x,{1.,0,4}),1e-5);

    p.clear();
    for(size_t i=0; i < 100; ++i)
       EXPECT_NEAR(p.gen(),0, std::numeric_limits<double>::epsilon());
}

TEST(pdf,modifier)
{
    test_pdf p = test_pdf();
    DST::Math::MaskedArray<double> a(100);

    p.add_pdf(static_cast<DST::Math::pdf_function>(&test_pdf::gauss_fcn),{0.75,0,1});
    EXPECT_THROW(p.setRange(3,3), std::invalid_argument);
    EXPECT_THROW(p.setRange(3,-3), std::invalid_argument);
}

TEST(pdf,normalization)
{
    test_pdf p = test_pdf();
    DST::Math::MaskedArray<double> a(100);

    p.add_pdf(static_cast<DST::Math::pdf_function>(&test_pdf::gauss_fcn),{10.75,0,1});
    p.add_pdf(static_cast<DST::Math::pdf_function>(&test_pdf::gauss_fcn),{5.25,0,4});
    
    p.setRange(-10,10);

    for(double x=-10; x <= 10; x+=0.1)
        EXPECT_NEAR(p(x),p.gauss_fcn(x,{10.75,0,1})+p.gauss_fcn(x,{5.25,0,4}),1e-5);

    size_t nStep = 100000;
    EXPECT_NEAR(p.getIntegral(nStep),10.75+5.25,1./(2.*static_cast<double>(nStep)));

    EXPECT_NEAR(p.NormalisationFactor(),1.,1./(2.*static_cast<double>(nStep)));

    p.Normalize(nStep);
    EXPECT_NEAR(p.NormalisationFactor(),1./(10.75+5.25),1./(2.*static_cast<double>(nStep)));
    EXPECT_NEAR(p.getIntegral(nStep),1.,1./(2.*static_cast<double>(nStep)));
}

TEST(pdf,MPV_and_derivatives)
{
    test_pdf p = test_pdf();
    DST::Math::MaskedArray<double> a(100);

    p.add_pdf(static_cast<DST::Math::pdf_function>(&test_pdf::gauss_fcn),{2,0.5,0.25});
    p.add_pdf(static_cast<DST::Math::pdf_function>(&test_pdf::gauss_fcn),{0.5,-0.5,0.15});
    
    p.setRange(-5,5);

    for(double x=-5; x <= 5; x+=0.01)
    {
        EXPECT_NEAR(p(x),p.gauss_fcn(x,{2,0.5,0.25})+p.gauss_fcn(x,{0.5,-0.5,0.15}),1e-5);
        EXPECT_NEAR(p.dfdx(x,1e-8),-0.5*2*((x-0.5)/(0.25*0.25)*p.gauss_fcn(x,{2,0.5,0.25})+(x+0.5)/(0.15*0.15)*p.gauss_fcn(x,{0.5,-0.5,0.15})),1e-5);
    }

    EXPECT_NEAR(p.getMPV(),0.5,1e-10)<<p(p.getMPV())<<" vs "<<p(0.5);
}

TEST(pdf,IO)
{
    test_pdf p = test_pdf();

    double mVal = 0.5;
    double sVal = 1.25;
    p.add_pdf(static_cast<DST::Math::pdf_function>(&test_pdf::gauss_fcn),{4,mVal,sVal});
    p.setRange(mVal-5*sVal,mVal+5*sVal);

    EXPECT_NO_THROW(p.dump(50));

    p.Normalize(100000);
    EXPECT_NO_THROW(p.dump(50));

    p.setRange(mVal-sVal,mVal+sVal);
    EXPECT_NO_THROW(p.dump(50));

}

//TRANSFORM GENERATED RANDOM SAMPLE INTO A UNIFORM CDF DISTRIBUTION
double cdf_lookup(const DST::Math::cdf_t& xs, double x)
{
    if (x <= xs.front().first) return 0.0;
    if (x >= xs.back().first ) return 1.0;
    
    DST::Math::cdf_t::const_iterator it = xs.cbegin();
    while(it != xs.cend() && it->first < x)
        ++it;
    
    size_t i = static_cast<size_t>(it - xs.cbegin());

    const double x0 = xs[i-1].first,  x1 = xs[i].first;
    const double F0 = xs[i-1].second, F1 = xs[i].second;
    const double t = (x - x0) / (x1 - x0);
    return F0 + t * (F1 - F0);

}

TEST(pdf,random_generation)
{
    test_pdf p = test_pdf();
    size_t nGen = 10000;

    p.add_pdf(static_cast<DST::Math::pdf_function>(&test_pdf::gauss_fcn),{1,0.5,0.25});
    
    p.setRange(-5,5);
    p.Normalize(100000);

    DST::Math::cdf_t xs = p.cdf(100);

    double xMin = xs.begin()->first; 
    double xMax = xs.rbegin()->first;
    //NORMALIZE RANGE OF CDF TO [0,1]
    for(DST::Math::cdf_t::iterator it = xs.begin(); it != xs.end(); ++it)
    {
        it->first -= xMin;
        it->first /= (xMax-xMin);
    }

    std::vector<double> histCdf;
    histCdf.reserve(xs.size());

    for(size_t i=0; i < nGen; ++i)
    {
        histCdf.push_back(cdf_lookup(xs,(p.gen()-xMin)/(xMax-xMin)));
    }
    
     // KS on Uniform(0,1) using your existing helper
     constexpr double alpha = 1e-4;
     double pval = DST::Math::ks_uniform_pvalue(histCdf, 0.0, 1.0);
     EXPECT_GT(pval, alpha);
}

TEST(pdf,error_handling)
{
    test_pdf p = test_pdf();

    EXPECT_THROW(p.getIntegral(), std::runtime_error);
    EXPECT_THROW(p.getIntegral(0,1), std::runtime_error);
    EXPECT_THROW(p.Normalize(1000), std::runtime_error);

    p.add_pdf(static_cast<DST::Math::pdf_function>(&test_pdf::gauss_fcn),{1.,0,1});
    EXPECT_THROW(p.setRange(1,1), std::invalid_argument);
    EXPECT_THROW(p.setRange(1,-1), std::invalid_argument);
    EXPECT_NO_THROW(p.setRange(-3,3));

    EXPECT_THROW(p.getIntegral(3,10), std::invalid_argument);
    EXPECT_THROW(p.getIntegral(-10,3), std::invalid_argument);
    EXPECT_THROW(p.getIntegral(3,-3), std::invalid_argument);
    EXPECT_THROW(p.getIntegral(3,-3,0), std::invalid_argument);
    EXPECT_THROW(p.getIntegral(0), std::invalid_argument);
    EXPECT_NO_THROW(p.getIntegral(-3,3));
    EXPECT_NO_THROW(p.getIntegral(-1,1));

    EXPECT_THROW(p.Normalize(0), std::invalid_argument);
    

    EXPECT_THROW(p.cdf(0), std::invalid_argument);
    EXPECT_THROW(p.cdf(100,0), std::invalid_argument);
    EXPECT_THROW(p.cdf(0,100), std::invalid_argument);
    EXPECT_NO_THROW(p.cdf(100,100));

    EXPECT_THROW(p.brent(-3,-0.5,1e-4), std::invalid_argument);
    EXPECT_THROW(p.brent(3,-3,1e-4), std::invalid_argument);
    EXPECT_NO_THROW(p.brent(-0.5,0.5,1e-4));
}

TEST(normal_distribution,inheritance)
{
    DST::Math::normal_distribution n = DST::Math::normal_distribution(5,1);
    n.setRange(-10,10);

    EXPECT_NEAR(n.getMPV(),5,1e-10);
    EXPECT_NEAR(n.getMPV(),n.DST::Math::pdf::getMPV(),1e-10);

    n.Mean(2.5);
    EXPECT_NEAR(n.getMPV(),2.5,1e-10);
    EXPECT_NEAR(n.getMPV(),n.DST::Math::pdf::getMPV(),1e-10);

    n.Sigma(0.5);
    EXPECT_NEAR(n.getMPV(),2.5,1e-10);
    EXPECT_NEAR(n.getMPV(),n.DST::Math::pdf::getMPV(),1e-10);

    n.setRange(-1,5);

    EXPECT_NEAR(n.getIntegral(-1,5,100000),1.0,1e-6);
    EXPECT_NEAR(n.getIntegral(-1,5,100000),n.getIntegral(0),1e-6);
    EXPECT_NEAR(n.getIntegral(-1,5,100000),n.DST::Math::pdf::getIntegral(-1,5,100000),1e-6);
    EXPECT_NEAR(n.getIntegral(),n.DST::Math::pdf::getIntegral(10000),1e-6);
    
    for(double ix = -1; ix <= 5; ix += 0.1)
    {
        double val = (1./(0.5*sqrt(2.*DST::Math::MathCore::Pi())))*exp(-0.5*((ix-2.5)/0.5)*((ix-2.5)/0.5));
        EXPECT_NEAR(n(ix),val,1e-5);
        EXPECT_NEAR(n.dfdx(ix,1e-8),-0.5*2*((ix-2.5)/(0.5*0.5))*val,1e-6);
        EXPECT_NEAR(n.dfdx(ix,1e-8),n.DST::Math::pdf::dfdx(ix,1e-8),1e-6);
    }

    DST::Math::cdf_t ccdf = n.cdf(1000);
    DST::Math::cdf_t pcdf = n.DST::Math::pdf::cdf(1000);
    
    for(size_t i=0; i < ccdf.size(); ++i)
    {
        EXPECT_NEAR(ccdf[i].first, pcdf[i].first, 1e-6);
        EXPECT_NEAR(ccdf[i].second, pcdf[i].second, 1e-6);
        EXPECT_NEAR(ccdf[i].second, 0.5*(std::erf((ccdf[i].first-2.5)/(0.5*sqrt(2)))-std::erf((-10-2.5)/(0.5*sqrt(2)))),1e-6);
    }

    double norm=n.NormalisationFactor();
    n.Normalize(100000);
    EXPECT_NEAR(n.NormalisationFactor(),norm,1e-6);

    n.dump(50);
}

TEST(normal_distribution,random_generation)
{
    DST::Math::normal_distribution p = DST::Math::normal_distribution(5,1);
    p.setRange(-10,10);
    size_t nGen = 10000;

    DST::Math::cdf_t xs = p.cdf(100);

    double xMin = xs.begin()->first; 
    double xMax = xs.rbegin()->first;
    //NORMALIZE RANGE OF CDF TO [0,1]
    for(DST::Math::cdf_t::iterator it = xs.begin(); it != xs.end(); ++it)
    {
        it->first -= xMin;
        it->first /= (xMax-xMin);
    }

    std::vector<double> histCdf;
    histCdf.reserve(xs.size());

    for(size_t i=0; i < nGen; ++i)
    {
        histCdf.push_back(cdf_lookup(xs,(p.gen()-xMin)/(xMax-xMin)));
    }
    
     // KS on Uniform(0,1) using your existing helper
     constexpr double alpha = 1e-4;
     double pval = DST::Math::ks_uniform_pvalue(histCdf, 0.0, 1.0);
     EXPECT_GT(pval, alpha);
}

TEST(normal_distribution,erroHandling)
{
    DST::Math::normal_distribution n = DST::Math::normal_distribution(5,1);
    n.setRange(-10,10);

    EXPECT_NO_THROW(n.getIntegral(-10,10,0));
    EXPECT_NO_THROW(n.getIntegral(0));
    EXPECT_THROW(n.getIntegral(10, -10), std::invalid_argument);
}

TEST(log_normal_distribution,inheritance)
{
    DST::Math::log_normal_distribution n = DST::Math::log_normal_distribution(5,1);
    double _PI2_ = 2*DST::Math::MathCore::Pi();

    EXPECT_NEAR(n.getMPV(),exp(log(5)-1*1),1e-10);

    n.Mean(2.5);
    EXPECT_NEAR(n.getMPV(),exp(log(2.5)-1*1),1e-10);

    n.Sigma(0.5);
    EXPECT_NEAR(n.getMPV(),exp(log(2.5)-0.5*0.5),1e-10);

    double minVal = 1e-6;
    double maxVal = 2.5 + 280.*0.5;

    n.setRange(minVal,maxVal);

    EXPECT_NEAR(n.getIntegral(minVal,maxVal,100000),1.0,1e-6);
    EXPECT_NEAR(n.getIntegral(minVal,maxVal,100000),n.getIntegral(0),1e-6);
    EXPECT_NEAR(n.getIntegral(minVal,maxVal,100000),n.DST::Math::pdf::getIntegral(minVal,maxVal,100000),1e-6);
    EXPECT_NEAR(n.getIntegral(),n.DST::Math::pdf::getIntegral(10000),1e-6);
    
    double dx=(maxVal - minVal)/10.;
    for(double ix = minVal; ix <= maxVal; ix += dx)
    {
        double val =1./(std::sqrt(0.5*0.5*_PI2_)*ix) * std::exp(-0.5*std::pow(std::log(ix)-2.5,2.)/(0.5*0.5));
        EXPECT_NEAR(n(ix),val,1e-5);
        EXPECT_NEAR(n.dfdx(ix,1e-8),-( 1 + (std::log(ix)-2.5)/(0.5*0.5) )*val/ix,1e-6);
        EXPECT_NEAR(n.dfdx(ix,1e-8),n.DST::Math::pdf::dfdx(ix,1e-8),1e-6);
    }

    DST::Math::cdf_t ccdf = n.cdf(50);
    DST::Math::cdf_t pcdf = n.DST::Math::pdf::cdf(50);
    
    for(size_t i=0; i < ccdf.size(); ++i)
    {
        EXPECT_NEAR(ccdf[i].first, pcdf[i].first, 1e-6);
        EXPECT_NEAR(ccdf[i].second, pcdf[i].second, 1e-6);
        EXPECT_NEAR(ccdf[i].second, 0.5*(std::erf((std::log(ccdf[i].first)-2.5)/(0.5*sqrt(2)))-std::erf((0-2.5)/(0.5*sqrt(2)))),1e-6);
    }

    double norm=n.NormalisationFactor();
    n.Normalize(100000);
    EXPECT_NEAR(n.NormalisationFactor(),norm,1e-6);

    n.dump(100);
}

TEST(log_normal_distribution,erroHandling)
{
    DST::Math::log_normal_distribution n = DST::Math::log_normal_distribution(5,1);
    EXPECT_THROW(n.setRange(-10,10), std::invalid_argument);

    n.setRange(1e-6,10);

    EXPECT_NO_THROW(n.getIntegral(1e-6,10,0));
    EXPECT_THROW(n.getIntegral(0,10,0), std::invalid_argument);
    EXPECT_NO_THROW(n.getIntegral(0));
    EXPECT_THROW(n.getIntegral(10, -10), std::invalid_argument);
}

TEST(log_normal_distribution,random_generation)
{
    DST::Math::log_normal_distribution p = DST::Math::log_normal_distribution(2.5,0.5);
    double minVal = 1e-6;
    double maxVal = 2.5 + 280.*0.5;

    p.setRange(minVal,maxVal);

    size_t nGen = 10000;

    DST::Math::cdf_t xs = p.cdf(100);

    double xMin = xs.begin()->first; 
    double xMax = xs.rbegin()->first;
    //NORMALIZE RANGE OF CDF TO [0,1]
    for(DST::Math::cdf_t::iterator it = xs.begin(); it != xs.end(); ++it)
    {
        it->first -= xMin;
        it->first /= (xMax-xMin);
    }

    std::vector<double> histCdf;
    histCdf.reserve(xs.size());

    for(size_t i=0; i < nGen; ++i)
    {
        histCdf.push_back(cdf_lookup(xs,(p.gen()-xMin)/(xMax-xMin)));
    }
    
     // KS on Uniform(0,1) using your existing helper
     constexpr double alpha = 1e-4;
     double pval = DST::Math::ks_uniform_pvalue(histCdf, 0.0, 1.0);
     EXPECT_GT(pval, alpha);
}

TEST(chi2_distribution,erroHandling)
{
    EXPECT_THROW(DST::Math::chi2_distribution(-5), std::invalid_argument);

    DST::Math::chi2_distribution n = DST::Math::chi2_distribution(5);
    EXPECT_THROW(n.setRange(-10,10), std::invalid_argument);

    n.setRange(1e-6,10);

    EXPECT_NO_THROW(n.getIntegral(1e-5,10,0));
    EXPECT_THROW(n.getIntegral(0,10,0),std::invalid_argument);
    EXPECT_NO_THROW(n.getIntegral(0));
    EXPECT_THROW(n.getIntegral(10, -10), std::invalid_argument);
}

TEST(chi2_distribution,inheritance)
{
    DST::Math::chi2_distribution n = DST::Math::chi2_distribution(5);


    EXPECT_NEAR(n.getMPV(),5-2,1e-10);
    EXPECT_NEAR(n.getMPV(),n.DST::Math::pdf::getMPV(),1e-6);

    n.NDF(2);
    EXPECT_NEAR(n.getMPV(),0,1e-10);
    EXPECT_NEAR(n.getMPV(),n.DST::Math::pdf::getMPV(),1e-6);

    n.NDF(8);
    EXPECT_NEAR(n.getMPV(),8-2,1e-10);
    EXPECT_NEAR(n.getMPV(),n.DST::Math::pdf::getMPV(),1e-6);

    double minVal = 1e-3;
    double maxVal = 3*n.get95Q();

    n.setRange(minVal,maxVal);

    EXPECT_NEAR(n.getIntegral(minVal,maxVal,100000),1.0,1e-6);
    EXPECT_NEAR(n.getIntegral(minVal,maxVal,100000),n.getIntegral(0),1e-6);
    EXPECT_NEAR(n.getIntegral(minVal,maxVal,100000),n.DST::Math::pdf::getIntegral(minVal,maxVal,100000),1e-6);
    
    double dx=(maxVal - minVal)/10.;
    for(double ix = minVal; ix <= maxVal; ix += dx)
    {
        double val =  1./(pow(2.,8./2.)*DST::Math::MathCore::Gamma(8./2.)) * pow(ix,(8./2.)-1.) * exp(-ix/2.);
        EXPECT_NEAR(n(ix),val,1e-5);
        EXPECT_NEAR(n.dfdx(ix,1e-8),((0.5*8.-1)/ix-1./2.)*val,1e-6);
        EXPECT_NEAR(n.dfdx(ix,1e-8),n.DST::Math::pdf::dfdx(ix,1e-8),1e-6);
    }

    DST::Math::cdf_t ccdf = n.cdf(50);
    DST::Math::cdf_t pcdf = n.DST::Math::pdf::cdf(50);

    
    for(size_t i=0; i < ccdf.size(); ++i)
    {
        EXPECT_NEAR(ccdf[i].first, pcdf[i].first, 1e-6);
        EXPECT_NEAR(ccdf[i].second, pcdf[i].second, 1e-6);
        EXPECT_NEAR(ccdf[i].second, DST::Math::MathCore::Gamma(8./2.,ccdf[i].first/2.) ,1e-6);
    }

    double norm=n.NormalisationFactor();
    n.Normalize(100000);
    EXPECT_NEAR(n.NormalisationFactor(),norm,1e-6);

    n.dump(100);
}

TEST(chi2_distribution,random_generation)
{
    DST::Math::chi2_distribution p = DST::Math::chi2_distribution(5);
    double minVal = 1e-3;
    double maxVal = 3*p.get95Q();

    p.setRange(minVal,maxVal);

    size_t nGen = 10000;

    DST::Math::cdf_t xs = p.cdf(100);

    double xMin = xs.begin()->first; 
    double xMax = xs.rbegin()->first;
    //NORMALIZE RANGE OF CDF TO [0,1]
    for(DST::Math::cdf_t::iterator it = xs.begin(); it != xs.end(); ++it)
    {
        it->first -= xMin;
        it->first /= (xMax-xMin);
    }

    std::vector<double> histCdf;
    histCdf.reserve(xs.size());

    for(size_t i=0; i < nGen; ++i)
    {
        histCdf.push_back(cdf_lookup(xs,(p.gen()-xMin)/(xMax-xMin)));
    }
    
     // KS on Uniform(0,1) using your existing helper
     constexpr double alpha = 1e-4;
     double pval = DST::Math::ks_uniform_pvalue(histCdf, 0.0, 1.0);
     EXPECT_GT(pval, alpha);
}

TEST(Cball_distribution,inheritance)
{
    double alpha = 10;
    double n = 3;
    double mean = 2.5;
    double sigma = 1.25;
    DST::Math::Cball_distribution p = DST::Math::Cball_distribution(alpha,n,mean,sigma);

    EXPECT_NEAR(p.getMPV(),mean,1e-10);

    alpha = 1.0;
    p.Alpha(alpha);
    EXPECT_NEAR(p.getMPV(),mean,1e-10);

    n=3;
    p.N(n);
    EXPECT_NEAR(p.getMPV(),mean,1e-10);

    mean = 1.5;
    p.Mean(mean);
    EXPECT_NEAR(p.getMPV(),mean,1e-10);

    sigma = 1.;
    p.Sigma(sigma);
    EXPECT_NEAR(p.getMPV(),mean,1e-10);

    double minVal = mean-30*sigma;
    double maxVal = mean+5*sigma;

    p.setRange(minVal,maxVal);

   // EXPECT_NEAR(p.getIntegral(minVal,maxVal,100000),1.0,1e-6);
   // EXPECT_NEAR(p.getIntegral(),1,1e-6);

    auto Cball_fn = [alpha,n,mean,sigma](double x)
    {
        double A = pow(n/std::abs(alpha),n) * exp(-0.5*alpha*alpha);
        double B = n/std::abs(alpha) - std::abs(alpha);
        double C = n/std::abs(alpha) * 1./( n - 1.) * exp(-0.5*alpha*alpha);
        double D = std::sqrt(DST::Math::MathCore::Pi()/2.) * (1.+std::erf(alpha/std::sqrt(2.)));
        double N = 1./sigma * 1./(C + D);

        double t = (x - mean) / sigma;
                    
        return  N* ((t <= -alpha) ? A*pow(B - t, -n) : exp(-0.5 * t*t));

    };

    auto Cball_dfn = [alpha,n,mean,sigma](double x)
    {
        double A = pow(n/std::abs(alpha),n) * exp(-0.5*alpha*alpha);
        double B = n/std::abs(alpha) - std::abs(alpha);
        double C = n/std::abs(alpha) * 1./( n - 1.) * exp(-0.5*alpha*alpha);
        double D = std::sqrt(DST::Math::MathCore::Pi()/2.) * (1.+std::erf(alpha/std::sqrt(2.)));
        double N = 1./sigma * 1./(C + D);

        double t = (x - mean) / sigma;
                    
        return  N* ((t <= -alpha) ? A*n/sigma*pow(B - t, -(n+1)) : -t/sigma*exp(-0.5 * t*t));

    };
    
    double dx=(maxVal - minVal)/10.;
    for(double ix = minVal; ix <= maxVal; ix += dx)
    {
        EXPECT_NEAR(p(ix),Cball_fn(ix),1e-5);
        EXPECT_NEAR(p.dfdx(ix,1e-8),Cball_dfn(ix),1e-6);
        EXPECT_NEAR(p.dfdx(ix,1e-8),p.DST::Math::pdf::dfdx(ix,1e-8),1e-6);
    }

    DST::Math::cdf_t ccdf = p.cdf(50);
    EXPECT_NEAR(ccdf.back().second,1.0,1e-6);
    EXPECT_NEAR(ccdf.front().second,0.0,1e-6);

    double norm=p.NormalisationFactor();
    EXPECT_NEAR(norm,1.0,1e-6);
    p.Normalize(100000);
    EXPECT_NEAR(p.NormalisationFactor(),norm,1e-6);

    p.dump(200);
}

TEST(Cball_distribution,random_generation)
{
    double alph = 10;
    double n = 3;
    double mean = 2.5;
    double sigma = 1.25;
    DST::Math::Cball_distribution p = DST::Math::Cball_distribution(alph,n,mean,sigma);

    double minVal = mean-10*sigma;
    double maxVal = mean+5*sigma;

    p.setRange(minVal,maxVal);

    size_t nGen = 10000;

    DST::Math::cdf_t xs = p.cdf(100);

    double xMin = xs.begin()->first; 
    double xMax = xs.rbegin()->first;
    //NORMALIZE RANGE OF CDF TO [0,1]
    for(DST::Math::cdf_t::iterator it = xs.begin(); it != xs.end(); ++it)
    {
        it->first -= xMin;
        it->first /= (xMax-xMin);
    }

    std::vector<double> histCdf;
    histCdf.reserve(xs.size());

    for(size_t i=0; i < nGen; ++i)
    {
        histCdf.push_back(cdf_lookup(xs,(p.gen()-xMin)/(xMax-xMin)));
    }
    
     // KS on Uniform(0,1) using your existing helper
     constexpr double alpha = 1e-4;
     double pval = DST::Math::ks_uniform_pvalue(histCdf, 0.0, 1.0);
     EXPECT_GT(pval, alpha);
}