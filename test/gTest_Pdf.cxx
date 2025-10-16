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
    test_pdf() : DST::Math::pdf() {};
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

    p.Normalize(nStep);
    EXPECT_NEAR(p.getNorme(),1./(10.75+5.25),1./(2.*static_cast<double>(nStep)));
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