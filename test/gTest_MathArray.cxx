#include <gtest/gtest.h>
#include <DSTmath/DSTarray.h>
#include <random>
#include <functional>
#include <iostream>

using namespace DST::Math;

TEST(math_array, test_double_ctor)
{

    double lower_bound = -10;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    const MaskedArray<double> ma;
    EXPECT_EQ(ma.size(), 0);
    EXPECT_EQ(ma.mask().size(), 0);
    EXPECT_EQ(ma.data().size(), 0);

    const MaskedArray<double> ma1(10);
    EXPECT_EQ(ma1.size(), 10);
    EXPECT_EQ(ma1.mask().size() , 10);
    EXPECT_EQ(ma1.data().size(), 10);

    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1.data()[i], 0,1e-10);
        EXPECT_EQ(ma1.mask()[i], false);
    }

    const MaskedArray<double> ma2(10, 1.0);
    EXPECT_EQ(ma2.size(), 10);
    EXPECT_EQ(ma2.mask().size() , 10);
    EXPECT_EQ(ma2.data().size(), 10);

    for (size_t i = 0; i < ma2.size(); i++)
    {
        EXPECT_NEAR(ma2.data()[i], 1.0,1e-10);
        EXPECT_EQ(ma2.mask()[i], false);
    }

    std::valarray<double> tdata(100);
    std::valarray<bool> tmask(100);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        double val = unif(re);
        tdata[i] = val;
        tmask[i] = (val > 0)?true:false;
    }

    const MaskedArray<double> ma3(tdata, tmask);
    EXPECT_EQ(ma3.size(), tdata.size());
    EXPECT_EQ(ma3.mask().size() , tdata.size());
    EXPECT_EQ(ma3.data().size(), tdata.size());

    for (size_t i = 0; i < ma3.size(); i++)
    {
        EXPECT_NEAR(ma3.data()[i], tdata[i],1e-10);
        EXPECT_EQ(ma3.mask()[i], tmask[i]);
    }

    const MaskedArray<double> ma4(tdata);
    EXPECT_EQ(ma4.size(), tdata.size());
    EXPECT_EQ(ma4.mask().size(), tdata.size());
    EXPECT_EQ(ma4.data().size(), tdata.size());

    for (size_t i = 0; i < ma4.size(); i++)
    {
        EXPECT_NEAR(ma4.data()[i], tdata[i],1e-10);
        EXPECT_EQ(ma4.mask()[i], false);
    }

    std::vector<double> vdata(100);
    std::vector<bool> vmask(100);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        double val = unif(re);
        vdata[i] = val;
        vmask[i] = (val <= 0)?true:false;
    }

    const MaskedArray<double> ma5(vdata, vmask);
    EXPECT_EQ(ma5.size(), vdata.size());
    EXPECT_EQ(ma5.mask().size() , vdata.size());
    EXPECT_EQ(ma5.data().size(), vdata.size());

    for (size_t i = 0; i < ma5.size(); i++)
    {
        EXPECT_NEAR(ma5.data()[i], vdata[i],1e-10);
        EXPECT_EQ(ma5.mask()[i], vmask[i]);
    }

    const MaskedArray<double> ma6(vdata);
    EXPECT_EQ(ma6.size(), vdata.size());
    EXPECT_EQ(ma6.mask().size() , vdata.size());
    EXPECT_EQ(ma6.data().size(), vdata.size());

    for (size_t i = 0; i < ma6.size(); i++)
    {
        EXPECT_NEAR(ma6.data()[i], vdata[i],1e-10);
        EXPECT_EQ(ma6.mask()[i], false);
    }

    const MaskedArray<double> ma7(ma5);
    EXPECT_EQ(ma7.size(), ma5.size());
    EXPECT_EQ(ma7.mask().size() , ma5.size());
    EXPECT_EQ(ma7.data().size(), ma5.size());

    for (size_t i = 0; i < ma7.size(); i++)
    {
        EXPECT_NEAR(ma7.data()[i], ma5.data()[i],1e-10);
        EXPECT_EQ(ma7.mask()[i], ma5.mask()[i]);
    }

    EXPECT_ANY_THROW(MaskedArray<double>(tdata, std::valarray<bool>(tdata.size()-1)));
    EXPECT_ANY_THROW(MaskedArray<double>(tdata, std::valarray<bool>(tdata.size()+1)));

    EXPECT_ANY_THROW(MaskedArray<double>(vdata, std::vector<bool>(vdata.size()-1)));
    EXPECT_ANY_THROW(MaskedArray<double>(vdata, std::vector<bool>(vdata.size()+1)));
}

TEST(math_array, test_int32_ctor)
{

    int32_t lower_bound = -10;
    int32_t upper_bound =  10;
    std::uniform_int_distribution<int32_t> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    const MaskedArray<int32_t> ma;
    EXPECT_EQ(ma.size(), 0);
    EXPECT_EQ(ma.mask().size(), 0);
    EXPECT_EQ(ma.data().size(), 0);

    const MaskedArray<int32_t> ma1(10);
    EXPECT_EQ(ma1.size(), 10);
    EXPECT_EQ(ma1.mask().size() , 10);
    EXPECT_EQ(ma1.data().size(), 10);

    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_EQ(ma1.data()[i], 0);
        EXPECT_EQ(ma1.mask()[i], false);
    }

    const MaskedArray<int32_t> ma2(10, 1.0);
    EXPECT_EQ(ma2.size(), 10);
    EXPECT_EQ(ma2.mask().size() , 10);
    EXPECT_EQ(ma2.data().size(), 10);

    for (size_t i = 0; i < ma2.size(); i++)
    {
        EXPECT_EQ(ma2.data()[i], 1.0);
        EXPECT_EQ(ma2.mask()[i], false);
    }

    std::valarray<int32_t> tdata(100);
    std::valarray<bool>    tmask(100);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        int32_t val = unif(re);
        tdata[i] = val;
        tmask[i] = (val > 0)?true:false;
    }
    
    const MaskedArray<int32_t> ma3(tdata, tmask);
    EXPECT_EQ(ma3.size(), tdata.size());
    EXPECT_EQ(ma3.mask().size() , tdata.size());
    EXPECT_EQ(ma3.data().size(), tdata.size());

    for (size_t i = 0; i < ma3.size(); i++)
    {
        EXPECT_EQ(ma3.data()[i], tdata[i]);
        EXPECT_EQ(ma3.mask()[i], tmask[i]);
    }

    const MaskedArray<int32_t> ma4(tdata);
    EXPECT_EQ(ma4.size(), tdata.size());
    EXPECT_EQ(ma4.mask().size(), tdata.size());
    EXPECT_EQ(ma4.data().size(), tdata.size());

    for (size_t i = 0; i < ma4.size(); i++)
    {
        EXPECT_EQ(ma4.data()[i], tdata[i]);
        EXPECT_EQ(ma4.mask()[i], false);
    }

    std::vector<int32_t> vdata(100);
    std::vector<bool> vmask(100);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        int32_t val = unif(re);
        vdata[i] = val;
        vmask[i] = (val <= 0)?true:false;
    }

    const MaskedArray<int32_t> ma5(vdata, vmask);
    EXPECT_EQ(ma5.size(), vdata.size());
    EXPECT_EQ(ma5.mask().size() , vdata.size());
    EXPECT_EQ(ma5.data().size(), vdata.size());

    for (size_t i = 0; i < ma5.size(); i++)
    {
        EXPECT_EQ(ma5.data()[i], vdata[i]);
        EXPECT_EQ(ma5.mask()[i], vmask[i]);
    }

    const MaskedArray<int32_t> ma6(vdata);
    EXPECT_EQ(ma6.size(), vdata.size());
    EXPECT_EQ(ma6.mask().size() , vdata.size());
    EXPECT_EQ(ma6.data().size(), vdata.size());

    for (size_t i = 0; i < ma6.size(); i++)
    {
        EXPECT_EQ(ma6.data()[i], vdata[i]);
        EXPECT_EQ(ma6.mask()[i], false);
    }

    const MaskedArray<int32_t> ma7(ma5);
    EXPECT_EQ(ma7.size(), ma5.size());
    EXPECT_EQ(ma7.mask().size() , ma5.size());
    EXPECT_EQ(ma7.data().size(), ma5.size());

    for (size_t i = 0; i < ma7.size(); i++)
    {
        EXPECT_EQ(ma7.data()[i], ma5.data()[i]);
        EXPECT_EQ(ma7.mask()[i], ma5.mask()[i]);
    }

    EXPECT_ANY_THROW(MaskedArray<int32_t>(tdata, std::valarray<bool>(tdata.size()-1)));
    EXPECT_ANY_THROW(MaskedArray<int32_t>(tdata, std::valarray<bool>(tdata.size()+1)));

    EXPECT_ANY_THROW(MaskedArray<int32_t>(vdata, std::vector<bool>(vdata.size()-1)));
    EXPECT_ANY_THROW(MaskedArray<int32_t>(vdata, std::vector<bool>(vdata.size()+1)));
}

TEST(math_array, test_uint32_ctor)
{

    uint32_t lower_bound = 0;
    uint32_t upper_bound =  100;
    std::uniform_int_distribution<int32_t> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    const MaskedArray<uint32_t> ma;
    EXPECT_EQ(ma.size(), 0);
    EXPECT_EQ(ma.mask().size(), 0);
    EXPECT_EQ(ma.data().size(), 0);

    const MaskedArray<uint32_t> ma1(10);
    EXPECT_EQ(ma1.size(), 10);
    EXPECT_EQ(ma1.mask().size() , 10);
    EXPECT_EQ(ma1.data().size(), 10);

    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_EQ(ma1.data()[i], 0);
        EXPECT_EQ(ma1.mask()[i], false);
    }

    const MaskedArray<uint32_t> ma2(10, 1.0);
    EXPECT_EQ(ma2.size(), 10);
    EXPECT_EQ(ma2.mask().size() , 10);
    EXPECT_EQ(ma2.data().size(), 10);

    for (size_t i = 0; i < ma2.size(); i++)
    {
        EXPECT_EQ(ma2.data()[i], 1.0);
        EXPECT_EQ(ma2.mask()[i], false);
    }

    std::valarray<uint32_t> tdata(100);
    std::valarray<bool> tmask(100);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        uint32_t val = unif(re);
        tdata[i] = val;
        tmask[i] = (val > 5)?true:false;
    }

    const MaskedArray<uint32_t> ma3(tdata, tmask);
    EXPECT_EQ(ma3.size(), tdata.size());
    EXPECT_EQ(ma3.mask().size() , tdata.size());
    EXPECT_EQ(ma3.data().size(), tdata.size());

    for (size_t i = 0; i < ma3.size(); i++)
    {
        EXPECT_EQ(ma3.data()[i], tdata[i]);
        EXPECT_EQ(ma3.mask()[i], tmask[i]);
    }

    const MaskedArray<uint32_t> ma4(tdata);
    EXPECT_EQ(ma4.size(), tdata.size());
    EXPECT_EQ(ma4.mask().size(), tdata.size());
    EXPECT_EQ(ma4.data().size(), tdata.size());

    for (size_t i = 0; i < ma4.size(); i++)
    {
        EXPECT_EQ(ma4.data()[i], tdata[i]);
        EXPECT_EQ(ma4.mask()[i], false);
    }

    std::vector<uint32_t> vdata(100);
    std::vector<bool> vmask(100);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        uint32_t val = unif(re);
        vdata[i] = val;
        vmask[i] = (val <= 5)?true:false;
    }

    const MaskedArray<uint32_t> ma5(vdata, vmask);
    EXPECT_EQ(ma5.size(), vdata.size());
    EXPECT_EQ(ma5.mask().size() , vdata.size());
    EXPECT_EQ(ma5.data().size(), vdata.size());

    for (size_t i = 0; i < ma5.size(); i++)
    {
        EXPECT_EQ(ma5.data()[i], vdata[i]);
        EXPECT_EQ(ma5.mask()[i], vmask[i]);
    }

    const MaskedArray<uint32_t> ma6(vdata);
    EXPECT_EQ(ma6.size(), vdata.size());
    EXPECT_EQ(ma6.mask().size() , vdata.size());
    EXPECT_EQ(ma6.data().size(), vdata.size());

    for (size_t i = 0; i < ma6.size(); i++)
    {
        EXPECT_EQ(ma6.data()[i], vdata[i]);
        EXPECT_EQ(ma6.mask()[i], false);
    }

    const MaskedArray<uint32_t> ma7(ma5);
    EXPECT_EQ(ma7.size(), ma5.size());
    EXPECT_EQ(ma7.mask().size() , ma5.size());
    EXPECT_EQ(ma7.data().size(), ma5.size());

    for (size_t i = 0; i < ma7.size(); i++)
    {
        EXPECT_EQ(ma7.data()[i], ma5.data()[i]);
        EXPECT_EQ(ma7.mask()[i], ma5.mask()[i]);
    }

    EXPECT_ANY_THROW(MaskedArray<uint32_t>(tdata, std::valarray<bool>(tdata.size()-1)));
    EXPECT_ANY_THROW(MaskedArray<uint32_t>(tdata, std::valarray<bool>(tdata.size()+1)));

    EXPECT_ANY_THROW(MaskedArray<uint32_t>(vdata, std::vector<bool>(vdata.size()-1)));
    EXPECT_ANY_THROW(MaskedArray<uint32_t>(vdata, std::vector<bool>(vdata.size()+1)));
}

TEST(math_array,test_assignement)
{

    double lower_bound = 0;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    std::valarray<double> tdata(10);
    std::valarray<bool>   tmask(10);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        double val = unif(re);
        tdata[i] = val;
        tmask[i] = (val < 0)?true:false;
    }

    const MaskedArray<double> ma(tdata, tmask);
    EXPECT_EQ(ma.size(), tdata.size());
    EXPECT_EQ(ma.mask().size() , tdata.size());
    EXPECT_EQ(ma.data().size(), tdata.size());

    for (size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma.data()[i], tdata[i],1e-10);
        EXPECT_EQ(ma.mask()[i], tmask[i]);
    }

    MaskedArray<double> ma1 = ma;
    EXPECT_EQ(ma1.size(), ma.size());
    EXPECT_EQ(ma1.mask().size() , ma.data().size());
    EXPECT_EQ(ma1.data().size(), ma.mask().size());

    for (size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma1.data()[i], ma.data()[i],1e-10);
        EXPECT_EQ(ma1.mask()[i], ma.mask()[i]);
    }

    size_t start=0;
    size_t size=tdata.size()/2;
    size_t stride=2;

    std::slice_array<double> sa = tdata[std::slice(start,size,stride)];

    ma1 = sa;
    EXPECT_EQ(ma1.size(), tdata.size()/2);

    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1.data()[i], tdata[start+i*stride],1e-10);
        EXPECT_EQ  (ma1.mask()[i], false);
    }

    start=0;
    size_t size_a=tdata.size()/3;
    size_t size_b=tdata.size()/4;
    size_t stride_a=2;
    size_t stride_b=1;

    std::gslice_array<double> gsa = tdata[std::gslice(start,{size_a,size_b},{stride_a,stride_b})];

    ma1 = gsa;
    EXPECT_EQ(ma1.size(), size_a*size_b);

    size_t n = 0;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        n = start + (i/size_b%size_a)*stride_a + (i%size_b)*stride_b;

        EXPECT_NEAR(ma1.data()[i], tdata[n],1e-10);
        EXPECT_EQ  (ma1.mask()[i], false);
    }

    std::mask_array<double> msa = tdata[tdata > 0];

    ma1 = msa;
    std::valarray<double> tdata_masked = tdata[tdata > 0];
    EXPECT_EQ(ma1.size(), tdata_masked.size());

    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1.data()[i], tdata_masked[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], false);
    }

    std::valarray<size_t> idx(tdata.size()/2);
    for (size_t i = 0; i < idx.size(); i++)
        idx[i] = i*2;
    
    std::indirect_array<double> isa= tdata[idx];
    ma1 = isa;
    EXPECT_EQ(ma1.size(), tdata.size()/2);

    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1.data()[i], tdata[i*2],1e-10);
        EXPECT_EQ  (ma1.mask()[i], false);
    }

    ma1 = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.};
    EXPECT_EQ(ma1.size(), 20);

    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1.data()[i], static_cast<double>(i+1),1e-10);
        EXPECT_EQ  (ma1.mask()[i], false);
    }
}

TEST(math_array,test_accessor)
{
    double lower_bound = -10;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    std::valarray<double> tdata(500);
    std::valarray<bool>   tmask(500);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        double val = unif(re);
        tdata[i] = val;
        tmask[i] = (val < 0)?true:false;
    }

    MaskedArray<double> ma(tdata,tmask);
    EXPECT_EQ(ma.size(), tdata.size());
    EXPECT_EQ(ma.mask().size() , tdata.size());

    size_t cnt = ma.count();
    EXPECT_EQ(cnt, ((std::valarray<bool>)tmask[!tmask]).size());

    for (size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma[i], ma.data()[i],1e-10);
        EXPECT_NEAR(ma.data()[i], tdata[i],1e-10);
        EXPECT_EQ(ma.mask()[i], tmask[i]);
    }

    ma.data() *= -2.;
    ma.mask() = !(ma.mask());

    cnt = ma.count();
    EXPECT_EQ(cnt, ((std::valarray<bool>)tmask[tmask]).size());

    for (size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma[i], ma.data()[i],1e-10);
        EXPECT_NEAR(ma.data()[i], tdata[i]*-2.,1e-10);
        EXPECT_EQ(ma.mask()[i], !tmask[i]);
    }

    ma.data() /= -2.;
    ma.mask() = !(ma.mask());

    cnt = ma.count();
    EXPECT_EQ(cnt, ((std::valarray<bool>)tmask[!tmask]).size());

    for (size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma[i], ma.data()[i],1e-10);
        EXPECT_NEAR(ma.data()[i], tdata[i],1e-10);
        EXPECT_EQ(ma.mask()[i], tmask[i]);
    }

    for (size_t i = 0; i < ma.size(); i+=2)
    {
        ma[i] *= 2.;
    }

    for (size_t i = 0; i < ma.size(); i++)
    {
        const double& val = ma[i];
        EXPECT_NEAR(ma[i], ma.data()[i],1e-10);
        EXPECT_NEAR(ma.data()[i], tdata[i]*((!(i%2))?2.:1.),1e-10);
        EXPECT_NEAR(val, tdata[i]*((!(i%2))?2.:1.),1e-10);
        EXPECT_EQ(ma.mask()[i], tmask[i]);
    }

    ma.data()*=0;
    ma.data()+=tdata;

    for (size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma.data()[i], tdata[i],1e-10);
    }

    size_t start   = 0;
    size_t size_a  = tdata.size()/4;
    size_t stride_a= 4;

    ma.data()[std::slice(start, size_a, stride_a)] = 2.;
    for (size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma[i], ((!(i%stride_a))?2.:tdata[i]),1e-10);
    }

    ma.data() = tdata;

    ma.data()[ma.mask()] = 0.;
    for (size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma[i],(tdata[i] < 0)?0:tdata[i],1e-10);
    }
}

TEST(math_array,test_modifier)
{
    double lower_bound = -10;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    std::valarray<double> tdata(500);
    std::valarray<bool>   tmask(500);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        double val = unif(re);
        tdata[i] = val;
        tmask[i] = (val < 0)?true:false;
    }

    MaskedArray<double> ma(tdata,tmask);
    
    ma.mask() &= false;
    for(size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_EQ(ma.mask()[i], false);
    }

    ma.mask() |= true;
    for(size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_EQ(ma.mask()[i], true);
    }

    ma.mask()[ma.data() >= 0] = false;
    for(size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_EQ(ma.mask()[i],(ma[i] < 0));
    }

    std::valarray<double> sdata(ma.size()/2);
    for(size_t i = 0; i < sdata.size(); i++)
    {
        double val = unif(re);
        sdata[i] = val;
    }

    ma.set_data(sdata);
    EXPECT_EQ(ma.size(), sdata.size());
    EXPECT_EQ(ma.data().size(), sdata.size());
    EXPECT_EQ(ma.mask().size(), sdata.size());

    for(size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma[i], sdata[i],1e-10);
    }

    tmask &=false;
    tmask[tdata < 0.5] = true;
    ma.set_data(tdata,tmask);

    EXPECT_EQ(ma.size(), tdata.size());
    EXPECT_EQ(ma.data().size(), tdata.size());
    EXPECT_EQ(ma.mask().size(), tdata.size());

    for(size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma[i], tdata[i],1e-10);
        EXPECT_EQ(ma.mask()[i], tdata[i]<0.5);
    }

    ma.mask() &= false;
    for(size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_EQ(ma.mask()[i], false);
    }

    tmask &=false;
    tmask[tdata > 15.] = true;
    ma.set_mask(tmask);
    EXPECT_EQ(ma.size(), tdata.size());
    EXPECT_EQ(ma.data().size(), tdata.size());
    EXPECT_EQ(ma.mask().size(), tdata.size());

    for(size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma[i], tdata[i],1e-10);
        EXPECT_EQ(ma.mask()[i], tdata[i]>15.);
    }

    for(size_t i = 0; i < ma.size(); i++)
    {
        ma.set_mask(i, tdata[i] < -0.5);
    }

    for(size_t i = 0; i < ma.size(); i++)
    {
        EXPECT_NEAR(ma[i], tdata[i],1e-10);
        EXPECT_EQ(ma.mask()[i], tdata[i]<-0.5);
    }
}

TEST(math_array,test_swap)
{
    double lower_bound = -10;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    std::valarray<double> tdata(500);
    std::valarray<bool>   tmask(500);

    std::valarray<double> sdata(50);
    std::valarray<bool>   smask(50);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        double val = unif(re);
        tdata[i] = val;
        tmask[i] = (val < 0)?true:false;
    }

    for (size_t i = 0; i < sdata.size(); i++)
    {
        double val = unif(re);
        sdata[i] = val;
        smask[i] = (val >= 0.5)?true:false;
    }

    MaskedArray<double> ma1(tdata,tmask);
    MaskedArray<double> ma2(sdata,smask);

    EXPECT_EQ(ma1.size(), tdata.size());
    EXPECT_EQ(ma2.size(), sdata.size());

    ma1.swap(ma2);
    
    EXPECT_EQ(ma1.size(), sdata.size());
    EXPECT_EQ(ma1.data().size(), sdata.size());
    EXPECT_EQ(ma1.mask().size(), smask.size());
    
    EXPECT_EQ(ma2.size(), tdata.size());
    EXPECT_EQ(ma2.data().size(), tdata.size());
    EXPECT_EQ(ma2.mask().size(), tmask.size());

    for(size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], sdata[i],1e-10);
        EXPECT_EQ(ma1.mask()[i], smask[i]);
    }

    for(size_t i = 0; i < ma2.size(); i++)
    {
        EXPECT_NEAR(ma2[i], tdata[i],1e-10);
        EXPECT_EQ(ma2.mask()[i], tmask[i]);
    }

}

TEST(math_array,test_dbl_statistics)
{
    double lower_bound = -10;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::uniform_real_distribution<double> unif2(0,1);
    std::default_random_engine re;

    std::valarray<double> tdata(5000);
    std::valarray<bool>   tmask(5000);

    std::valarray<double> sdata(5000);
    std::valarray<bool>   smask(5000);

    double min_val= std::numeric_limits<double>::max();
    double max_val=-std::numeric_limits<double>::max();
    double sum_val =0;
    double sum_val2=0;
    double cnt_val =0;
    double cnt_cov=0;
    double cov_val =0;
    size_t i_min = 0;
    size_t i_max = 0;

    for (size_t i = 0; i < tdata.size(); i++)
    {
        double val1 = unif(re);
        tdata[i] = val1;
        tmask[i] = (unif2(re) < 0.5)?true:false;

        double val2 = unif(re);
        sdata[i] = val2;
        smask[i] = (unif2(re) > 0.5)?true:false;

        if (!tmask[i])
        {
            i_min = (val1 <= min_val)?i:i_min;
            i_max = (val1 >= max_val)?i:i_max;

            min_val = (val1 <= min_val)?val1:min_val;
            max_val = (val1 >= max_val)?val1:max_val;
            sum_val += val1;
            sum_val2+= val1*val1;
            cnt_val++;

            if (!smask[i])
            {
                cov_val += val1*val2;
                cnt_cov++;
            }
        }
    }

    cov_val /= cnt_cov;
    double mean_val = sum_val/cnt_val;
    double mean_sqares_val = sum_val2/cnt_val;
    double var_val = (mean_sqares_val - mean_val*mean_val)*cnt_val/(cnt_val-1);
    
    double err_var_val = 0;
    for (size_t i = 0; i < tdata.size(); i++)
    {
        if (!tmask[i])
        {
            err_var_val += (tdata[i]-mean_val)*(tdata[i]-mean_val)*(tdata[i]-mean_val)*(tdata[i]-mean_val);
        }
    }

    err_var_val /= cnt_val;
    err_var_val = sqrt((err_var_val-(cnt_val-3)/(cnt_val-1)*var_val*var_val)/cnt_val);

    std::valarray<double> mval = tdata[!tmask];
    std::sort(std::begin(mval),std::end(mval));
    double med_val = (mval.size()%2)?mval[mval.size()/2]:(mval[mval.size()/2-1]+mval[mval.size()/2])/2;
    std::valarray<double> mad_vect(((std::valarray<double>)tdata[!tmask]).size());

    size_t j = 0;
    for (size_t i = 0; i < tdata.size(); i++)
    {
        if (!tmask[i])
        {
            mad_vect[j] = std::abs(tdata[i]-med_val);
            j++;
        }
    }

    std::sort(std::begin(mad_vect),std::end(mad_vect));
    double mad_val = (mad_vect.size()%2)?mad_vect[mad_vect.size()/2]:(mad_vect[mad_vect.size()/2-1]+mad_vect[mad_vect.size()/2])/2;



    MaskedArray<double> ma(tdata,tmask);
    MaskedArray<double> ma2(sdata,smask);

    EXPECT_NEAR(ma.sum(), sum_val,1e-10);
    EXPECT_NEAR(ma.min(), min_val,1e-10);
    EXPECT_NEAR(ma.max(), max_val,1e-10);
    EXPECT_NEAR(ma.mean(), mean_val,1e-10);
    EXPECT_NEAR(ma.sum_of_square(), sum_val2,1e-10);
    EXPECT_NEAR(ma.mean_of_square(), mean_sqares_val,1e-10);
    EXPECT_NEAR(ma.quadratic_mean(), sqrt(mean_sqares_val),1e-10);
    EXPECT_NEAR(ma.variance(), var_val,1e-10);
    EXPECT_NEAR(ma.variance_error(), err_var_val,1e-10);
    EXPECT_NEAR(ma.covariance(ma2), cov_val,1e-10);
    EXPECT_NEAR(ma.median(), med_val,1e-10);
    EXPECT_NEAR(ma.nmad(), mad_val,1e-10);
    
    EXPECT_EQ(ma.find_minimum(), i_min);
    EXPECT_EQ(ma.find_maximum(), i_max);
}

TEST(math_array,test_int_statistics)
{
    int32_t lower_bound = -10;
    int32_t upper_bound =  10;
    std::uniform_int_distribution<int32_t> unif(lower_bound,upper_bound);
    std::uniform_real_distribution<double> unif2(0,1);
    std::default_random_engine re;

    std::valarray<int32_t> tdata(500);
    std::valarray<bool>   tmask(500);

    std::valarray<int32_t> sdata(500);
    std::valarray<bool>   smask(500);

    int32_t min_val= std::numeric_limits<int32_t>::max();
    int32_t max_val=-std::numeric_limits<int32_t>::max();
    int32_t sum_val =0;
    int32_t sum_val2=0;
    int32_t cnt_val =0;
    double  cnt_cov =0;
    double  cov_val =0;
    size_t i_min    = 0;
    size_t i_max    = 0;

    for (size_t i = 0; i < tdata.size(); i++)
    {
        int32_t val1 = unif(re);
        tdata[i] = val1;
        tmask[i] = (unif2(re) < -2)?true:false;

        int32_t val2 = unif(re);
        sdata[i] = val2;
        smask[i] = (unif2(re) > 5)?true:false;

        if (!tmask[i])
        {
            i_min = (val1 < min_val)?i:i_min;
            i_max = (val1 > max_val)?i:i_max;

            min_val = (val1 <= min_val)?val1:min_val;
            max_val = (val1 >= max_val)?val1:max_val;
            sum_val += val1;
            sum_val2+= val1*val1;
            cnt_val++;

            if (!smask[i])
            {
                cov_val += val1*val2;
                cnt_cov++;
            }
        }
    }

    cov_val /= cnt_cov;
    double mean_val = static_cast<double>(sum_val)/static_cast<double>(cnt_val);
    double mean_sqares_val = static_cast<double>(sum_val2)/static_cast<double>(cnt_val);
    double var_val = (static_cast<double>(mean_sqares_val) - mean_val*mean_val)*static_cast<double>(cnt_val)/static_cast<double>(cnt_val-1);
    
    double err_var_val = 0;
    for (size_t i = 0; i < tdata.size(); i++)
    {
        if (!tmask[i])
        {
            err_var_val += (static_cast<double>(tdata[i])-mean_val)*(static_cast<double>(tdata[i])-mean_val)*(static_cast<double>(tdata[i])-mean_val)*(static_cast<double>(tdata[i])-mean_val);
        }
    }

    err_var_val /= static_cast<double>(cnt_val);
    err_var_val = sqrt((err_var_val-static_cast<double>(cnt_val-3)/static_cast<double>(cnt_val-1)*var_val*var_val)/static_cast<double>(cnt_val));

    std::valarray<int32_t> mval = tdata[!tmask];
    std::sort(std::begin(mval),std::end(mval));

    int32_t med_val = (mval.size()%2)?mval[mval.size()/2]:(mval[mval.size()/2-1]+mval[mval.size()/2])/2;
    std::valarray<int32_t> mad_vect(((std::valarray<int32_t>)tdata[!tmask]).size());

    size_t j = 0;
    for (size_t i = 0; i < tdata.size(); i++)
    {
        if (!tmask[i])
        {
            mad_vect[j] = std::abs(tdata[i]-med_val);
            j++;
        }
    }

    std::sort(std::begin(mad_vect),std::end(mad_vect));
    int32_t mad_val = (mad_vect.size()%2)?mad_vect[mad_vect.size()/2]:(mad_vect[mad_vect.size()/2-1]+mad_vect[mad_vect.size()/2])/2;



    MaskedArray<int32_t> ma(tdata,tmask);
    MaskedArray<int32_t> ma2(sdata,smask);

    EXPECT_EQ(ma.sum(), sum_val);
    EXPECT_EQ(ma.min(), min_val);
    EXPECT_EQ(ma.max(), max_val);
    EXPECT_NEAR(ma.mean(), mean_val,1e-10);
    EXPECT_NEAR(ma.sum_of_square(), sum_val2,1e-10);
    EXPECT_NEAR(ma.mean_of_square(), mean_sqares_val,1e-10);
    EXPECT_NEAR(ma.quadratic_mean(), sqrt(mean_sqares_val),1e-10);
    EXPECT_NEAR(ma.variance(), var_val,1e-10);
    EXPECT_NEAR(ma.variance_error(), err_var_val,1e-10);
    EXPECT_NEAR(ma.covariance(ma2), cov_val,1e-10);
    EXPECT_NEAR(ma.median(), med_val,1e-10);
    EXPECT_NEAR(ma.nmad(), mad_val,1e-10);
    
    EXPECT_EQ(ma.find_minimum(), i_min);
    EXPECT_EQ(ma.find_maximum(), i_max);
}

TEST(math_array,test_math_operators)
{
    double lower_bound = -10;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::uniform_real_distribution<double> unif2(0,1);
    std::default_random_engine re;

    std::valarray<double> tdata(5000);
    std::valarray<bool>   tmask(5000);

    std::valarray<double> sdata(5000);
    std::valarray<bool>   smask(5000);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        double val1 = unif(re);
        tdata[i] = val1;
        tmask[i] = (unif2(re) < -2)?true:false;

        double val2 = unif(re);
        sdata[i] = val2;
        smask[i] = (unif2(re) > 5)?true:false;
    }

    MaskedArray<double> ma1(tdata,tmask);
    MaskedArray<double> ma2(sdata,smask);

    ma1 += ma2;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]+sdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]||smask[i]);
    }
    ma1.set_data(tdata,tmask);
    
    ma1 -= ma2;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]-sdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]||smask[i]);
    }
    ma1.set_data(tdata,tmask);

    ma1*=ma2;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]*sdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]||smask[i]);
    }
    ma1.set_data(tdata,tmask);

    ma1/=ma2;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]/sdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]||smask[i]);
    }

    ma1.set_data(tdata,tmask);

    ma1+=2.;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]+2.,1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]);
    }
    ma1.set_data(tdata,tmask);

    ma1-=2.;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]-2.,1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]);
    }
    ma1.set_data(tdata,tmask);

    ma1*=2.;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]*2.,1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]);
    }
    ma1.set_data(tdata,tmask);

    ma1/=2.;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]/2.,1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]);
    }
    ma1.set_data(tdata,tmask);

    ma1+=sdata;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]+sdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]);
    }
    ma1.set_data(tdata,tmask);

    ma1-=sdata;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]-sdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]);
    }
    ma1.set_data(tdata,tmask);

    ma1*=sdata;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]*sdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]);
    }
    ma1.set_data(tdata,tmask);

    ma1/=sdata;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i]/sdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]);
    }
    ma1.set_data(tdata,tmask);
}

TEST(math_array,test_logic_transformation)
{
    double lower_bound = -10;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::uniform_real_distribution<double> unif2(0,1);
    std::default_random_engine re;

    std::valarray<double> tdata(5000);
    std::valarray<bool>   tmask(5000);

    std::valarray<double> sdata(5000);
    std::valarray<bool>   smask(5000);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        double val1 = unif(re);
        tdata[i] = val1;
        tmask[i] = (unif2(re) < -2)?true:false;

        double val2 = unif(re);
        sdata[i] = val2;
        smask[i] = (unif2(re) > 5)?true:false;
    }

    MaskedArray<double> ma1(tdata,tmask);
    MaskedArray<double> ma2(sdata,smask);

    ma1 |= ma2;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]||smask[i]);
    }
    ma1.set_mask(tmask);

    ma1 &= ma2;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]&&smask[i]);
    }
    ma1.set_mask(tmask);

    ma1 ^= ma2;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]^smask[i]);
    }
    ma1.set_mask(tmask);

    ma1 |= true;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]||true);
    }
    ma1.set_mask(tmask);

    ma1 |= false;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]||false);
    }
    ma1.set_mask(tmask);

    ma1 &= true;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]&&true);
    }
    ma1.set_mask(tmask);

    ma1 &= false;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]&&false);
    }
    ma1.set_mask(tmask);

    ma1 ^= true;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]^true);
    }
    ma1.set_mask(tmask);

    ma1 ^= false;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]^false);
    }
    ma1.set_mask(tmask);

    ma1 |= sdata;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]||smask[i]);
    }
    ma1.set_mask(tmask);

    ma1 &= sdata;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]&&smask[i]);
    }
    ma1.set_mask(tmask);

    ma1 ^= sdata;
    for (size_t i = 0; i < ma1.size(); i++)
    {
        EXPECT_NEAR(ma1[i], tdata[i],1e-10);
        EXPECT_EQ  (ma1.mask()[i], tmask[i]^smask[i]);
    }
}

TEST(math_array,test_cout)
{
    int16_t lower_bound = -10;

    std::valarray<int16_t> tdata(21);
    std::valarray<bool>   tmask(21);
    std::string str = std::string();

    for (size_t i = 0; i < tdata.size(); i++)
    {
        tdata[i] = lower_bound + i;
        tmask[i] = (i < 5)?true:false;
        str+=" "+std::to_string(tdata[i])+"["+((tmask[i])?"\033[32mtrue\033[0m":"\033[31mfalse\033[0m")+"]";
    } 

    testing::internal::CaptureStdout();

    std::cout<<MaskedArray<int16_t>(tdata,tmask);
    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_FALSE(output.compare(str)) <<"Value of: output.compare(str)"<<std::endl<<
        "  Actual: "<<output.compare(str)<<std::endl<<
        "Expected: false"<<std::endl<<
        output<<std::endl<<
        str<<std::endl;

}

TEST(math_array, test_equality_operator)
{

    int32_t lower_bound = -10;
    int32_t upper_bound =  10;
    std::uniform_int_distribution<int32_t> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    std::valarray<int32_t> tdata(500);
    std::valarray<bool>    tmask(500);

    std::valarray<int32_t> sdata(500);
    std::valarray<bool>    smask(500);

    std::valarray<int32_t> vdata(500);
    std::valarray<bool>    vmask(500);

    for (size_t i = 0; i < tdata.size(); i++)
    {
        tdata[i] = unif(re);
        tmask[i] = (tdata[i] < 0)?true:false;

        sdata[i] = unif(re);
        smask[i] = (sdata[i] < 0)?true:false;

        vdata[i] = tdata[i];
        vmask[i] = smask[i];
    }

    const MaskedArray<int32_t> ma1(tdata,tmask);
    const MaskedArray<int32_t> ma2(sdata,smask);
    const MaskedArray<int32_t> ma3(vdata,vmask);

    std::valarray<bool> test_mask1 = (ma1 == ma1);
    for(size_t i = 0; i < test_mask1.size(); i++)
    {
        EXPECT_EQ(test_mask1[i], true);
    }

    std::valarray<bool> test_mask2 = (ma1 != ma1);
    for(size_t i = 0; i < test_mask2.size(); i++)
    {
        EXPECT_EQ(test_mask2[i], false);
    }

    std::valarray<bool> test_mask3 = (ma1 == ma2);
    for(size_t i = 0; i < test_mask3.size(); i++)
    {
        EXPECT_EQ(test_mask3[i], (tdata[i] == sdata[i])&& (tmask[i] == smask[i]));
    }

    std::valarray<bool> test_mask4 = (ma1 != ma2);
    for(size_t i = 0; i < test_mask4.size(); i++)
    {
        EXPECT_EQ(test_mask4[i], (tdata[i] != sdata[i]) && (tmask[i] != smask[i]));
    }

    std::valarray<bool> test_mask5 = (ma1 == ma3);
    for(size_t i = 0; i < test_mask5.size(); i++)
    {
        EXPECT_EQ(test_mask5[i], (tdata[i] == vdata[i]) && (vmask[i] == tmask[i]));
    }

    std::valarray<bool> test_mask6 = (ma1 != ma3);
    for(size_t i = 0; i < test_mask6.size(); i++)
    {
        EXPECT_EQ(test_mask6[i], (tdata[i] != vdata[i])&&(vmask[i] != tmask[i]));
    }

    std::valarray<bool> test_mask7 = (ma1 != ma1[243]);
    for(size_t i = 0; i < test_mask7.size(); i++)
    {
        EXPECT_EQ(test_mask7[i], (ma1[i] != ma1[243])&&!tmask[i]);
    }

    std::valarray<bool> test_mask8 = (ma1 == ma1[326]);
    for(size_t i = 0; i < test_mask8.size(); i++)
    {
        EXPECT_EQ(test_mask8[i], (ma1[i] == ma1[326])&&!tmask[i]);
    }

    std::valarray<bool> test_mask9 = (ma1[243] != ma1);
    for(size_t i = 0; i < test_mask9.size(); i++)
    {
        EXPECT_EQ(test_mask9[i], (ma1[i] != ma1[243])&&!tmask[i]);
        EXPECT_EQ(test_mask7[i], test_mask9[i]);
    }

    std::valarray<bool> test_mask10 = (ma1[326] == ma1);
    for(size_t i = 0; i < test_mask10.size(); i++)
    {
        EXPECT_EQ(test_mask10[i], (ma1[i] == ma1[326])&&!tmask[i]);
        EXPECT_EQ(test_mask10[i], test_mask8[i]);
    }

}