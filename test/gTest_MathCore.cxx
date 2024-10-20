#include <gtest/gtest.h>
#include <DSTmath/DSTmath.h>

using namespace DST::Math;


TEST(math_core_test, test_sign)
{

#if __cplusplus >= 199711L
    ASSERT_EQ(sgn<uint16_t>((uint16_t) 125,std::is_signed<uint16_t>())==1,true)
        <<"sgn((uint16_t)125) is positive ?!"<<sgn<uint16_t>((uint16_t) 125,std::is_signed<uint16_t>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<uint16_t>((uint16_t) 0,std::is_signed<uint16_t>())==0,true)
        <<"sgn((uint16_t)125) is positive ?!"<<sgn<uint16_t>((uint16_t) 0,std::is_signed<uint16_t>())<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<uint32_t>((uint32_t) 125,std::is_signed<uint32_t>())==1,true)
        <<"sgn((uint32_t)125) is positive ?!"<<sgn<uint32_t>((uint32_t) 125,std::is_signed<uint32_t>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<uint32_t>((uint32_t) 0,std::is_signed<uint32_t>())==0,true)
        <<"sgn((uint32_t)125) is positive ?!"<<sgn<uint32_t>((uint32_t) 0,std::is_signed<uint32_t>())<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<uint64_t>((uint64_t) 125,std::is_signed<uint64_t>())==1,true)
        <<"sgn((uint64_t)125) is positive ?!"<<sgn<uint64_t>((uint64_t) 125,std::is_signed<uint64_t>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<uint64_t>((uint64_t) 0,std::is_signed<uint64_t>())==0,true)
        <<"sgn((uint64_t)125) is positive ?!"<<sgn<uint64_t>((uint64_t) 0,std::is_signed<uint64_t>())<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<int16_t>((int16_t) 125,std::is_signed<int16_t>())==1,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) 125,std::is_signed<int16_t>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int16_t>((int16_t) 0,std::is_signed<int16_t>())==0,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) 0,std::is_signed<int16_t>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int16_t>((int16_t) -125,std::is_signed<int16_t>())==-1,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) -125,std::is_signed<int16_t>())<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<int32_t>((int32_t) 125,std::is_signed<int32_t>())==1,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) 125,std::is_signed<int32_t>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int32_t>((int32_t) 0,std::is_signed<int32_t>())==0,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) 0,std::is_signed<int32_t>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int32_t>((int32_t) -125,std::is_signed<int32_t>())==-1,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) -125,std::is_signed<int32_t>())<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<int64_t>((int64_t) 125,std::is_signed<int64_t>())==1,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) 125,std::is_signed<int64_t>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int64_t>((int64_t) 0,std::is_signed<int64_t>())==0,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) 0,std::is_signed<int64_t>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int64_t>((int64_t) -125,std::is_signed<int64_t>())==-1,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) -125,std::is_signed<int64_t>())<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<float>((float) 125,std::is_signed<float>())==1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 125,std::is_signed<float>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<float>((float) 0,std::is_signed<float>())==0,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 0,std::is_signed<float>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<float>((float) -125,std::is_signed<float>())==-1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) -125,std::is_signed<float>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<float>((float) 1e-6,std::is_signed<float>())==1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 1e-6,std::is_signed<float>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<float>((float) 0,std::is_signed<float>())==0,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 0,std::is_signed<float>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<float>((float) -1e-6,std::is_signed<float>())==-1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) -1e-6,std::is_signed<float>())<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<double>((double) 125,std::is_signed<double>())==1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 125,std::is_signed<double>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<double>((double) 0,std::is_signed<double>())==0,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 0,std::is_signed<double>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<double>((double) -125,std::is_signed<double>())==-1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) -125,std::is_signed<double>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<double>((double) 1e-6,std::is_signed<double>())==1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 1e-6,std::is_signed<double>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<double>((double) 0,std::is_signed<double>())==0,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 0,std::is_signed<double>())<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<double>((double) -1e-6,std::is_signed<double>())==-1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) -1e-6,std::is_signed<float>())<<" ["<<__LINE__<<"]"<<std::endl;

#endif

    ASSERT_EQ(sgn<uint16_t>((uint16_t) 125)==1,true)
        <<"sgn((uint16_t)125) is positive ?!"<<sgn<uint16_t>((uint16_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<uint16_t>((uint16_t) 0)==0,true)
        <<"sgn((uint16_t)125) is positive ?!"<<sgn<uint16_t>((uint16_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((uint16_t) 125)==1,true)
        <<"sgn((uint16_t)125) is positive ?!"<<sgn((uint16_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((uint16_t) 0)==0,true)
        <<"sgn((uint16_t)125) is positive ?!"<<sgn((uint16_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<uint32_t>((uint32_t) 125)==1,true)
        <<"sgn((uint32_t)125) is positive ?!"<<sgn<uint32_t>((uint32_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<uint32_t>((uint32_t) 0)==0,true)
        <<"sgn((uint32_t)125) is positive ?!"<<sgn<uint32_t>((uint32_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((uint32_t) 125)==1,true)
        <<"sgn((uint32_t)125) is positive ?!"<<sgn((uint32_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((uint32_t) 0)==0,true)
        <<"sgn((uint32_t)125) is positive ?!"<<sgn((uint32_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<uint64_t>((uint64_t) 125)==1,true)
        <<"sgn((uint64_t)125) is positive ?!"<<sgn<uint64_t>((uint64_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<uint64_t>((uint64_t) 0)==0,true)
        <<"sgn((uint64_t)125) is positive ?!"<<sgn<uint64_t>((uint64_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((uint64_t) 125)==1,true)
        <<"sgn((uint64_t)125) is positive ?!"<<sgn((uint64_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((uint64_t) 0)==0,true)
        <<"sgn((uint64_t)125) is positive ?!"<<sgn((uint64_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<int16_t>((int16_t) 125)==1,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int16_t>((int16_t) 0)==0,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int16_t>((int16_t) -125)==-1,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) -125)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<int32_t>((int32_t) 125)==1,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int32_t>((int32_t) 0)==0,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int32_t>((int32_t) -125)==-1,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) -125)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<int64_t>((int64_t) 125)==1,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int64_t>((int64_t) 0)==0,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<int64_t>((int64_t) -125)==-1,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) -125)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<float>((float) 125)==1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<float>((float) 0)==0,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<float>((float) -125)==-1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) -125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<float>((float) 1e-6)==1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 1e-6)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<float>((float) 0)==0,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<float>((float) -1e-6)==-1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) -1e-6)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn<double>((double) 125)==1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<double>((double) 0)==0,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<double>((double) -125)==-1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) -125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<double>((double) 1e-6)==1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 1e-6)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<double>((double) 0)==0,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn<double>((double) -1e-6)==-1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) -1e-6)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn((int16_t) 125)==1,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn((int16_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((int16_t) 0)==0,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn((int16_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((int16_t) -125)==-1,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn((int16_t) -125)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn((int32_t) 125)==1,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn((int32_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((int32_t) 0)==0,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn((int32_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((int32_t) -125)==-1,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn((int32_t) -125)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn((int64_t) 125)==1,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn((int64_t) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((int64_t) 0)==0,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn((int64_t) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((int64_t) -125)==-1,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn((int64_t) -125)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn((float) 125)==1,true)
        <<"sgn((float)125) is positive ?!"<<sgn((float) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((float) 0)==0,true)
        <<"sgn((float)125) is positive ?!"<<sgn((float) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((float) -125)==-1,true)
        <<"sgn((float)125) is positive ?!"<<sgn((float) -125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((float) 1e-6)==1,true)
        <<"sgn((float)125) is positive ?!"<<sgn((float) 1e-6)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((float) 0)==0,true)
        <<"sgn((float)125) is positive ?!"<<sgn((float) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((float) -1e-6)==-1,true)
        <<"sgn((float)125) is positive ?!"<<sgn((float) -1e-6)<<" ["<<__LINE__<<"]"<<std::endl;

    ASSERT_EQ(sgn((double) 125)==1,true)
        <<"sgn((double)125) is positive ?!"<<sgn((double) 125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((double) 0)==0,true)
        <<"sgn((double)125) is positive ?!"<<sgn((double) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((double) -125)==-1,true)
        <<"sgn((double)125) is positive ?!"<<sgn((double) -125)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((double) 1e-6)==1,true)
        <<"sgn((double)125) is positive ?!"<<sgn((double) 1e-6)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((double) 0)==0,true)
        <<"sgn((double)125) is positive ?!"<<sgn((double) 0)<<" ["<<__LINE__<<"]"<<std::endl;
    ASSERT_EQ(sgn((double) -1e-6)==-1,true)
        <<"sgn((double)125) is positive ?!"<<sgn((double) -1e-6)<<" ["<<__LINE__<<"]"<<std::endl;
}

TEST(math_core_test, test_pi)
{
    ASSERT_EQ(MathCore::Pi()==std::acos(-1),true)
        <<"Pi function doesn't return "<<std::acos(-1)<<" ["<<__LINE__<<"]"<<std::endl;
}

TEST(math_core_test, gamma_function)
{
    EXPECT_DOUBLE_EQ(Func::Gamma(1),1)
        <<"Gamma(1) function doesn't return "<<1<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_DOUBLE_EQ(Func::Gamma(2),1)
        <<"Gamma(2) function doesn't return "<<1<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_DOUBLE_EQ(Func::Gamma(3),2)
        <<"Gamma(3) function doesn't return "<<2<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_DOUBLE_EQ(Func::Gamma(4),6)
        <<"Gamma(4) function doesn't return "<<6<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_DOUBLE_EQ(Func::Gamma(5),24)
        <<"Gamma(5) function doesn't return "<<24<<" ["<<__LINE__<<"]"<<std::endl;

    EXPECT_DOUBLE_EQ(Func::Gamma(1./2.), std::sqrt(MathCore::Pi()))
        <<"Gamma(1/2) function doesn't return "<<std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_DOUBLE_EQ(Func::Gamma(3./2),0.5*std::sqrt(MathCore::Pi()))
        <<"Gamma(3/2) function doesn't return "<<0.5*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_DOUBLE_EQ(Func::Gamma(5./2),3./4.*std::sqrt(MathCore::Pi()))
        <<"Gamma(3/4) function doesn't return "<<3./4.*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_DOUBLE_EQ(Func::Gamma(7./2),15./8.*std::sqrt(MathCore::Pi()))
        <<"Gamma(7/2) function doesn't return "<<15./8.*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;

    EXPECT_DOUBLE_EQ(Func::Gamma(-1./2.), -2*std::sqrt(MathCore::Pi()))
        <<"Gamma(-1/2) function doesn't return "<<-2*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_DOUBLE_EQ(Func::Gamma(-3./2.),4./3.*std::sqrt(MathCore::Pi()))
        <<"Gamma(-3/2) function doesn't return "<<4./3.*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_DOUBLE_EQ(Func::Gamma(-5./2),-8./15.*std::sqrt(MathCore::Pi()))
        <<"Gamma(-5/2) function doesn't return "<<-8./15.*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;

    EXPECT_DOUBLE_EQ(Func::Gamma(1.,std::numeric_limits<double>::max()),1)
        <<"Gamma(1,+infinity) function doesn't return "<<1<<" ["<<__LINE__<<"]"<<std::endl;

    for (int i = 0; i < 100; i++)
    for (int j = 0; j < 100; j++)
    {
        EXPECT_DOUBLE_EQ(Func::Gamma(-1*static_cast<double>(i),static_cast<double>(j)),1.)
            <<"Gamma("<<-i<<","<<j<<") function doesn't return "<<1<<" ["<<__LINE__<<"]"<<std::endl;
    }

    for (int i = 1; i < 100; i++)
    for (int j = 0; j < 100; j++)
    {
        EXPECT_DOUBLE_EQ(Func::Gamma(static_cast<double>(i),-1*static_cast<double>(j)),0.)
            <<"Gamma("<<i<<","<<-j<<") function doesn't return "<<0<<" ["<<__LINE__<<"]"<<std::endl;
    }

}