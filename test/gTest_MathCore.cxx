#include <gtest/gtest.h>
#include <DSTmath/DSTmath.h>

using namespace DST::Math;


TEST(math_core_test, test_sign)
{

#if __cplusplus >= 199711L
    ASSERT_EQ(sgn<uint16_t>((uint16_t) 125,std::is_signed<uint16_t>())==1,true)
        <<"sgn((uint16_t)125) is positive ?!"<<sgn<uint16_t>((uint16_t) 125,std::is_signed<uint16_t>())<<std::endl;
    ASSERT_EQ(sgn<uint16_t>((uint16_t) 0,std::is_signed<uint16_t>())==0,true)
        <<"sgn((uint16_t)125) is positive ?!"<<sgn<uint16_t>((uint16_t) 0,std::is_signed<uint16_t>())<<std::endl;

    ASSERT_EQ(sgn<uint32_t>((uint32_t) 125,std::is_signed<uint32_t>())==1,true)
        <<"sgn((uint32_t)125) is positive ?!"<<sgn<uint32_t>((uint32_t) 125,std::is_signed<uint32_t>())<<std::endl;
    ASSERT_EQ(sgn<uint32_t>((uint32_t) 0,std::is_signed<uint32_t>())==0,true)
        <<"sgn((uint32_t)125) is positive ?!"<<sgn<uint32_t>((uint32_t) 0,std::is_signed<uint32_t>())<<std::endl;

    ASSERT_EQ(sgn<uint64_t>((uint64_t) 125,std::is_signed<uint64_t>())==1,true)
        <<"sgn((uint64_t)125) is positive ?!"<<sgn<uint64_t>((uint64_t) 125,std::is_signed<uint64_t>())<<std::endl;
    ASSERT_EQ(sgn<uint64_t>((uint64_t) 0,std::is_signed<uint64_t>())==0,true)
        <<"sgn((uint64_t)125) is positive ?!"<<sgn<uint64_t>((uint64_t) 0,std::is_signed<uint64_t>())<<std::endl;

    ASSERT_EQ(sgn<int16_t>((int16_t) 125,std::is_signed<int16_t>())==1,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) 125,std::is_signed<int16_t>())<<std::endl;
    ASSERT_EQ(sgn<int16_t>((int16_t) 0,std::is_signed<int16_t>())==0,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) 0,std::is_signed<int16_t>())<<std::endl;
    ASSERT_EQ(sgn<int16_t>((int16_t) -125,std::is_signed<int16_t>())==-1,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) -125,std::is_signed<int16_t>())<<std::endl;

    ASSERT_EQ(sgn<int32_t>((int32_t) 125,std::is_signed<int32_t>())==1,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) 125,std::is_signed<int32_t>())<<std::endl;
    ASSERT_EQ(sgn<int32_t>((int32_t) 0,std::is_signed<int32_t>())==0,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) 0,std::is_signed<int32_t>())<<std::endl;
    ASSERT_EQ(sgn<int32_t>((int32_t) -125,std::is_signed<int32_t>())==-1,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) -125,std::is_signed<int32_t>())<<std::endl;

    ASSERT_EQ(sgn<int64_t>((int64_t) 125,std::is_signed<int64_t>())==1,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) 125,std::is_signed<int64_t>())<<std::endl;
    ASSERT_EQ(sgn<int64_t>((int64_t) 0,std::is_signed<int64_t>())==0,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) 0,std::is_signed<int64_t>())<<std::endl;
    ASSERT_EQ(sgn<int64_t>((int64_t) -125,std::is_signed<int64_t>())==-1,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) -125,std::is_signed<int64_t>())<<std::endl;

    ASSERT_EQ(sgn<float>((float) 125,std::is_signed<float>())==1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 125,std::is_signed<float>())<<std::endl;
    ASSERT_EQ(sgn<float>((float) 0,std::is_signed<float>())==0,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 0,std::is_signed<float>())<<std::endl;
    ASSERT_EQ(sgn<float>((float) -125,std::is_signed<float>())==-1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) -125,std::is_signed<float>())<<std::endl;
    ASSERT_EQ(sgn<float>((float) 1e-6,std::is_signed<float>())==1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 1e-6,std::is_signed<float>())<<std::endl;
    ASSERT_EQ(sgn<float>((float) 0,std::is_signed<float>())==0,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 0,std::is_signed<float>())<<std::endl;
    ASSERT_EQ(sgn<float>((float) -1e-6,std::is_signed<float>())==-1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) -1e-6,std::is_signed<float>())<<std::endl;

    ASSERT_EQ(sgn<double>((double) 125,std::is_signed<double>())==1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 125,std::is_signed<double>())<<std::endl;
    ASSERT_EQ(sgn<double>((double) 0,std::is_signed<double>())==0,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 0,std::is_signed<double>())<<std::endl;
    ASSERT_EQ(sgn<double>((double) -125,std::is_signed<double>())==-1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) -125,std::is_signed<double>())<<std::endl;
    ASSERT_EQ(sgn<double>((double) 1e-6,std::is_signed<double>())==1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 1e-6,std::is_signed<double>())<<std::endl;
    ASSERT_EQ(sgn<double>((double) 0,std::is_signed<double>())==0,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 0,std::is_signed<double>())<<std::endl;
    ASSERT_EQ(sgn<double>((double) -1e-6,std::is_signed<double>())==-1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) -1e-6,std::is_signed<float>())<<std::endl;

#endif

    ASSERT_EQ(sgn<uint16_t>((uint16_t) 125)==1,true)
        <<"sgn((uint16_t)125) is positive ?!"<<sgn<uint16_t>((uint16_t) 125)<<std::endl;
    ASSERT_EQ(sgn<uint16_t>((uint16_t) 0)==0,true)
        <<"sgn((uint16_t)125) is positive ?!"<<sgn<uint16_t>((uint16_t) 0)<<std::endl;

    ASSERT_EQ(sgn<uint32_t>((uint32_t) 125)==1,true)
        <<"sgn((uint32_t)125) is positive ?!"<<sgn<uint32_t>((uint32_t) 125)<<std::endl;
    ASSERT_EQ(sgn<uint32_t>((uint32_t) 0)==0,true)
        <<"sgn((uint32_t)125) is positive ?!"<<sgn<uint32_t>((uint32_t) 0)<<std::endl;

    ASSERT_EQ(sgn<uint64_t>((uint64_t) 125)==1,true)
        <<"sgn((uint64_t)125) is positive ?!"<<sgn<uint64_t>((uint64_t) 125)<<std::endl;
    ASSERT_EQ(sgn<uint64_t>((uint64_t) 0)==0,true)
        <<"sgn((uint64_t)125) is positive ?!"<<sgn<uint64_t>((uint64_t) 0)<<std::endl;

    ASSERT_EQ(sgn<int16_t>((int16_t) 125)==1,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) 125)<<std::endl;
    ASSERT_EQ(sgn<int16_t>((int16_t) 0)==0,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) 0)<<std::endl;
    ASSERT_EQ(sgn<int16_t>((int16_t) -125)==-1,true)
        <<"sgn((int16_t)125) is positive ?!"<<sgn<int16_t>((int16_t) -125)<<std::endl;

    ASSERT_EQ(sgn<int32_t>((int32_t) 125)==1,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) 125)<<std::endl;
    ASSERT_EQ(sgn<int32_t>((int32_t) 0)==0,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) 0)<<std::endl;
    ASSERT_EQ(sgn<int32_t>((int32_t) -125)==-1,true)
        <<"sgn((int32_t)125) is positive ?!"<<sgn<int32_t>((int32_t) -125)<<std::endl;

    ASSERT_EQ(sgn<int64_t>((int64_t) 125)==1,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) 125)<<std::endl;
    ASSERT_EQ(sgn<int64_t>((int64_t) 0)==0,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) 0)<<std::endl;
    ASSERT_EQ(sgn<int64_t>((int64_t) -125)==-1,true)
        <<"sgn((int64_t)125) is positive ?!"<<sgn<int64_t>((int64_t) -125)<<std::endl;

    ASSERT_EQ(sgn<float>((float) 125)==1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 125)<<std::endl;
    ASSERT_EQ(sgn<float>((float) 0)==0,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 0)<<std::endl;
    ASSERT_EQ(sgn<float>((float) -125)==-1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) -125)<<std::endl;
    ASSERT_EQ(sgn<float>((float) 1e-6)==1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 1e-6)<<std::endl;
    ASSERT_EQ(sgn<float>((float) 0)==0,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) 0)<<std::endl;
    ASSERT_EQ(sgn<float>((float) -1e-6)==-1,true)
        <<"sgn((float)125) is positive ?!"<<sgn<float>((float) -1e-6)<<std::endl;

    ASSERT_EQ(sgn<double>((double) 125)==1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 125)<<std::endl;
    ASSERT_EQ(sgn<double>((double) 0)==0,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 0)<<std::endl;
    ASSERT_EQ(sgn<double>((double) -125)==-1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) -125)<<std::endl;
    ASSERT_EQ(sgn<double>((double) 1e-6)==1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 1e-6)<<std::endl;
    ASSERT_EQ(sgn<double>((double) 0)==0,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) 0)<<std::endl;
    ASSERT_EQ(sgn<double>((double) -1e-6)==-1,true)
        <<"sgn((double)125) is positive ?!"<<sgn<double>((double) -1e-6)<<std::endl;
}