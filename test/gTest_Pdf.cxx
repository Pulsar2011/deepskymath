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

TEST(pdf,ctor)
{
    DST::Math::pdf p = DST::Math::pdf();

    DST::Math::MaskedArray<double> a(100);

    for(size_t i=0; i < 100; ++i)
        a.data()[i] = p.gen();
    
    EXPECT_NEAR(a.sum(),0, std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(a.mean(),0, std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(a.variance(),0, std::numeric_limits<double>::epsilon());
}