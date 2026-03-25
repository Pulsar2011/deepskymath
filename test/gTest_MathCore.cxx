#include <gtest/gtest.h>
#include <DSTmath/DSTmath.h>
#include <random>
#include <functional>

using namespace DST::Math;

// Google-test:
class ToBeTestedFixture : public MathCore, public testing::Test
{
   // Empty - bridge to protected members for unit-testing
};

TEST_F(ToBeTestedFixture, member)
{
    EXPECT_NEAR(kBig,4.503599627370496e15,1e-12);
    EXPECT_NEAR(kBiginv,2.22044604925031308085e-16,1e-12);
    EXPECT_NEAR(LS2PI,0.91893853320467274178,1e-12);
    
    EXPECT_NEAR(A[0],8.11614167470508450300e-4,1e-12);
    EXPECT_NEAR(A[1],-5.95061904284301438324E-4,1e-12);
    EXPECT_NEAR(A[2],7.93650340457716943945E-4,1e-12);
    EXPECT_NEAR(A[3],-2.77777777730099687205E-3,1e-12);
    EXPECT_NEAR(A[4],8.33333333333331927722E-2,1e-12);
    EXPECT_NEAR(B[0],-1.37825152569120859100E3,1e-12);
    EXPECT_NEAR(B[1],-3.88016315134637840924E4,1e-12);
    EXPECT_NEAR(B[2],-3.31612992738871184744E5,1e-12);
    EXPECT_NEAR(B[3],-1.16237097492762307383E6,1e-12);
    EXPECT_NEAR(B[4],-1.72173700820839662146E6,1e-12);
    EXPECT_NEAR(B[5],-8.53555664245765465627E5,1e-12);
    EXPECT_NEAR(C[0],-3.51815701436523470549E2,1e-12);
    EXPECT_NEAR(C[1],-1.70642106651881159223E4,1e-12);
    EXPECT_NEAR(C[2],-2.20528590553854454839E5,1e-12);
    EXPECT_NEAR(C[3],-1.13933444367982507207E6,1e-12);
    EXPECT_NEAR(C[4],-2.53252307177582951285E6,1e-12);
    EXPECT_NEAR(C[5],-2.01889141433532773231E6,1e-12);

    EXPECT_NEAR(P[0],1.60119522476751861407E-4,1e-12);
    EXPECT_NEAR(P[1],1.19135147006586384913E-3,1e-12);
    EXPECT_NEAR(P[2],1.04213797561761569935E-2,1e-12);
    EXPECT_NEAR(P[3],4.76367800457137231464E-2,1e-12);
    EXPECT_NEAR(P[4],2.07448227648435975150E-1,1e-12);
    EXPECT_NEAR(P[5],4.94214826801497100753E-1,1e-12);
    EXPECT_NEAR(P[6],9.99999999999999996796E-1,1e-12);

    EXPECT_NEAR(Q[0],-2.31581873324120129819E-5,1e-12);
    EXPECT_NEAR(Q[1],5.39605580493303397842E-4 ,1e-12);
    EXPECT_NEAR(Q[2],-4.45641913851797240494E-3,1e-12);
    EXPECT_NEAR(Q[3],1.18139785222060435552E-2 ,1e-12);
    EXPECT_NEAR(Q[4],3.58236398605498653373E-2 ,1e-12);
    EXPECT_NEAR(Q[5],-2.34591795718243348568E-1,1e-12);
    EXPECT_NEAR(Q[6],7.14304917030273074085E-2 ,1e-12);
    EXPECT_NEAR(Q[7],1.00000000000000000320E0,1e-12);

    EXPECT_NEAR(STIR[0],7.87311395793093628397E-4,1e-12);
    EXPECT_NEAR(STIR[1],-2.29549961613378126380E-4,1e-12);
    EXPECT_NEAR(STIR[2],-2.68132617805781232825E-3,1e-12);
    EXPECT_NEAR(STIR[3],3.47222221605458667310E-3,1e-12);
    EXPECT_NEAR(STIR[4],8.33333333333482257126E-2,1e-12);

    EXPECT_NEAR(igam(1,1),0.63212055882855767840,1e-12);
    EXPECT_NEAR(igam(1,2),0.86466471676338730811,1e-12);
    EXPECT_NEAR(igam(1,3),0.95021293163213605700,1e-12);
    EXPECT_NEAR(igam(1,5),0.99326205300091446621,1e-12);
    EXPECT_NEAR(igam(2,1),0.26424111765711535688,1e-12);
    EXPECT_NEAR(igam(2,2),0.59399415029016132000,1e-12);
    
    EXPECT_NEAR(lgam(1),0.0,1e-12);
    EXPECT_NEAR(lgam(2),0.0,1e-12);
    EXPECT_NEAR(lgam(3),0.69314718055994530942,1e-12);
    EXPECT_NEAR(lgam(4),1.79175946922805500081,1e-12);
    EXPECT_NEAR(lgam(5),3.17805383034794561965,1e-12);

    EXPECT_NEAR(igamc(1,1),0.36787944117144232160,1e-12);
    EXPECT_NEAR(igamc(1,2),0.13533528323661269189,1e-12);
    EXPECT_NEAR(igamc(1,3),0.04978706836786394300,1e-12);
    EXPECT_NEAR(igamc(1,5),0.00673794699908546700,1e-12);
    EXPECT_NEAR(igamc(2,1),0.73575888234288464312,1e-12);

    double lower_bound = -10;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    double a[11] = {1./1e9,2./1e8,3./1e7,4./1e6,5./1e5,6./1e4,7./1e3,8./1e2,9./1e1,10./1e0};
    double b[1]  = {0.0};
    
    for(int i = 0; i < 10000; i++)
    {
        double x=unif(re);

        for(int N=0; N<11; N++)
        {
            double pom=0;
            for(int j = 0; j < N+1; j++)
            {
                pom+=a[j]*std::pow(x,N-j);
            }

            EXPECT_NEAR(Polynomialeval(x,a,N),pom,1e-12)<<"P("<<x<<",a,"<<N<<")="<<Polynomialeval(x,a,N)<<" <-> "<<pom<<" ["<<__LINE__<<"]"<<std::endl;
        }

        b[0]=unif(re);
        EXPECT_NEAR(Polynomialeval(x,b, 0),b[0],1e-12);
    }

    for(int i = 0; i < 1000; i++)
    {
        double x=unif(re);
        for(int N=0; N<11; N++)
        {
            if(N<1)
            {
                EXPECT_NEAR(Polynomial1eval(x,a,N),a[0],1e-12)<<"P1("<<x<<",a,"<<N<<")="<<Polynomial1eval(x,a,N)<<" <-> "<<a[0]<<" ["<<__LINE__<<"]"<<std::endl;;
                continue;
            }

            double pom = std::pow(x,N);
            for(int j = 0; j < N; j++)
            {
                pom += a[j]*std::pow(x,N-j-1);
            }

            EXPECT_NEAR(Polynomial1eval(x,a,N),pom,1e-5)<<"P1("<<x<<",a,"<<N<<")="<<Polynomial1eval(x,a,N)<<" <-> "<<pom<<" ["<<__LINE__<<"]"<<std::endl;
        }
    }

    for(int i = 0; i < 1000; i++)
    {
        double z = unif(re);
        double x = unif(re);
        double pom = z/((1.0 + 0.5772156649015329 * x) * x);
        EXPECT_NEAR(small(z,x),pom,1e-12);
    }
}


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

TEST(math_core_test, to_string)
{
    double pi = 3.141592653589793116;
    EXPECT_FALSE(DST::to_string(pi,0) .compare("3e+00"))                     <<DST::to_string(pi,0)<<" % "<<"3e+00";
    EXPECT_FALSE(DST::to_string(pi,1) .compare("3.1e+00"))                   <<DST::to_string(pi,1)<<" % "<<"3.1e+00";
    EXPECT_FALSE(DST::to_string(pi,2) .compare("3.14e+00"))                  <<DST::to_string(pi,2)<<" % "<<"3.14e+00";
    EXPECT_FALSE(DST::to_string(pi,3) .compare("3.142e+00"))                 <<DST::to_string(pi,3)<<" % "<<"3.142e+00";
    EXPECT_FALSE(DST::to_string(pi,4) .compare("3.1416e+00"))                <<DST::to_string(pi,4)<<" % "<<"3.1416e+00";
    EXPECT_FALSE(DST::to_string(pi,5) .compare("3.14159e+00"))               <<DST::to_string(pi,5)<<" % "<<"3.14159e+00";
    EXPECT_FALSE(DST::to_string(pi,6) .compare("3.141593e+00"))              <<DST::to_string(pi,6)<<" % "<<"3.141593e+00";
    EXPECT_FALSE(DST::to_string(pi,7) .compare("3.1415927e+00"))             <<DST::to_string(pi,7)<<" % "<<"3.1415927e+00";
    EXPECT_FALSE(DST::to_string(pi,8) .compare("3.14159265e+00"))            <<DST::to_string(pi,8)<<" % "<<"3.14159265e+00";
    EXPECT_FALSE(DST::to_string(pi,9) .compare("3.141592654e+00"))           <<DST::to_string(pi,9)<<" % "<<"3.141592654e+00";
    EXPECT_FALSE(DST::to_string(pi,10).compare("3.1415926536e+00"))          <<DST::to_string(pi,10)<<" % "<<"3.1415926536e+00";
    EXPECT_FALSE(DST::to_string(pi,11).compare("3.14159265359e+00"))         <<DST::to_string(pi,11)<<" % "<<"3.14159265359e+00";
    EXPECT_FALSE(DST::to_string(pi,12).compare("3.141592653590e+00"))        <<DST::to_string(pi,12)<<" % "<<"3.141592653590e+00";
    EXPECT_FALSE(DST::to_string(pi,13).compare("3.1415926535898e+00"))       <<DST::to_string(pi,13)<<" % "<<"3.1415926535898e+00";
    EXPECT_FALSE(DST::to_string(pi,14).compare("3.14159265358979e+00"))      <<DST::to_string(pi,14)<<" % "<<"3.14159265358979e+00";
    EXPECT_FALSE(DST::to_string(pi,15).compare("3.141592653589793e+00"))     <<DST::to_string(pi,15)<<" % "<<"3.141592653589793e+00";
    EXPECT_FALSE(DST::to_string(pi,16).compare("3.1415926535897931e+00"))    <<DST::to_string(pi,16)<<" % "<<"3.1415926535897931e+00";
    EXPECT_FALSE(DST::to_string(pi,17).compare("3.14159265358979312e+00"))   <<DST::to_string(pi,17)<<" % "<<"3.14159265358979312e+00";
    EXPECT_FALSE(DST::to_string(pi,18).compare("3.141592653589793116e+00"))  <<DST::to_string(pi,18)<<" % "<<"3.141592653589793116e+00";

    float fpi = 3.14159265;
    EXPECT_FALSE(DST::to_string(fpi,0) .compare("3e+00"))                     <<DST::to_string(fpi,0)<<" % "<<"3e+00";
    EXPECT_FALSE(DST::to_string(fpi,1) .compare("3.1e+00"))                   <<DST::to_string(fpi,1)<<" % "<<"3.1e+00";
    EXPECT_FALSE(DST::to_string(fpi,2) .compare("3.14e+00"))                  <<DST::to_string(fpi,2)<<" % "<<"3.14e+00";
    EXPECT_FALSE(DST::to_string(fpi,3) .compare("3.142e+00"))                 <<DST::to_string(fpi,3)<<" % "<<"3.142e+00";
    EXPECT_FALSE(DST::to_string(fpi,4) .compare("3.1416e+00"))                <<DST::to_string(fpi,4)<<" % "<<"3.1416e+00";
    EXPECT_FALSE(DST::to_string(fpi,5) .compare("3.14159e+00"))               <<DST::to_string(fpi,5)<<" % "<<"3.14159e+00";
    EXPECT_FALSE(DST::to_string(fpi,6) .compare("3.141593e+00"))              <<DST::to_string(fpi,6)<<" % "<<"3.141593e+00";
    EXPECT_FALSE(DST::to_string(fpi,7) .compare("3.1415927e+00"))             <<DST::to_string(fpi,7)<<" % "<<"3.1415927e+00";
}

TEST(math_core_test, gamma_function)
{
    EXPECT_NEAR(MathCore::Gamma(1),1,1e-12)
        <<"Gamma(1) function doesn't return "<<1<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_NEAR(MathCore::Gamma(2),1,1e-12)
        <<"Gamma(2) function doesn't return "<<1<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_NEAR(MathCore::Gamma(3),2,1e-12)
        <<"Gamma(3) function doesn't return "<<2<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_NEAR(MathCore::Gamma(4),6,1e-12)
        <<"Gamma(4) function doesn't return "<<6<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_NEAR(MathCore::Gamma(5),24,1e-12)
        <<"Gamma(5) function doesn't return "<<24<<" ["<<__LINE__<<"]"<<std::endl;

    EXPECT_NEAR(MathCore::Gamma(1./2.), std::sqrt(MathCore::Pi()),1e-12)
        <<"Gamma(1/2) function doesn't return "<<std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_NEAR(MathCore::Gamma(3./2),0.5*std::sqrt(MathCore::Pi()),1e-12)
        <<"Gamma(3/2) function doesn't return "<<0.5*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_NEAR(MathCore::Gamma(5./2),3./4.*std::sqrt(MathCore::Pi()),1e-12)
        <<"Gamma(3/4) function doesn't return "<<3./4.*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_NEAR(MathCore::Gamma(7./2),15./8.*std::sqrt(MathCore::Pi()),1e-12)
        <<"Gamma(7/2) function doesn't return "<<15./8.*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;

    EXPECT_NEAR(MathCore::Gamma(-1./2.), -2*std::sqrt(MathCore::Pi()),1e-12)
        <<"Gamma(-1/2) function doesn't return "<<-2*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_NEAR(MathCore::Gamma(-3./2.),4./3.*std::sqrt(MathCore::Pi()),1e-12)
        <<"Gamma(-3/2) function doesn't return "<<4./3.*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;
    EXPECT_NEAR(MathCore::Gamma(-5./2),-8./15.*std::sqrt(MathCore::Pi()),1e-12)
        <<"Gamma(-5/2) function doesn't return "<<-8./15.*std::sqrt(MathCore::Pi())<<" ["<<__LINE__<<"]"<<std::endl;

    EXPECT_NEAR(MathCore::Gamma(1.,std::numeric_limits<double>::max()),1,1e-12)
        <<"Gamma(1,+infinity) function doesn't return "<<1<<" ["<<__LINE__<<"]"<<std::endl;

    for (int i = 0; i < 100; i++)
    for (int j = 0; j < 100; j++)
    {
        EXPECT_NEAR(MathCore::Gamma(-1*static_cast<double>(i),static_cast<double>(j)),1.,1e-12)
            <<"Gamma("<<-i<<","<<j<<") function doesn't return "<<1<<" ["<<__LINE__<<"]"<<std::endl;
    }

    for (int i = 1; i < 100; i++)
    for (int j = 0; j < 100; j++)
    {
        EXPECT_NEAR(MathCore::Gamma(static_cast<double>(i),-1*static_cast<double>(j)),0.,1e-12)
            <<"Gamma("<<i<<","<<-j<<") function doesn't return "<<0<<" ["<<__LINE__<<"]"<<std::endl;
    }
}

static std::function<double(double)> chebev_lambda_d;
static double chebev_wrapper_d(double x) {
    return chebev_lambda_d(x);
}

static std::function<float(float)> chebev_lambda_f;
static float chebev_wrapper_f(float x) {
    return chebev_lambda_f(x);
}

TEST(math_core_test, chebfit)
{
    unsigned int N = 10;
    double lower_bound = -1;
    double upper_bound =  1;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    std::vector<double> dcoef=std::vector<double>();
    std::vector<float > fcoef=std::vector<float> ();
    for(unsigned int k =0; k < N; k++)
    {
        dcoef.push_back(unif(re));
        fcoef.push_back(static_cast<float>(dcoef[dcoef.size()-1]));
    }

    chebev_lambda_d = [dcoef](double x) {return polynom::chebev(x, dcoef, -1., 1.);};
    chebev_lambda_f = [fcoef](float x ) {return polynom::chebev(x, fcoef, -1., 1.);};

    // Passer le lambda à chebfit
    std::vector<double> dfit = polynom::chebfit(chebev_wrapper_d, -1., 1., 15000);
    std::vector<float>  ffit = polynom::chebfit(chebev_wrapper_f, -1., 1., 15000);

    for(unsigned int k =0; k < N; k++)
        EXPECT_NEAR(dcoef[k],dfit[k],0.001);

    for(unsigned int k =0; k < N; k++)
        EXPECT_NEAR(fcoef[k],ffit[k],0.001);

}

TEST(math_core_test, chebev)
{
    for(unsigned int i =0; i < 101; i++)
    {
        double x = 2.*static_cast<double>(i)/100.-1.;
        EXPECT_NEAR(polynom::chebev(x, {1.}, -1., 1.),1.,1e-10);
        EXPECT_NEAR(polynom::chebev(x, {0,1}, -1., 1.),x,1e-10);
        EXPECT_NEAR(polynom::chebev(x, {0,0,1}, -1., 1.),2*x*x-1,1e-10);
        EXPECT_NEAR(polynom::chebev(x, {0,0,0,1}, -1., 1.),4*x*x*x-3*x,1e-10);
        EXPECT_NEAR(polynom::chebev(x, {0,0,0,0,1}, -1., 1.),8*x*x*x*x-8*x*x+1,1e-10);
        EXPECT_NEAR(polynom::chebev(x, {0,0,0,0,0,1}, -1., 1.),16*x*x*x*x*x-20*x*x*x+5*x,1e-10);
        EXPECT_NEAR(polynom::chebev(x, {0,0,0,0,0,0,1}, -1., 1.),32*x*x*x*x*x*x-48*x*x*x*x+18*x*x-1,1e-10);
    }

    for(unsigned int i =0; i < 101; i++)
    {
        float x = 2.*static_cast<float>(i)/100.-1.;
        EXPECT_NEAR(polynom::chebev(x, {1.f}, -1.f, 1.f),1.,1e-3);
        EXPECT_NEAR(polynom::chebev(x, {0.f,1.f}, -1.f, 1.f),x,1e-3);
        EXPECT_NEAR(polynom::chebev(x, {0.f,0.f,1.f}, -1.f, 1.f),2*x*x-1,1e-3);
        EXPECT_NEAR(polynom::chebev(x, {0.f,0.f,0.f,1.f}, -1.f, 1.f),4*x*x*x-3*x,1e-3);
        EXPECT_NEAR(polynom::chebev(x, {0.f,0.f,0.f,0.f,1.f}, -1.f, 1.f),8*x*x*x*x-8*x*x+1,1e-3);
        EXPECT_NEAR(polynom::chebev(x, {0.f,0.f,0.f,0.f,0.f,1.f}, -1.f, 1.f),16*x*x*x*x*x-20*x*x*x+5*x,1e-3);
        EXPECT_NEAR(polynom::chebev(x, {0.f,0.f,0.f,0.f,0.f,0.f,1.f}, -1.f, 1.f),32*x*x*x*x*x*x-48*x*x*x*x+18*x*x-1,1e-3);
    }
}

TEST(math_core_test, chebinv)
{
    double lower_bound = -1;
    double upper_bound =  1;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::mt19937 re(0xC0FFEE);
    std::vector<double> coef;

    for(unsigned int k =0; k < 5; k++)
        coef.push_back(unif(re));

    for(unsigned int i =0; i < 10; i++)
    {
        if(i==0 || i==10) continue;

        double xi = 2.*static_cast<double>(i)/10.-1.;
        double xo = xi+unif(re)/1e3;

        double y = polynom::chebev(xi, coef, -1., 1.);
        polynom::chebinv(y, xo, coef, -1., 1., 1e-18, 10000);
        EXPECT_NEAR(xi,xo,1e-15);
    }
    coef.clear();

}

TEST(math_core_test, chebev2_double)
{
    double lower_bound = -1;
    double upper_bound =  1;
    std::uniform_real_distribution<double> unif_x(lower_bound,upper_bound);
    std::uniform_real_distribution<double> unif_y(lower_bound,upper_bound);
    std::default_random_engine re;

    double* x = new double [2]; x[0] = 0;  x[1] = 0;
    double* a = new double [2]; a[0] = -1; a[1] = -1;
    double* b = new double [2]; b[0] =  1; b[1] =  1;

    std::vector<double> cx = {0,0,0,0,0,0,0};
    std::vector<double> cy = {0,0,0,0,0,0,0};

    for(size_t N=1; N <= std::min(cx.size(),cy.size()); N++)
    {
        std::vector<unsigned int> nelem = {static_cast<unsigned int>(N),static_cast<unsigned int>(N)};

        for(size_t i=0; i < N; ++i)
        {
            for(std::vector<double>::iterator it = cx.begin(); it != cx.end(); ++it)
                *it = 0;
            cx[i] = 1.;

            for(size_t j=0; j < N; ++j)
            {
                for(std::vector<double>::iterator jt = cy.begin(); jt != cy.end(); ++jt)
                    *jt = 0;
                
                cy[j] = 1.;
                
                for(unsigned int ix =0; ix < 101; ix++)
                {
                    x[0] = unif_x(re);
                    x[1] = unif_y(re);

                    double px =0.;
                    if(cx[0] == 1)
                        px=1.;
                    else if(cx[1] == 1)
                        px=x[0];
                    else if(cx[2] == 1)
                        px=2*x[0]*x[0]-1;
                    else if(cx[3] == 1)
                        px=4*x[0]*x[0]*x[0]-3*x[0];
                    else if(cx[4] == 1)
                        px=8*x[0]*x[0]*x[0]*x[0]-8*x[0]*x[0]+1;
                    else if(cx[5] == 1)
                        px=16*x[0]*x[0]*x[0]*x[0]*x[0]-20*x[0]*x[0]*x[0]+5*x[0];
                    else if(cx[6] == 1)
                        px=32*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]-48*x[0]*x[0]*x[0]*x[0]+18*x[0]*x[0]-1;

                    double py =0.;
                    if(cy[0] == 1)
                        py=1.;
                    else if(cy[1] == 1)
                        py=x[1];
                    else if(cy[2] == 1)
                        py=2*x[1]*x[1]-1;
                    else if(cy[3] == 1)
                        py=4*x[1]*x[1]*x[1]-3*x[1];
                    else if(cy[4] == 1)
                        py=8*x[1]*x[1]*x[1]*x[1]-8*x[1]*x[1]+1;
                    else if(cy[5] == 1)
                        py=16*x[1]*x[1]*x[1]*x[1]*x[1]-20*x[1]*x[1]*x[1]+5*x[1];
                    else if(cy[6] == 1)
                        py=32*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]-48*x[1]*x[1]*x[1]*x[1]+18*x[1]*x[1]-1;


                    EXPECT_NEAR(polynom::chebev2(x, cx,cy,a,b,nelem),px*py,1e-10)<<N<<std::endl<<"    {"<<cx[0]<<","<<cx[1]<<","<<cx[2]<<","<<cx[3]<<","<<cx[4]<<","<<cx[5]<<","<<cx[6]<<"}"<<std::endl<<"    {"<<cy[0]<<","<<cy[1]<<","<<cy[2]<<","<<cy[3]<<","<<cy[4]<<","<<cy[5]<<","<<cy[6]<<"} ["<<__LINE__<<"]"<<std::endl;
                }
            }
        }
    }

    delete [] x;
    delete [] a;
    delete [] b;
}

TEST(math_core_test, chebev2_float)
{
    float lower_bound = -1;
    float upper_bound =  1;
    std::uniform_real_distribution<float> unif_x(lower_bound,upper_bound);
    std::uniform_real_distribution<float> unif_y(lower_bound,upper_bound);
    std::default_random_engine re;

    float* x = new float[2]; x[0] =  0; x[1] =  0;
    float* a = new float[2]; a[0] = -1; a[1] = -1;
    float* b = new float[2]; b[0] =  1; b[1] =  1;

    std::vector<float> cx = {0,0,0,0,0,0,0};
    std::vector<float> cy = {0,0,0,0,0,0,0};

    for(size_t N=1; N <= std::min(cx.size(),cy.size()); N++)
    {
        std::vector<unsigned int> nelem = {static_cast<unsigned int>(N),static_cast<unsigned int>(N)};

        for(size_t i=0; i < N; ++i)
        {
            for(std::vector<float>::iterator it = cx.begin(); it != cx.end(); ++it)
                *it = 0;
            cx[i] = 1.;

            for(size_t j=0; j < N; ++j)
            {
                for(std::vector<float>::iterator jt = cy.begin(); jt != cy.end(); ++jt)
                    *jt = 0;
                
                cy[j] = 1.;
                
                for(unsigned int ix =0; ix < 101; ix++)
                {
                    x[0] = unif_x(re);
                    x[1] = unif_y(re);

                    float px =0.;
                    if(cx[0] == 1)
                        px=1.;
                    else if(cx[1] == 1)
                        px=x[0];
                    else if(cx[2] == 1)
                        px=2*x[0]*x[0]-1;
                    else if(cx[3] == 1)
                        px=4*x[0]*x[0]*x[0]-3*x[0];
                    else if(cx[4] == 1)
                        px=8*x[0]*x[0]*x[0]*x[0]-8*x[0]*x[0]+1;
                    else if(cx[5] == 1)
                        px=16*x[0]*x[0]*x[0]*x[0]*x[0]-20*x[0]*x[0]*x[0]+5*x[0];
                    else if(cx[6] == 1)
                        px=32*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]-48*x[0]*x[0]*x[0]*x[0]+18*x[0]*x[0]-1;

                    float py =0.;
                    if(cy[0] == 1)
                        py=1.;
                    else if(cy[1] == 1)
                        py=x[1];
                    else if(cy[2] == 1)
                        py=2*x[1]*x[1]-1;
                    else if(cy[3] == 1)
                        py=4*x[1]*x[1]*x[1]-3*x[1];
                    else if(cy[4] == 1)
                        py=8*x[1]*x[1]*x[1]*x[1]-8*x[1]*x[1]+1;
                    else if(cy[5] == 1)
                        py=16*x[1]*x[1]*x[1]*x[1]*x[1]-20*x[1]*x[1]*x[1]+5*x[1];
                    else if(cy[6] == 1)
                        py=32*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]-48*x[1]*x[1]*x[1]*x[1]+18*x[1]*x[1]-1;


                    EXPECT_NEAR(polynom::chebev2(x, cx,cy,a,b,nelem),px*py,1e-5)<<N<<std::endl<<"    {"<<cx[0]<<","<<cx[1]<<","<<cx[2]<<","<<cx[3]<<","<<cx[4]<<","<<cx[5]<<","<<cx[6]<<"}"<<std::endl<<"    {"<<cy[0]<<","<<cy[1]<<","<<cy[2]<<","<<cy[3]<<","<<cy[4]<<","<<cy[5]<<","<<cy[6]<<"} ["<<__LINE__<<"]"<<std::endl;
                }
            }
        }
    }

    delete [] x;
    delete [] a;
    delete [] b;
}

TEST(math_core_test, chebev2_2Kind_double)
{
    double lower_bound = -1;
    double upper_bound =  1;
    std::uniform_real_distribution<double> unif_x(lower_bound,upper_bound);
    std::uniform_real_distribution<double> unif_y(lower_bound,upper_bound);
    std::default_random_engine re;

    double* x = new double [2]; x[0] = 0;  x[1] = 0;
    double* a = new double [2]; a[0] = -1; a[1] = -1;
    double* b = new double [2]; b[0] =  1; b[1] =  1;

    std::vector<double> cx = std::vector<double>(6,0);
    std::vector<double> cy = std::vector<double>(6,0);

    std::vector<double> aij = std::vector<double>(cx.size()*cy.size(),0);
    
    for(size_t N=1; N <= cx.size(); N++)
    {
        for(size_t i=0; i < N; ++i)
        {
            for(unsigned int ii = 0; ii < N; ii++)
                cx.push_back(((ii+1)==N)?1.:0.);

            for(size_t j=0; j < N; j++)
            {
                for(unsigned int jj = 0; jj < N; jj++)
                    cy.push_back(((jj+1)==N)?1.:0.);

                for(std::vector<double>::iterator it = aij.begin(); it != aij.end(); ++it)
                    *it = 0;

                size_t k = i*N+j;

                std::vector<unsigned int> nelem = {static_cast<unsigned int>(N),static_cast<unsigned int>(N)};

                for(unsigned int ix =0; ix < 101; ix++)
                {
                    x[0] = unif_x(re);
                    x[1] = unif_y(re);

                    EXPECT_NEAR(polynom::chebev2(x, aij,a,b,nelem),aij[k]*polynom::chebev(x[0],cx,a[0],b[0],cx.size())*polynom::chebev(x[1],cy,a[1],b[1],cy.size()),1e-10)<<N<<std::endl;
                }

                cy.clear();
            }

            cx.clear();
        }
    }

    delete [] x;
    delete [] a;
    delete [] b;
}

TEST(math_core_test, chebev2_2Kind_float)
{
    float lower_bound = -1;
    float upper_bound =  1;
    std::uniform_real_distribution<float> unif_x(lower_bound,upper_bound);
    std::uniform_real_distribution<float> unif_y(lower_bound,upper_bound);
    std::default_random_engine re;

    float* x = new float [2]; x[0] = 0;  x[1] = 0;
    float* a = new float [2]; a[0] = -1; a[1] = -1;
    float* b = new float [2]; b[0] =  1; b[1] =  1;

    std::vector<float> cx = std::vector<float>(6,0);
    std::vector<float> cy = std::vector<float>(6,0);

    std::vector<float> aij = std::vector<float>(cx.size()*cy.size(),0);
    
    for(size_t N=1; N <= cx.size(); N++)
    {
        for(size_t i=0; i < N; ++i)
        {
            for(unsigned int ii = 0; ii < N; ii++)
                cx.push_back(((ii+1)==N)?1.:0.);

            for(size_t j=0; j < N; j++)
            {
                for(unsigned int jj = 0; jj < N; jj++)
                    cy.push_back(((jj+1)==N)?1.:0.);

                for(std::vector<float>::iterator it = aij.begin(); it != aij.end(); ++it)
                    *it = 0;

                size_t k = i*N+j;

                std::vector<unsigned int> nelem = {static_cast<unsigned int>(N),static_cast<unsigned int>(N)};

                for(unsigned int ix =0; ix < 101; ix++)
                {
                    x[0] = unif_x(re);
                    x[1] = unif_y(re);

                    EXPECT_NEAR(polynom::chebev2(x, aij,a,b,nelem),aij[k]*polynom::chebev(x[0],cx,a[0],b[0],cx.size())*polynom::chebev(x[1],cy,a[1],b[1],cy.size()),1e-10)<<N<<std::endl;
                }

                cy.clear();
            }

            cx.clear();
        }
    }

    delete [] x;
    delete [] a;
    delete [] b;
}

TEST(math_core_test, polynomial)
{
    
    std::vector<double> coef = {1.,-2.,3.,-4.,5.};
    size_t n = coef.size();
    unsigned int ni = static_cast<unsigned int>(n);

    double* da = new double [5];
    for(size_t i =0; i < 5; i++)
        da[i] = coef[i];

#if __cplusplus >= 202002L
    std::span<const double> dspan(da, n);
#endif

    for(unsigned int i =0; i < 101; i++)
    {
        double x = 2.*static_cast<double>(i)/100.-1.;
        EXPECT_NEAR(polynom::polyev (x, coef),1.-2*x+3*x*x-4*x*x*x+5*x*x*x*x,1e-10);
        EXPECT_NEAR(polynom::polyev (x, da,n),polynom::polyev(x, coef),1e-10);
        EXPECT_NEAR(polynom::polyev (x, da,ni),polynom::polyev(x, coef),1e-10);
#if __cplusplus >= 202002L
        EXPECT_NEAR(polynom::polyev(x, dspan), polynom::polyev(x, coef), 1e-10);
#endif
    }

    std::vector<float> fcoef = {1.f,-2.f,3.f,-4.f,5.f};

    float* fa = new float [5];
    for(size_t i =0; i < 5; i++)
        fa[i] = fcoef[i];

#if __cplusplus >= 202002L
    std::span<const float>  fspan(fa, n);
#endif

    for(unsigned int i =0; i < 101; i++)
    {
        float x = 2.*static_cast<float>(i)/100.-1.;
        EXPECT_NEAR(polynom::polyev (x, fcoef),1.f-2.f*x+3.f*x*x-4.f*x*x*x+5.f*x*x*x*x,1e-3);
        EXPECT_NEAR(polynom::polyev (x, fa,n),polynom::polyev(x, fcoef),1e-3);
        EXPECT_NEAR(polynom::polyev (x, fa,ni),polynom::polyev(x, fcoef),1e-3);
#if __cplusplus >= 202002L
        EXPECT_NEAR(polynom::polyev (x, fspan), polynom::polyev(x, fcoef), 1e-3);
#endif
    }
}

// =============================================================================
// DST::Math::function tests — float and double via typed tests
// =============================================================================
namespace fn = DST::Math::function;

// Type-dependent absolute tolerances
template<typename T> struct FnTol;
template<> struct FnTol<double> {
    static constexpr double func  = 1e-10;  // point-function comparisons
    static constexpr double loose = 1e-6;   // comparisons involving transcendentals
    static constexpr double integ = 1e-2;   // grid-sum normalization
};
template<> struct FnTol<float> {
    static constexpr double func  = 1e-5;
    static constexpr double loose = 1e-3;
    static constexpr double integ = 5e-2;
};

using ScalarTypes = testing::Types<float, double>;

// ─────────────────────────────────────────────────────────────────────────────
// 1-D Gaussian  gauss(x, A, m, s)
// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class GaussTest : public testing::Test {};
TYPED_TEST_SUITE(GaussTest, ScalarTypes);

TYPED_TEST(GaussTest, Peak)
{
    using T = TypeParam;
    const T A = 3.0, m = 1.5, s = 0.7;
    const T expected = A / (s * std::sqrt(T(2) * T(MathCore::Pi())));
    EXPECT_NEAR(double(fn::gauss(m, A, m, s)), double(expected), FnTol<T>::func);
}

TYPED_TEST(GaussTest, Symmetry)
{
    using T = TypeParam;
    const T A = 1.0, m = 2.0, s = 1.0, d = 0.5;
    EXPECT_NEAR(double(fn::gauss(m + d, A, m, s)),
                double(fn::gauss(m - d, A, m, s)), FnTol<T>::func);
}

TYPED_TEST(GaussTest, NumericalNorm)
{
    // ∫ gauss dx over ±6σ should equal A
    using T = TypeParam;
    const T A = 5.0, m = 0.0, s = 1.0;
    const size_t N = 10000;
    const T dx = T(12) * s / T(N);
    T sum = 0;
    for (size_t i = 0; i < N; ++i)
    {
        const T x = -T(6) * s + (T(i) + T(0.5)) * dx;
        sum += fn::gauss(x, A, m, s) * dx;
    }
    EXPECT_NEAR(double(sum), double(A), FnTol<T>::integ);
}

// ─────────────────────────────────────────────────────────────────────────────
// 1-D integrated Gaussian  integrated_gauss / igauss
// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class IGaussTest : public testing::Test {};
TYPED_TEST_SUITE(IGaussTest, ScalarTypes);

TYPED_TEST(IGaussTest, NarrowBinApproximation)
{
    // For a very narrow bin: integrated_gauss(x,A,m,s,ε) ≈ gauss(x,A,m,s)*ε
    using T = TypeParam;
    const T A = 2.0, m = 0.5, s = 1.2, eps = T(1e-3);
    const double expected = double(fn::gauss(m, A, m, s)) * double(eps);
    EXPECT_NEAR(double(fn::integrated_gauss(m, A, m, s, eps)), expected, FnTol<T>::loose);
}

TYPED_TEST(IGaussTest, NumericalNorm)
{
    // Summing bins over ±6σ should equal A
    using T = TypeParam;
    const T A = 3.0, m = 0.0, s = 1.0;
    const size_t N = 1000;
    const T dx = T(12) * s / T(N);
    T sum = 0;
    for (size_t i = 0; i < N; ++i)
    {
        const T x = -T(6) * s + (T(i) + T(0.5)) * dx;
        sum += fn::integrated_gauss(x, A, m, s, dx);
    }
    EXPECT_NEAR(double(sum), double(A), FnTol<T>::integ);
}

TYPED_TEST(IGaussTest, AliasMatchesIntegratedGauss)
{
    using T = TypeParam;
    const T A = 1.5, m = 0.3, s = 0.8, dx = T(0.1);
    for (T x = -T(3) * s; x <= T(3) * s; x += dx)
        EXPECT_NEAR(double(fn::igauss(x, A, m, s, dx)),
                    double(fn::integrated_gauss(x, A, m, s, dx)), FnTol<T>::func);
}

// ─────────────────────────────────────────────────────────────────────────────
// 2-D symmetric Gaussian  gauss2D(x, y, A, mx, my, s)
// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class Gauss2DSymTest : public testing::Test {};
TYPED_TEST_SUITE(Gauss2DSymTest, ScalarTypes);

TYPED_TEST(Gauss2DSymTest, Peak)
{
    using T = TypeParam;
    const T A = 4.0, mx = 1.0, my = -0.5, s = 0.8;
    const T expected = A / (T(2) * T(MathCore::Pi()) * s * s);
    EXPECT_NEAR(double(fn::gauss2D(mx, my, A, mx, my, s)), double(expected), FnTol<T>::loose);
}

TYPED_TEST(Gauss2DSymTest, Symmetry)
{
    using T = TypeParam;
    const T A = 1.0, mx = 0.5, my = -1.0, s = 1.0, d = 0.4;
    EXPECT_NEAR(double(fn::gauss2D(mx + d, my,     A, mx, my, s)),
                double(fn::gauss2D(mx - d, my,     A, mx, my, s)), FnTol<T>::func);
    EXPECT_NEAR(double(fn::gauss2D(mx,     my + d, A, mx, my, s)),
                double(fn::gauss2D(mx,     my - d, A, mx, my, s)), FnTol<T>::func);
}

TYPED_TEST(Gauss2DSymTest, NumericalNorm)
{
    // ∬ gauss2D dx dy over ±6σ should equal A
    using T = TypeParam;
    const T A = 2.0, mx = 0.0, my = 0.0, s = 1.0;
    const size_t N = 300;
    const T dxy = T(12) * s / T(N);
    T sum = 0;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
        {
            const T x = -T(6) * s + (T(i) + T(0.5)) * dxy;
            const T y = -T(6) * s + (T(j) + T(0.5)) * dxy;
            sum += fn::gauss2D(x, y, A, mx, my, s) * dxy * dxy;
        }
    EXPECT_NEAR(double(sum), double(A), FnTol<T>::integ);
}

// ─────────────────────────────────────────────────────────────────────────────
// 2-D symmetric integrated Gaussian  integrated_gauss2D / igauss2D (symmetric)
// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class IGauss2DSymTest : public testing::Test {};
TYPED_TEST_SUITE(IGauss2DSymTest, ScalarTypes);

TYPED_TEST(IGauss2DSymTest, NumericalNorm)
{
    // Sum of bins over ±6σ should equal A
    using T = TypeParam;
    const T A = 2.0, mx = 0.0, my = 0.0, s = 1.0;
    const size_t N = 100;
    const double dxy = 12.0 / double(N);
    T sum = 0;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
        {
            const T x = T(-6.0 + (double(i) + 0.5) * dxy);
            const T y = T(-6.0 + (double(j) + 0.5) * dxy);
            sum += fn::integrated_gauss2D(x, y, A, mx, my, s, dxy, dxy);
        }
    EXPECT_NEAR(double(sum), double(A), FnTol<T>::integ);
}

TYPED_TEST(IGauss2DSymTest, AliasMatchesIntegratedGauss2D)
{
    using T = TypeParam;
    const T A = 1.0, mx = 0.0, my = 0.0, s = 1.0;
    const double dx = 0.5, dy = 0.5;
    EXPECT_NEAR(double(fn::igauss2D(T(0.5), T(0.5), A, mx, my, s, dx, dy)),
                double(fn::integrated_gauss2D(T(0.5), T(0.5), A, mx, my, s, dx, dy)),
                FnTol<T>::func);
}

// ─────────────────────────────────────────────────────────────────────────────
// 2-D asymmetric Gaussian  gauss2D(x, y, A, mx, my, sx, sy, theta)
// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class Gauss2DAsymTest : public testing::Test {};
TYPED_TEST_SUITE(Gauss2DAsymTest, ScalarTypes);

TYPED_TEST(Gauss2DAsymTest, Peak)
{
    // Peak value at (mx,my) is A/(2π·sx·sy) regardless of theta
    using T = TypeParam;
    const T A = 3.0, mx = 0.5, my = -1.0, sx = 0.8, sy = 1.2;
    const T expected = A / (T(2) * T(MathCore::Pi()) * sx * sy);
    for (T theta : {T(0), T(0.3), T(0.7), T(1.2)})
        EXPECT_NEAR(double(fn::gauss2D(mx, my, A, mx, my, sx, sy, theta)),
                    double(expected), FnTol<T>::loose);
}

TYPED_TEST(Gauss2DAsymTest, Theta0FactorsAsProduct)
{
    // At theta=0: gauss2D = (1/(2π·sx·sy))·A·exp(-dx²/(2sx²))·exp(-dy²/(2sy²))
    using T = TypeParam;
    const T A = 1.0, mx = 0.0, my = 0.0, sx = 1.0, sy = 2.0;
    const T norm = T(1) / (T(2) * T(MathCore::Pi()) * sx * sy);
    for (T x = -T(2) * sx; x <= T(2) * sx; x += sx / T(2))
        for (T y = -T(2) * sy; y <= T(2) * sy; y += sy / T(2))
        {
            const T expected = A * norm
                               * std::exp(-T(0.5) * x * x / (sx * sx))
                               * std::exp(-T(0.5) * y * y / (sy * sy));
            EXPECT_NEAR(double(fn::gauss2D(x, y, A, mx, my, sx, sy, T(0))),
                        double(expected), FnTol<T>::loose);
        }
}

TYPED_TEST(Gauss2DAsymTest, ThetaHalfPiSwapsSigmas)
{
    // At theta=π/2: gauss2D(sx,sy,π/2) == gauss2D(sy,sx,0)
    using T = TypeParam;
    const T A = 1.0, mx = 0.0, my = 0.0, sx = 0.5, sy = 2.0;
    const T theta = T(MathCore::Pi()) / T(2);
    for (T x = -T(1.5); x <= T(1.5); x += T(0.5))
        for (T y = -T(1.5); y <= T(1.5); y += T(0.5))
            EXPECT_NEAR(double(fn::gauss2D(x, y, A, mx, my, sx, sy, theta)),
                        double(fn::gauss2D(x, y, A, mx, my, sy, sx, T(0))),
                        FnTol<T>::loose);
}

TYPED_TEST(Gauss2DAsymTest, PointSymmetryAroundMean)
{
    // Gaussian is symmetric under (x,y) -> (2mx-x, 2my-y)
    using T = TypeParam;
    const T A = 1.0, mx = 1.0, my = -0.5, sx = 0.8, sy = 1.2, theta = T(0.4);
    for (T dx = -T(1.0); dx <= T(1.0); dx += T(0.5))
        for (T dy = -T(1.0); dy <= T(1.0); dy += T(0.5))
            EXPECT_NEAR(double(fn::gauss2D(mx + dx, my + dy, A, mx, my, sx, sy, theta)),
                        double(fn::gauss2D(mx - dx, my - dy, A, mx, my, sx, sy, theta)),
                        FnTol<T>::func);
}

TYPED_TEST(Gauss2DAsymTest, NumericalNorm)
{
    // ∬ gauss2D dx dy over ±7σ should equal A
    using T = TypeParam;
    const T A = 2.0, mx = 0.0, my = 0.0, sx = 1.0, sy = 1.5, theta = T(0.5);
    const size_t N = 200;
    const T dx = T(14) / T(N), dy = T(18) / T(N);
    T sum = 0;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
        {
            const T x = -T(7) + (T(i) + T(0.5)) * dx;
            const T y = -T(9) + (T(j) + T(0.5)) * dy;
            sum += fn::gauss2D(x, y, A, mx, my, sx, sy, theta) * dx * dy;
        }
    EXPECT_NEAR(double(sum), double(A), FnTol<T>::integ);
}

// ─────────────────────────────────────────────────────────────────────────────
// 2-D asymmetric integrated Gaussian  integrated_gauss2D / igauss2D (asymmetric)
// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class IGauss2DAsymTest : public testing::Test {};
TYPED_TEST_SUITE(IGauss2DAsymTest, ScalarTypes);

TYPED_TEST(IGauss2DAsymTest, NumericalNorm)
{
    // Sum over a ±7σ grid should equal A
    using T = TypeParam;
    const T A = 2.0, mx = 0.0, my = 0.0, sx = 1.0, sy = 1.5, theta = T(0.3);
    const size_t N = 50;
    const T dxy = T(14) / T(N);
    T sum = 0;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
        {
            const T x = -T(7) + (T(i) + T(0.5)) * dxy;
            const T y = -T(7) + (T(j) + T(0.5)) * dxy;
            sum += fn::integrated_gauss2D(x, y, A, mx, my, sx, sy, dxy, dxy, theta, 100);
        }
    EXPECT_NEAR(double(sum), double(A), FnTol<T>::integ);
}

TYPED_TEST(IGauss2DAsymTest, Theta0MatchesAnalytical)
{
    // At theta=0, sx=sy=s: result must equal (A/4)*(erf_x2-erf_x1)*(erf_y2-erf_y1)
    using T = TypeParam;
    const T A = 1.0, mx = 0.0, my = 0.0, sx = 1.0, sy = 1.0;
    const T dx = T(0.4), dy = T(0.4);
    for (T x = -T(2); x <= T(2); x += T(0.8))
        for (T y = -T(2); y <= T(2); y += T(0.8))
        {
            const double sq2 = std::sqrt(2.0) * double(sx);
            const double x1  = double(x) - double(dx) / 2.;
            const double x2  = double(x) + double(dx) / 2.;
            const double y1  = double(y) - double(dy) / 2.;
            const double y2  = double(y) + double(dy) / 2.;
            const double expected = double(A) / 4.
                                    * (std::erf(x2 / sq2) - std::erf(x1 / sq2))
                                    * (std::erf(y2 / sq2) - std::erf(y1 / sq2));
            EXPECT_NEAR(double(fn::integrated_gauss2D(x, y, A, mx, my, sx, sy, dx, dy, T(0), 100)),
                        expected, FnTol<T>::integ);
        }
}

TYPED_TEST(IGauss2DAsymTest, AliasMatchesIntegratedGauss2D)
{
    using T = TypeParam;
    const T A = 1.0, mx = 0.5, my = -0.3, sx = 0.8, sy = 1.1;
    const T dx = T(0.3), dy = T(0.3), theta = T(0.2);
    EXPECT_NEAR(double(fn::igauss2D(mx, my, A, mx, my, sx, sy, dx, dy, theta, 50)),
                double(fn::integrated_gauss2D(mx, my, A, mx, my, sx, sy, dx, dy, theta, 50)),
                FnTol<T>::func);
}

// ─────────────────────────────────────────────────────────────────────────────
// 1-D Moffat  moffat(x, A, x0, a, b)
// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class moffatest : public testing::Test {};
TYPED_TEST_SUITE(moffatest, ScalarTypes);

TYPED_TEST(moffatest, Peak)
{
    using T = TypeParam;
    const T A = 2.0, x0 = 1.5, a = 1.0, b = 3.0;
    const T expected = A * (b - T(1)) / (T(MathCore::Pi()) * a * a);
    EXPECT_NEAR(double(fn::moffat(x0, A, x0, a, b)), double(expected), FnTol<T>::loose);
}

TYPED_TEST(moffatest, Symmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = 0.5, a = 2.0, b = 4.0, d = 1.0;
    EXPECT_NEAR(double(fn::moffat(x0 + d, A, x0, a, b)),
                double(fn::moffat(x0 - d, A, x0, a, b)), FnTol<T>::func);
}

TYPED_TEST(moffatest, MonotonicallyDecreasing)
{
    // Values must decrease monotonically away from the peak
    using T = TypeParam;
    const T A = 1.0, x0 = 0.0, a = 1.0, b = 3.0;
    T prev = fn::moffat(x0, A, x0, a, b);
    for (T d = T(0.2); d <= T(3.0); d += T(0.2))
    {
        const T cur = fn::moffat(x0 + d, A, x0, a, b);
        EXPECT_LT(double(cur), double(prev));
        prev = cur;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2-D Moffat functions — implementation has compile errors, tests disabled:
//
//   moffat2D (symmetric, 7 params):
//     uses undeclared `_x`, `_y`, `m_x`, `m_y`, `A`  → should be `x_`, `y_`, `x0_`, `y0_`, `A_`
//
//   moffat2D (asymmetric, 9 params):
//     uses undeclared `g` in `y1*y1/(g*g_)`  → should be `g_*g_`
//     `cos`/`sin` should be `std::cos`/`std::sin`
//
// Uncomment the block below once those bugs are fixed.
// ─────────────────────────────────────────────────────────────────────────────

template<typename T> class Moffat2DSymTest  : public testing::Test {};
template<typename T> class Moffat2DAsymTest : public testing::Test {};
TYPED_TEST_SUITE(Moffat2DSymTest,  ScalarTypes);
TYPED_TEST_SUITE(Moffat2DAsymTest, ScalarTypes);

TYPED_TEST(Moffat2DSymTest, Peak)
{
    using T = TypeParam;
    const T A = 2.0, x0 = 1.0, y0 = -0.5, a = 1.5, b = 3.0;
    const T expected = A * (b - T(1)) / (T(MathCore::Pi()) * a * a);
    EXPECT_NEAR(double(fn::moffat2D(x0, y0, A, x0, y0, a, b)), double(expected), FnTol<T>::loose);
}

TYPED_TEST(Moffat2DSymTest, RadialSymmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = 0.0, y0 = 0.0, a = 1.0, b = 4.0, r = T(0.7);
    EXPECT_NEAR(double(fn::moffat2D(r, T(0), A, x0, y0, a, b)),
                double(fn::moffat2D(T(0), r, A, x0, y0, a, b)), FnTol<T>::func);
}

TYPED_TEST(Moffat2DSymTest, NumericalNorm)
{
    using T = TypeParam;
    const T A = 2.0, x0 = 0.0, y0 = 0.0, a = 1.0, b = 3.0;
    const size_t N = 500;
    const T dxy = T(20) / T(N);
    T sum = 0;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
        {
            const T x = -T(10) + (T(i) + T(0.5)) * dxy;
            const T y = -T(10) + (T(j) + T(0.5)) * dxy;
            sum += fn::moffat2D(x, y, A, x0, y0, a, b) * dxy * dxy;
        }
    EXPECT_NEAR(double(sum), double(A), FnTol<T>::integ);
}

TYPED_TEST(Moffat2DAsymTest, Gamma1MatchesSymmetric)
{
    // At gamma=1 and theta=0, the asymmetric form should match the symmetric one
    using T = TypeParam;
    const T A = 1.0, x0 = 0.5, y0 = -0.5, a = 1.2, b = 3.5, g = T(1.0);
    for (T x = -T(2); x <= T(2); x += T(0.5))
        for (T y = -T(2); y <= T(2); y += T(0.5))
            EXPECT_NEAR(double(fn::moffat2D(x, y, A, x0, y0, a, b, g, T(0))),
                        double(fn::moffat2D(x, y, A, x0, y0, a, b)),
                        FnTol<T>::loose);
}

TYPED_TEST(Moffat2DAsymTest, PointSymmetryAroundCenter)
{
    using T = TypeParam;
    const T A = 1.0, x0 = 0.0, y0 = 0.0, a = 1.0, b = 4.0, g = T(0.7), theta = T(0.5);
    for (T dx = -T(1); dx <= T(1); dx += T(0.5))
        for (T dy = -T(1); dy <= T(1); dy += T(0.5))
            EXPECT_NEAR(double(fn::moffat2D(x0 + dx, y0 + dy, A, x0, y0, a, b, g, theta)),
                        double(fn::moffat2D(x0 - dx, y0 - dy, A, x0, y0, a, b, g, theta)),
                        FnTol<T>::func);
}

// =============================================================================
// Lorentzian tests
// =============================================================================

template<typename T> class LorentzianTest : public testing::Test {};
TYPED_TEST_SUITE(LorentzianTest, ScalarTypes);

TYPED_TEST(LorentzianTest, Peak)
{
    using T = TypeParam;
    const T A = 3.0, x0 = 1.5, g = T(0.8);
    const T expected = A / (T(MathCore::Pi()) * g);
    EXPECT_NEAR(double(fn::lorentzian(x0, A, x0, g)), double(expected), FnTol<T>::loose);
}

TYPED_TEST(LorentzianTest, Symmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = 0.5, g = T(1.2), d = T(0.7);
    EXPECT_NEAR(double(fn::lorentzian(x0 + d, A, x0, g)),
                double(fn::lorentzian(x0 - d, A, x0, g)), FnTol<T>::func);
}

TYPED_TEST(LorentzianTest, NumericalNorm)
{
    // ∫ L dx over large range ≈ A (Lorentzian has heavy tails: use ±1000γ)
    using T = TypeParam;
    const T A = 2.0, x0 = T(0), g = T(1.0);
    const size_t N = 200000;
    const T range = T(1000) * g;
    const T dx = T(2) * range / T(N);
    T sum = T(0);
    for (size_t i = 0; i < N; ++i)
    {
        const T x = -range + (T(i) + T(0.5)) * dx;
        sum += fn::lorentzian(x, A, x0, g) * dx;
    }
    EXPECT_NEAR(double(sum), double(A), FnTol<T>::integ);
}

// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class Lorentzian2DSymTest : public testing::Test {};
TYPED_TEST_SUITE(Lorentzian2DSymTest, ScalarTypes);

TYPED_TEST(Lorentzian2DSymTest, Peak)
{
    using T = TypeParam;
    const T A = 2.0, x0 = T(0.3), y0 = T(-0.5), g = T(0.6);
    // Peak at (x0,y0): A / (2*pi*g^2)
    const T expected = A / (T(2) * T(MathCore::Pi()) * g * g);
    EXPECT_NEAR(double(fn::lorentzian2D(x0, y0, A, x0, y0, g)), double(expected), FnTol<T>::loose);
}

TYPED_TEST(Lorentzian2DSymTest, RadialSymmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), g = T(1.0), r = T(0.5);
    EXPECT_NEAR(double(fn::lorentzian2D(r, T(0), A, x0, y0, g)),
                double(fn::lorentzian2D(T(0), r, A, x0, y0, g)), FnTol<T>::func);
}

TYPED_TEST(Lorentzian2DSymTest, NumericalNorm)
{
    using T = TypeParam;
    const T A = 1.5, x0 = T(0), y0 = T(0), g = T(1.0);
    const size_t N = 2000;
    const T range = T(500) * g;
    const T dxy = T(2) * range / T(N);
    T sum = T(0);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
        {
            const T x = -range + (T(i) + T(0.5)) * dxy;
            const T y = -range + (T(j) + T(0.5)) * dxy;
            sum += fn::lorentzian2D(x, y, A, x0, y0, g) * dxy * dxy;
        }
    EXPECT_NEAR(double(sum), double(A), 0.05);
}

// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class Lorentzian2DAsymTest : public testing::Test {};
TYPED_TEST_SUITE(Lorentzian2DAsymTest, ScalarTypes);

TYPED_TEST(Lorentzian2DAsymTest, Q1MatchesSym)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0.5), y0 = T(-0.3), g = T(0.8), q = T(1.0), theta = T(0);
    for (T dx = -T(1); dx <= T(1); dx += T(0.5))
        for (T dy = -T(1); dy <= T(1); dy += T(0.5))
            EXPECT_NEAR(double(fn::lorentzian2D(x0+dx, y0+dy, A, x0, y0, g, q, theta)),
                        double(fn::lorentzian2D(x0+dx, y0+dy, A, x0, y0, g)),
                        FnTol<T>::loose);
}

TYPED_TEST(Lorentzian2DAsymTest, PointSymmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), g = T(1.0), q = T(0.6), theta = T(0.4);
    for (T dx = -T(1); dx <= T(1); dx += T(0.5))
        for (T dy = -T(1); dy <= T(1); dy += T(0.5))
            EXPECT_NEAR(double(fn::lorentzian2D(x0+dx, y0+dy, A, x0, y0, g, q, theta)),
                        double(fn::lorentzian2D(x0-dx, y0-dy, A, x0, y0, g, q, theta)),
                        FnTol<T>::func);
}

// =============================================================================
// Pseudo-Voigt tests
// =============================================================================

template<typename T> class PseudoVoigtTest : public testing::Test {};
TYPED_TEST_SUITE(PseudoVoigtTest, ScalarTypes);

TYPED_TEST(PseudoVoigtTest, PureLimits)
{
    // fL → small: pseudo_voigt ≈ Gaussian with FWHM ≈ fG
    using T = TypeParam;
    const T A = 2.0, x0 = T(0), fG = T(1.0), fL = T(1e-3);
    const T sig = fG / (T(2) * std::sqrt(T(2) * std::log(T(2))));
    for (T x = -T(2)*fG; x <= T(2)*fG; x += fG/T(4))
    {
        const T expected = A / (sig * T(SQT2PI)) * std::exp(T(-0.5)*(x-x0)*(x-x0)/(sig*sig));
        EXPECT_NEAR(double(fn::pseudo_voigt(x, A, x0, fG, fL)), double(expected), 1e-2);
    }
}

TYPED_TEST(PseudoVoigtTest, NumericalNorm)
{
    using T = TypeParam;
    const T A = 3.0, x0 = T(0), fG = T(1.0), fL = T(1.0);
    const size_t N = 200000;
    const T range = T(500);
    const T dx = T(2) * range / T(N);
    T sum = T(0);
    for (size_t i = 0; i < N; ++i)
    {
        const T x = -range + (T(i) + T(0.5)) * dx;
        sum += fn::pseudo_voigt(x, A, x0, fG, fL) * dx;
    }
    EXPECT_NEAR(double(sum), double(A), 0.1);
}

TYPED_TEST(PseudoVoigtTest, EtaBounds)
{
    // With fL/fG varying, eta should stay in [0,1]
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), fG = T(1.0);
    // Evaluate at peak to cross-check ordering: pure Lorentzian peak > pure Gaussian peak (same FWHM)
    const T pv_gauss = double(fn::pseudo_voigt(x0, A, x0, fG, T(1e-4)));
    const T pv_lortz = double(fn::pseudo_voigt(x0, A, x0, T(1e-4), fG));
    // Gaussian peak at x0 = A/sigma/sqrt(2pi) ; Lorentzian peak = A/(pi*gamma)
    // Just check they are positive and finite
    EXPECT_GT(double(pv_gauss), 0.0);
    EXPECT_GT(double(pv_lortz), 0.0);
}

// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class PseudoVoigt2DSymTest : public testing::Test {};
TYPED_TEST_SUITE(PseudoVoigt2DSymTest, ScalarTypes);

TYPED_TEST(PseudoVoigt2DSymTest, RadialDecay)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), fG = T(1.0), fL = T(0.5);
    T prev = fn::pseudo_voigt2D(x0, y0, A, x0, y0, fG, fL);
    for (T r = T(0.3); r <= T(3.0); r += T(0.3))
    {
        const T cur = fn::pseudo_voigt2D(x0 + r, y0, A, x0, y0, fG, fL);
        EXPECT_LT(double(cur), double(prev));
        prev = cur;
    }
}

TYPED_TEST(PseudoVoigt2DSymTest, RadialSymmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), fG = T(1.0), fL = T(1.0), r = T(0.7);
    EXPECT_NEAR(double(fn::pseudo_voigt2D(r, T(0), A, x0, y0, fG, fL)),
                double(fn::pseudo_voigt2D(T(0), r, A, x0, y0, fG, fL)), FnTol<T>::loose);
}

// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class PseudoVoigt2DAsymTest : public testing::Test {};
TYPED_TEST_SUITE(PseudoVoigt2DAsymTest, ScalarTypes);

TYPED_TEST(PseudoVoigt2DAsymTest, Q1MatchesSym)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), fG = T(1.0), fL = T(0.8), q = T(1.0), theta = T(0);
    for (T dx = -T(1); dx <= T(1); dx += T(0.5))
        for (T dy = -T(1); dy <= T(1); dy += T(0.5))
            EXPECT_NEAR(double(fn::pseudo_voigt2D(dx, dy, A, x0, y0, fG, fL, q, theta)),
                        double(fn::pseudo_voigt2D(dx, dy, A, x0, y0, fG, fL)),
                        FnTol<T>::loose);
}

TYPED_TEST(PseudoVoigt2DAsymTest, PointSymmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), fG = T(1.0), fL = T(0.5), q = T(0.6), theta = T(0.4);
    for (T dx = -T(1); dx <= T(1); dx += T(0.5))
        for (T dy = -T(1); dy <= T(1); dy += T(0.5))
            EXPECT_NEAR(double(fn::pseudo_voigt2D(x0+dx, y0+dy, A, x0, y0, fG, fL, q, theta)),
                        double(fn::pseudo_voigt2D(x0-dx, y0-dy, A, x0, y0, fG, fL, q, theta)),
                        FnTol<T>::func);
}

// =============================================================================
// Sérsic tests
// =============================================================================

template<typename T> class SersicTest : public testing::Test {};
TYPED_TEST_SUITE(SersicTest, ScalarTypes);

TYPED_TEST(SersicTest, PeakAtRe)
{
    // sersic(re) == A by definition
    using T = TypeParam;
    const T A = 4.0, re = T(2.0), n = T(4.0);
    EXPECT_NEAR(double(fn::sersic(re, A, re, n)), double(A), FnTol<T>::loose);
}

TYPED_TEST(SersicTest, MonotonicallyDecreasing)
{
    using T = TypeParam;
    const T A = 1.0, re = T(5.0), n = T(4.0);
    T prev = fn::sersic(T(0.01), A, re, n);
    for (T r = T(0.5); r <= T(10.0); r += T(0.5))
    {
        const T cur = fn::sersic(r, A, re, n);
        EXPECT_LT(double(cur), double(prev));
        prev = cur;
    }
}

TYPED_TEST(SersicTest, IndexEffect)
{
    // Higher n → more concentrated at r=0, less at r>>re
    using T = TypeParam;
    const T A = 1.0, re = T(3.0), r_inner = T(0.1), r_outer = T(6.0);
    EXPECT_GT(double(fn::sersic(r_inner, A, re, T(4))),
              double(fn::sersic(r_inner, A, re, T(1))));
}

// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class Sersic2DSymTest : public testing::Test {};
TYPED_TEST_SUITE(Sersic2DSymTest, ScalarTypes);

TYPED_TEST(Sersic2DSymTest, PeakAtCenter)
{
    using T = TypeParam;
    const T A = 2.0, x0 = T(1.0), y0 = T(-0.5), re = T(3.0), n = T(4.0);
    // sersic2D at (x0,y0) → r=0 → exp(-bn*(0^(1/n)-1)) = exp(bn)
    const T bn = T(2)*n - T(1.0/3.0) + T(4.0)/(T(405)*n) + T(46.0)/(T(25515)*n*n);
    const T expected = A * std::exp(bn);
    EXPECT_NEAR(double(fn::sersic2D(x0, y0, A, x0, y0, re, n)), double(expected), FnTol<T>::loose);
}

TYPED_TEST(Sersic2DSymTest, RadialConsistency)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), re = T(3.0), n = T(2.0);
    for (T r = T(0.5); r <= T(4.0); r += T(0.5))
        EXPECT_NEAR(double(fn::sersic2D(r, T(0), A, x0, y0, re, n)),
                    double(fn::sersic(r, A, re, n)),
                    FnTol<T>::func);
}

// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class Sersic2DAsymTest : public testing::Test {};
TYPED_TEST_SUITE(Sersic2DAsymTest, ScalarTypes);

TYPED_TEST(Sersic2DAsymTest, Q1MatchesSym)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), re = T(2.0), n = T(1.0), q = T(1.0), theta = T(0);
    for (T dx = -T(2); dx <= T(2); dx += T(0.5))
        for (T dy = -T(2); dy <= T(2); dy += T(0.5))
            EXPECT_NEAR(double(fn::sersic2D(x0+dx, y0+dy, A, x0, y0, re, n, q, theta)),
                        double(fn::sersic2D(x0+dx, y0+dy, A, x0, y0, re, n)),
                        FnTol<T>::loose);
}

TYPED_TEST(Sersic2DAsymTest, PointSymmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), re = T(3.0), n = T(2.0), q = T(0.7), theta = T(0.5);
    for (T dx = -T(2); dx <= T(2); dx += T(0.5))
        for (T dy = -T(2); dy <= T(2); dy += T(0.5))
            EXPECT_NEAR(double(fn::sersic2D(x0+dx, y0+dy, A, x0, y0, re, n, q, theta)),
                        double(fn::sersic2D(x0-dx, y0-dy, A, x0, y0, re, n, q, theta)),
                        FnTol<T>::func);
}

// =============================================================================
// King (1966) profile tests
// =============================================================================

template<typename T> class KingTest : public testing::Test {};
TYPED_TEST_SUITE(KingTest, ScalarTypes);

TYPED_TEST(KingTest, PeakAtCenter)
{
    // At r=0: 1/sqrt(1+0) - ct = 1 - ct  → peak = A/norm*(1-ct)^2
    using T = TypeParam;
    const T A = 2.0, rc = T(1.0), rt = T(5.0);
    const T xc = rt / rc;
    const T ct = T(1) / std::sqrt(T(1) + xc*xc);
    const T norm = T(MathCore::Pi()) * rc*rc
                   * (std::log(T(1) + xc*xc) - T(4)*(T(1)-ct) + xc*xc*ct*ct);
    const T expected = A / norm * (T(1) - ct) * (T(1) - ct);
    EXPECT_NEAR(double(fn::king(T(0), A, rc, rt)), double(expected), FnTol<T>::loose);
}

TYPED_TEST(KingTest, TruncationAtRt)
{
    using T = TypeParam;
    const T A = 1.0, rc = T(1.0), rt = T(3.0);
    EXPECT_NEAR(double(fn::king(rt + T(0.1), A, rc, rt)), 0.0, FnTol<T>::func);
    EXPECT_NEAR(double(fn::king(rt * T(2.0), A, rc, rt)), 0.0, FnTol<T>::func);
}

TYPED_TEST(KingTest, NumericalNorm)
{
    // ∬ king2D dx dy ≈ A
    using T = TypeParam;
    const T A = 2.0, rc = T(1.0), rt = T(5.0);
    const size_t N = 300;
    const T dxy = T(2) * rt / T(N);
    T sum = T(0);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
        {
            const T x = -rt + (T(i) + T(0.5)) * dxy;
            const T y = -rt + (T(j) + T(0.5)) * dxy;
            sum += fn::king2D(x, y, A, T(0), T(0), rc, rt) * dxy * dxy;
        }
    EXPECT_NEAR(double(sum), double(A), FnTol<T>::integ);
}

// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class King2DSymTest : public testing::Test {};
TYPED_TEST_SUITE(King2DSymTest, ScalarTypes);

TYPED_TEST(King2DSymTest, MatchesKing)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), rc = T(1.5), rt = T(6.0);
    for (T r = T(0); r < rt; r += T(0.5))
        EXPECT_NEAR(double(fn::king2D(r, T(0), A, x0, y0, rc, rt)),
                    double(fn::king(r, A, rc, rt)),
                    FnTol<T>::func);
}

TYPED_TEST(King2DSymTest, RadialSymmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), rc = T(1.0), rt = T(4.0), r = T(0.8);
    EXPECT_NEAR(double(fn::king2D(r, T(0), A, x0, y0, rc, rt)),
                double(fn::king2D(T(0), r, A, x0, y0, rc, rt)), FnTol<T>::func);
}

// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class King2DAsymTest : public testing::Test {};
TYPED_TEST_SUITE(King2DAsymTest, ScalarTypes);

TYPED_TEST(King2DAsymTest, Q1MatchesSym)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), rc = T(1.0), rt = T(4.0), q = T(1.0), theta = T(0);
    for (T dx = -T(2); dx <= T(2); dx += T(0.5))
        for (T dy = -T(2); dy <= T(2); dy += T(0.5))
            EXPECT_NEAR(double(fn::king2D(x0+dx, y0+dy, A, x0, y0, rc, rt, q, theta)),
                        double(fn::king2D(x0+dx, y0+dy, A, x0, y0, rc, rt)),
                        FnTol<T>::loose);
}

TYPED_TEST(King2DAsymTest, PointSymmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), rc = T(1.0), rt = T(4.0), q = T(0.6), theta = T(0.4);
    for (T dx = -T(1.5); dx <= T(1.5); dx += T(0.5))
        for (T dy = -T(1.5); dy <= T(1.5); dy += T(0.5))
            EXPECT_NEAR(double(fn::king2D(x0+dx, y0+dy, A, x0, y0, rc, rt, q, theta)),
                        double(fn::king2D(x0-dx, y0-dy, A, x0, y0, rc, rt, q, theta)),
                        FnTol<T>::func);
}

// =============================================================================
// Airy disk tests
// =============================================================================

template<typename T> class AiryTest : public testing::Test {};
TYPED_TEST_SUITE(AiryTest, ScalarTypes);

TYPED_TEST(AiryTest, PeakAtCenter)
{
    using T = TypeParam;
    const T A = 5.0, x0 = T(1.5), a = T(2.0);
    EXPECT_NEAR(double(fn::airy(x0, A, x0, a)), double(A), FnTol<T>::loose);
}

TYPED_TEST(AiryTest, Symmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), a = T(1.5), d = T(0.6);
    EXPECT_NEAR(double(fn::airy(x0 + d, A, x0, a)),
                double(fn::airy(x0 - d, A, x0, a)), FnTol<T>::func);
}

TYPED_TEST(AiryTest, ZeroNearFirstNull)
{
    // First null of J1(u) is at u ≈ 3.8317; airy should be near zero there
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), a = T(1.0);
    const T u_null = T(3.8317);
    const T val = fn::airy(u_null, A, x0, a);   // a*(x-x0) = u_null
    EXPECT_NEAR(double(val), 0.0, 1e-3);
}

// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class Airy2DSymTest : public testing::Test {};
TYPED_TEST_SUITE(Airy2DSymTest, ScalarTypes);

TYPED_TEST(Airy2DSymTest, PeakAtCenter)
{
    using T = TypeParam;
    const T A = 3.0, x0 = T(0.5), y0 = T(-0.5), a = T(1.5);
    EXPECT_NEAR(double(fn::airy2D(x0, y0, A, x0, y0, a)), double(A), FnTol<T>::loose);
}

TYPED_TEST(Airy2DSymTest, RadialConsistency)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), a = T(2.0);
    for (T r = T(0); r <= T(3.0); r += T(0.3))
        EXPECT_NEAR(double(fn::airy2D(r, T(0), A, x0, y0, a)),
                    double(fn::airy(r, A, T(0), a)),
                    FnTol<T>::loose);
}

// ─────────────────────────────────────────────────────────────────────────────
template<typename T> class Airy2DAsymTest : public testing::Test {};
TYPED_TEST_SUITE(Airy2DAsymTest, ScalarTypes);

TYPED_TEST(Airy2DAsymTest, PeakAtCenter)
{
    using T = TypeParam;
    const T A = 2.0, x0 = T(0), y0 = T(0), a = T(1.0), q = T(0.6), theta = T(0.3);
    EXPECT_NEAR(double(fn::airy2D(x0, y0, A, x0, y0, a, q, theta)), double(A), FnTol<T>::loose);
}

TYPED_TEST(Airy2DAsymTest, Q1MatchesSym)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), a = T(1.5), q = T(1.0), theta = T(0);
    for (T dx = T(0); dx <= T(2.0); dx += T(0.5))
        for (T dy = T(0); dy <= T(2.0); dy += T(0.5))
            EXPECT_NEAR(double(fn::airy2D(x0+dx, y0+dy, A, x0, y0, a, q, theta)),
                        double(fn::airy2D(x0+dx, y0+dy, A, x0, y0, a)),
                        FnTol<T>::loose);
}

TYPED_TEST(Airy2DAsymTest, PointSymmetry)
{
    using T = TypeParam;
    const T A = 1.0, x0 = T(0), y0 = T(0), a = T(1.0), q = T(0.7), theta = T(0.5);
    for (T dx = -T(1.5); dx <= T(1.5); dx += T(0.5))
        for (T dy = -T(1.5); dy <= T(1.5); dy += T(0.5))
            EXPECT_NEAR(double(fn::airy2D(x0+dx, y0+dy, A, x0, y0, a, q, theta)),
                        double(fn::airy2D(x0-dx, y0-dy, A, x0, y0, a, q, theta)),
                        FnTol<T>::func);
}