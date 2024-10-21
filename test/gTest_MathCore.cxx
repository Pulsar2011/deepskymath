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
    EXPECT_DOUBLE_EQ(kBig,4.503599627370496e15);
    EXPECT_DOUBLE_EQ(kBiginv,2.22044604925031308085e-16);
    EXPECT_DOUBLE_EQ(LS2PI,0.91893853320467274178);
    
    EXPECT_DOUBLE_EQ(A[0],8.11614167470508450300e-4);
    EXPECT_DOUBLE_EQ(A[1],-5.95061904284301438324E-4);
    EXPECT_DOUBLE_EQ(A[2],7.93650340457716943945E-4);
    EXPECT_DOUBLE_EQ(A[3],-2.77777777730099687205E-3);
    EXPECT_DOUBLE_EQ(A[4],8.33333333333331927722E-2);
    EXPECT_DOUBLE_EQ(B[0],-1.37825152569120859100E3);
    EXPECT_DOUBLE_EQ(B[1],-3.88016315134637840924E4);
    EXPECT_DOUBLE_EQ(B[2],-3.31612992738871184744E5);
    EXPECT_DOUBLE_EQ(B[3],-1.16237097492762307383E6);
    EXPECT_DOUBLE_EQ(B[4],-1.72173700820839662146E6);
    EXPECT_DOUBLE_EQ(B[5],-8.53555664245765465627E5);
    EXPECT_DOUBLE_EQ(C[0],-3.51815701436523470549E2);
    EXPECT_DOUBLE_EQ(C[1],-1.70642106651881159223E4);
    EXPECT_DOUBLE_EQ(C[2],-2.20528590553854454839E5);
    EXPECT_DOUBLE_EQ(C[3],-1.13933444367982507207E6);
    EXPECT_DOUBLE_EQ(C[4],-2.53252307177582951285E6);
    EXPECT_DOUBLE_EQ(C[5],-2.01889141433532773231E6);

    EXPECT_DOUBLE_EQ(P[0],1.60119522476751861407E-4);
    EXPECT_DOUBLE_EQ(P[1],1.19135147006586384913E-3);
    EXPECT_DOUBLE_EQ(P[2],1.04213797561761569935E-2);
    EXPECT_DOUBLE_EQ(P[3],4.76367800457137231464E-2);
    EXPECT_DOUBLE_EQ(P[4],2.07448227648435975150E-1);
    EXPECT_DOUBLE_EQ(P[5],4.94214826801497100753E-1);
    EXPECT_DOUBLE_EQ(P[6],9.99999999999999996796E-1);

    EXPECT_DOUBLE_EQ(Q[0],-2.31581873324120129819E-5);
    EXPECT_DOUBLE_EQ(Q[1],5.39605580493303397842E-4 );
    EXPECT_DOUBLE_EQ(Q[2],-4.45641913851797240494E-3);
    EXPECT_DOUBLE_EQ(Q[3],1.18139785222060435552E-2 );
    EXPECT_DOUBLE_EQ(Q[4],3.58236398605498653373E-2 );
    EXPECT_DOUBLE_EQ(Q[5],-2.34591795718243348568E-1);
    EXPECT_DOUBLE_EQ(Q[6],7.14304917030273074085E-2 );
    EXPECT_DOUBLE_EQ(Q[7],1.00000000000000000320E0);

    EXPECT_DOUBLE_EQ(STIR[0],7.87311395793093628397E-4);
    EXPECT_DOUBLE_EQ(STIR[1],-2.29549961613378126380E-4);
    EXPECT_DOUBLE_EQ(STIR[2],-2.68132617805781232825E-3);
    EXPECT_DOUBLE_EQ(STIR[3],3.47222221605458667310E-3);
    EXPECT_DOUBLE_EQ(STIR[4],8.33333333333482257126E-2);

    EXPECT_DOUBLE_EQ(igam(1,1),0.63212055882855767840);
    EXPECT_DOUBLE_EQ(igam(1,2),0.86466471676338730811);
    EXPECT_DOUBLE_EQ(igam(1,3),0.95021293163213605700);
    EXPECT_DOUBLE_EQ(igam(1,5),0.99326205300091446621);
    EXPECT_DOUBLE_EQ(igam(2,1),0.26424111765711535688);
    EXPECT_DOUBLE_EQ(igam(2,2),0.59399415029016132000);
    
    EXPECT_DOUBLE_EQ(lgam(1),0.0);
    EXPECT_DOUBLE_EQ(lgam(2),0.0);
    EXPECT_DOUBLE_EQ(lgam(3),0.69314718055994530942);
    EXPECT_DOUBLE_EQ(lgam(4),1.79175946922805500081);
    EXPECT_DOUBLE_EQ(lgam(5),3.17805383034794561965);

    EXPECT_DOUBLE_EQ(igamc(1,1),0.36787944117144232160);
    EXPECT_DOUBLE_EQ(igamc(1,2),0.13533528323661269189);
    EXPECT_DOUBLE_EQ(igamc(1,3),0.04978706836786394300);
    EXPECT_DOUBLE_EQ(igamc(1,5),0.00673794699908546700);
    EXPECT_DOUBLE_EQ(igamc(2,1),0.73575888234288464312);

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

            EXPECT_DOUBLE_EQ(Polynomialeval(x,a,N),pom)<<"P("<<x<<",a,"<<N<<")="<<Polynomialeval(x,a,N)<<" <-> "<<pom<<" ["<<__LINE__<<"]"<<std::endl;
        }

        b[0]=unif(re);
        EXPECT_DOUBLE_EQ(Polynomialeval(x,b, 0),b[0]);
    }

    for(int i = 0; i < 1000; i++)
    {
        double x=unif(re);
        for(int N=0; N<11; N++)
        {
            if(N<1)
            {
                EXPECT_DOUBLE_EQ(Polynomial1eval(x,a,N),a[0])<<"P1("<<x<<",a,"<<N<<")="<<Polynomial1eval(x,a,N)<<" <-> "<<a[0]<<" ["<<__LINE__<<"]"<<std::endl;;
                continue;
            }

            double pom = std::pow(x,N);
            for(int j = 0; j < N; j++)
            {
                pom += a[j]*std::pow(x,N-j-1);
            }

            EXPECT_FLOAT_EQ(Polynomial1eval(x,a,N),pom)<<"P1("<<x<<",a,"<<N<<")="<<Polynomial1eval(x,a,N)<<" <-> "<<pom<<" ["<<__LINE__<<"]"<<std::endl;
        }
    }

    for(int i = 0; i < 1000; i++)
    {
        double z = unif(re);
        double x = unif(re);
        double pom = z/((1.0 + 0.5772156649015329 * x) * x);
        EXPECT_DOUBLE_EQ(small(z,x),pom);
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
    std::default_random_engine re;
    std::vector<double> coef;

    for(unsigned int k =0; k < 5; k++)
        coef.push_back(unif(re));

    for(unsigned int i =0; i < 11; i++)
    {
        if(i==0 || i==10) continue;

        double xi = 2.*static_cast<double>(i)/10.-1.;
        double xo = xi+unif(re)/1e3;

        double y = polynom::chebev(xi, coef, -1., 1.);
        polynom::chebinv(y, xo, coef, -1., 1.);

        EXPECT_NEAR(xi,xo,std::numeric_limits<double>::epsilon()*1e6);
    }
}