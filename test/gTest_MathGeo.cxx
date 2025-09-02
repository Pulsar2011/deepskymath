#include <gtest/gtest.h>
#include <DSTmath/DSTgeo.h>
#include <random>
#include <functional>
#include <iostream>

using namespace DST::Math;

#define NTEST 1000

#pragma region - Test point class
TEST(points, ctor)
{
    point p0 = point();
    EXPECT_EQ(p0.size(), 0);

    std::uniform_real_distribution<double> unid(-100.,100.);
    std::uniform_int_distribution<int32_t> unif(1,100);
    std::default_random_engine re;

    for(size_t k=0; k<NTEST; k++)
    {
        size_t n = unif(re);
        point p1 = point(n);
        EXPECT_EQ(p1.size(), n);

        for(size_t i=0; i<n; i++)
        {
            EXPECT_NEAR(p1[i], 0., point::precision);
        }
    }

    double x = unid(re);
    double y = unid(re);
    double z = unid(re);
    double t = unid(re);
    double u = unid(re);
    double v = unid(re);
    
    point p3 = point(1,x);
    EXPECT_EQ(p3.size(), 1);
    EXPECT_NEAR(p3[0], x, point::precision);

    point cp3 = point(p3);
    EXPECT_EQ(cp3.size(), 1);
    EXPECT_NEAR(cp3[0], x, point::precision);

    point p4 = point(2,x,y);
    EXPECT_EQ(p4.size(), 2);
    EXPECT_NEAR(p4[0], x, point::precision);
    EXPECT_NEAR(p4[1], y, point::precision);

    point cp4 = point(p4);
    EXPECT_EQ(cp4.size(), 2);
    EXPECT_NEAR(cp4[0], x, point::precision);
    EXPECT_NEAR(cp4[1], y, point::precision);

    point p5 = point(3,x,y,z);
    EXPECT_EQ(p5.size(), 3);
    EXPECT_NEAR(p5[0], x, point::precision);
    EXPECT_NEAR(p5[1], y, point::precision);
    EXPECT_NEAR(p5[2], z, point::precision);

    point cp5 = point(p5);
    EXPECT_EQ(cp5.size(), 3);
    EXPECT_NEAR(cp5[0], x, point::precision);
    EXPECT_NEAR(cp5[1], y, point::precision);
    EXPECT_NEAR(cp5[2], z, point::precision);

    point p6 = point(4,x,y,z,t);
    EXPECT_EQ(p6.size(), 4);
    EXPECT_NEAR(p6[0], x, point::precision);
    EXPECT_NEAR(p6[1], y, point::precision);
    EXPECT_NEAR(p6[2], z, point::precision);
    EXPECT_NEAR(p6[3], t, point::precision);

    point cp6 = point(p6);
    EXPECT_EQ(cp6.size(), 4);
    EXPECT_NEAR(cp6[0], x, point::precision);
    EXPECT_NEAR(cp6[1], y, point::precision);
    EXPECT_NEAR(cp6[2], z, point::precision);
    EXPECT_NEAR(cp6[3], t, point::precision);

    point p7 = point(5,x,y,z,t,u);
    EXPECT_EQ(p7.size(), 5);
    EXPECT_NEAR(p7[0], x, point::precision);
    EXPECT_NEAR(p7[1], y, point::precision);
    EXPECT_NEAR(p7[2], z, point::precision);
    EXPECT_NEAR(p7[3], t, point::precision);
    EXPECT_NEAR(p7[4], u, point::precision);

    point cp7 = point(p7);
    EXPECT_EQ(cp7.size(), 5);
    EXPECT_NEAR(cp7[0], x, point::precision);
    EXPECT_NEAR(cp7[1], y, point::precision);
    EXPECT_NEAR(cp7[2], z, point::precision);
    EXPECT_NEAR(cp7[3], t, point::precision);
    EXPECT_NEAR(cp7[4], u, point::precision);

    point p8 = point(6,x,y,z,t,u,v);
    EXPECT_EQ(p8.size(), 6);
    EXPECT_NEAR(p8[0], x, point::precision);
    EXPECT_NEAR(p8[1], y, point::precision);
    EXPECT_NEAR(p8[2], z, point::precision);
    EXPECT_NEAR(p8[3], t, point::precision);
    EXPECT_NEAR(p8[4], u, point::precision);
    EXPECT_NEAR(p8[5], v, point::precision);

    point cp8 = point(p8);
    EXPECT_EQ(cp8.size(), 6);
    EXPECT_NEAR(cp8[0], x, point::precision);
    EXPECT_NEAR(cp8[1], y, point::precision);
    EXPECT_NEAR(cp8[2], z, point::precision);
    EXPECT_NEAR(cp8[3], t, point::precision);
    EXPECT_NEAR(cp8[4], u, point::precision);
    EXPECT_NEAR(cp8[5], v, point::precision);

    point p10 = point(3,0.,0.,0.);
    EXPECT_EQ(p10.size(), 3);
    EXPECT_NEAR(p10[0], 0., point::precision);
    EXPECT_NEAR(p10[1], 0., point::precision);
    EXPECT_NEAR(p10[2], 0., point::precision);

    p10 = p3;
    EXPECT_EQ(p10.size(), p3.size());
    EXPECT_NEAR(p10[0], x, point::precision);

    p10=p8;
    EXPECT_EQ(p10.size(), p8.size());
    EXPECT_NEAR(p10[0], x, point::precision);
    EXPECT_NEAR(p10[1], y, point::precision);
    EXPECT_NEAR(p10[2], z, point::precision);
    EXPECT_NEAR(p10[3], t, point::precision);
    EXPECT_NEAR(p10[4], u, point::precision);
    EXPECT_NEAR(p10[5], v, point::precision);
}

TEST(points, setpoints)
{
    point p  = point(3);
    EXPECT_EQ(p.size(), 3);

    point p0 = point(3);
    EXPECT_EQ(p0.size(), 3);

    point p1 = point(4);
    EXPECT_EQ(p1.size(), 4);

    point p2 = point(6);
    EXPECT_EQ(p2.size(), 6);

    std::uniform_real_distribution<double> unid(-100.,100.);
    std::default_random_engine re;

    for(size_t k=0; k<NTEST; k++)
    {
        double x = unid(re);
        double y = unid(re);
        double z = unid(re);
        double t = unid(re);
        
        p.SetPoints({z,y,x});
        EXPECT_NEAR(p[0], z, point::precision);
        EXPECT_NEAR(p[1], y, point::precision);
        EXPECT_NEAR(p[2], x, point::precision);

        p0.SetPoints({x,y,z});
        EXPECT_NEAR(p0[0], x, point::precision);
        EXPECT_NEAR(p0[1], y, point::precision);
        EXPECT_NEAR(p0[2], z, point::precision);

        p0.SetPoints({y,z,t,x});
        EXPECT_NEAR(p0[0], y, point::precision);
        EXPECT_NEAR(p0[1], z, point::precision);
        EXPECT_NEAR(p0[2], t, point::precision);

        p1.SetPoints({x,y,z,t});
        EXPECT_NEAR(p1[0], x, point::precision);
        EXPECT_NEAR(p1[1], y, point::precision);
        EXPECT_NEAR(p1[2], z, point::precision);
        EXPECT_NEAR(p1[3], t, point::precision);

        EXPECT_ANY_THROW(p2.SetPoints({x,y,z}));

        p2.SetPoints(std::vector<double>({x,y,z,t,x,y}));
        EXPECT_NEAR(p2[0], x, point::precision);
        EXPECT_NEAR(p2[1], y, point::precision);
        EXPECT_NEAR(p2[2], z, point::precision);
        EXPECT_NEAR(p2[3], t, point::precision);
        EXPECT_NEAR(p2[4], x, point::precision);
        EXPECT_NEAR(p2[5], y, point::precision);

        p2.SetPoints(std::vector<double>({x,x,x}));
        EXPECT_NEAR(p2[0], x, point::precision);
        EXPECT_NEAR(p2[1], x, point::precision);
        EXPECT_NEAR(p2[2], x, point::precision);
        EXPECT_NEAR(p2[3], t, point::precision);
        EXPECT_NEAR(p2[4], x, point::precision);
        EXPECT_NEAR(p2[5], y, point::precision);

        p2.SetPoints(std::vector<double>({y,y,y,y,y,y,y,y,y}));
        EXPECT_EQ(p2.size(), 6);
        EXPECT_NEAR(p2[0], y, point::precision);
        EXPECT_NEAR(p2[1], y, point::precision);
        EXPECT_NEAR(p2[2], y, point::precision);
        EXPECT_NEAR(p2[3], y, point::precision);
        EXPECT_NEAR(p2[4], y, point::precision);
        EXPECT_NEAR(p2[5], y, point::precision);

        EXPECT_ANY_THROW(p2.SetPoints(p0));
        
        p0.SetPoints(p);
        EXPECT_EQ(p0.size(), 3);
        EXPECT_NEAR(p0[0], p[0], point::precision);
        EXPECT_NEAR(p0[1], p[1], point::precision);
        EXPECT_NEAR(p0[2], p[2], point::precision);

        p1.SetPoint(0,0.);
        p1.SetPoint(1,0.);
        p1.SetPoint(2,0.);
        p1.SetPoint(3,0.);
        EXPECT_NEAR(p1[0], 0.0, point::precision);
        EXPECT_NEAR(p1[1], 0.0, point::precision);
        EXPECT_NEAR(p1[2], 0.0, point::precision);
        EXPECT_NEAR(p1[2], 0.0, point::precision);

        EXPECT_ANY_THROW(p1.SetPoint(4,0.));

        p1.SetX(x);
        p1.SetY(y);
        p1.SetZ(z);

        EXPECT_NEAR(p1[0], x, point::precision);
        EXPECT_NEAR(p1[1], y, point::precision);
        EXPECT_NEAR(p1[2], z, point::precision);
        EXPECT_NEAR(p1[3], 0.0, point::precision);

        p0.SetPoints(std::vector<double>({x,y,z}));
        EXPECT_NEAR(p0[0], x, point::precision);
        EXPECT_NEAR(p0[1], y, point::precision);
        EXPECT_NEAR(p0[2], z, point::precision);

        p0.SetPoints(std::vector<float>({static_cast<float>(x),static_cast<float>(y),static_cast<float>(z)}));
        EXPECT_NEAR(p0[0], x,1e-4);
        EXPECT_NEAR(p0[1], y,1e-4);
        EXPECT_NEAR(p0[2], z,1e-4);

        p0.SetPoints({static_cast<float>(z),static_cast<float>(x),static_cast<float>(t)});
        EXPECT_NEAR(p0[0], z,1e-4);
        EXPECT_NEAR(p0[1], x,1e-4);
        EXPECT_NEAR(p0[2], t,1e-4);
    }

}

TEST(point, accessor)
{
    point p = point(3);
    EXPECT_EQ(p.size(), 3);

    std::uniform_real_distribution<double> unid(-100.,100.);
    std::default_random_engine re;

    for(size_t k=0; k<NTEST; k++)
    {
        double x = unid(re);
        double y = unid(re);
        double z = unid(re);
        
        p.SetPoints({x,y,z});
        EXPECT_NEAR(p.X(), x, point::precision);
        EXPECT_NEAR(p.Y(), y, point::precision);
        EXPECT_NEAR(p.Z(), z, point::precision);
        EXPECT_ANY_THROW(p[3]);

        EXPECT_NEAR(p.R(), sqrt(x*x + y*y + z*z), point::precision);
        EXPECT_NEAR(p.Phi(),((atan2(y,x)>0)?atan2(y,x):2.*acos(-1) - atan2(y,x)), point::precision);
        EXPECT_NEAR(p.Theta(), acos(z/p.R()), point::precision);
    }
}

TEST(point, operator)
{
    point p1 = point(3);
    point p2 = point(3);
    point p4 = point(4);

    std::uniform_real_distribution<double> unid(-100.,100.);
    std::default_random_engine re;

    for(size_t k=0; k<NTEST; k++)
    {
        double x1 = unid(re);
        double y1 = unid(re);
        double z1 = unid(re);

        double x2 = unid(re);
        double y2 = unid(re);
        double z2 = unid(re);

        double t = unid(re);
        
        p1.SetPoints({x1,y1,z1});
        p2.SetPoints({x2,y2,z2});
        p4.SetPoints({x2,y2,z2,t});

        EXPECT_FALSE(p1 == p4);
        EXPECT_TRUE (p1 != p4);

        point cp1 = point(p1);
        point cp2 = point(p2);

        if( std::abs(x1 - x2) <= point::precision && std::abs(y1 - y2) <= point::precision && std::abs(z1 - z2) <= point::precision )
        {
            EXPECT_TRUE (p1 == p2) << " ("
                <<p1[0]<<","<<p1[1]<<","<<p1[2]<<") != ("
                <<p2[0]<<","<<p2[1]<<","<<p2[2]<<")";
            EXPECT_FALSE(p1 != p2) << " ("
                <<p1[0]<<","<<p1[1]<<","<<p1[2]<<") == ("
                <<p2[0]<<","<<p2[1]<<","<<p2[2]<<")";
        }
        else
        {
            EXPECT_FALSE(p1 == p2) << " ("
                <<p1[0]<<","<<p1[1]<<","<<p1[2]<<") != ("
                <<p2[0]<<","<<p2[1]<<","<<p2[2]<<")";
            EXPECT_TRUE(p1 != p2) << " ("
                <<p1[0]<<","<<p1[1]<<","<<p1[2]<<") == ("
                <<p2[0]<<","<<p2[1]<<","<<p2[2]<<")";
        }

        if(p1.R() > p2.R())
        {
            EXPECT_TRUE (p1 >  p2)<<p1.Dump()<<" >  "<<p2.Dump();
            EXPECT_FALSE(p1 <  p2)<<p1.Dump()<<" <  "<<p2.Dump();
            EXPECT_TRUE (p1 >= p2)<<p1.Dump()<<" >= "<<p2.Dump();
            EXPECT_FALSE(p1 <= p2)<<p1.Dump()<<" <= "<<p2.Dump();
        }
        else if (p1.Phi() > p2.Phi())
        {
            EXPECT_TRUE (p1 >  p2)<<p1.Dump()<<" >  "<<p2.Dump();
            EXPECT_FALSE(p1 <  p2)<<p1.Dump()<<" <  "<<p2.Dump();
            EXPECT_TRUE (p1 >= p2)<<p1.Dump()<<" >= "<<p2.Dump();
            EXPECT_FALSE(p1 <= p2)<<p1.Dump()<<" <= "<<p2.Dump();
        }
        else if (p1.Theta() > p2.Theta())
        {
            EXPECT_TRUE (p1 >  p2)<<p1.Dump()<<" >  "<<p2.Dump();
            EXPECT_FALSE(p1 <  p2)<<p1.Dump()<<" <  "<<p2.Dump();
            EXPECT_TRUE (p1 >= p2)<<p1.Dump()<<" >= "<<p2.Dump();
            EXPECT_FALSE(p1 <= p2)<<p1.Dump()<<" <= "<<p2.Dump();
        }
        else 
        {
            EXPECT_FALSE(p1 >  p2)<<p1.Dump()<<" >  "<<p2.Dump();
            EXPECT_TRUE (p1 <  p2)<<p1.Dump()<<" <  "<<p2.Dump();
            EXPECT_FALSE(p1 >= p2)<<p1.Dump()<<" >= "<<p2.Dump();
            EXPECT_TRUE (p1 <= p2)<<p1.Dump()<<" <= "<<p2.Dump();
        }



        if(p1.R() > p4.R())
        {
            EXPECT_TRUE (p1 >  p4)<<p1.Dump()<<" >  "<<p4.Dump();
            EXPECT_FALSE(p1 <  p4)<<p1.Dump()<<" <  "<<p4.Dump();
            EXPECT_TRUE (p1 >= p4)<<p1.Dump()<<" >= "<<p4.Dump();
            EXPECT_FALSE(p1 <= p4)<<p1.Dump()<<" <= "<<p4.Dump();
        }
        else if (p1.Phi() > p4.Phi())
        {
            EXPECT_TRUE (p1 >  p4)<<p1.Dump()<<" >  "<<p4.Dump();
            EXPECT_FALSE(p1 <  p4)<<p1.Dump()<<" <  "<<p4.Dump();
            EXPECT_TRUE (p1 >= p4)<<p1.Dump()<<" >= "<<p4.Dump();
            EXPECT_FALSE(p1 <= p4)<<p1.Dump()<<" <= "<<p4.Dump();
        }
        else if (p1.Theta() > p4.Theta())
        {
            EXPECT_TRUE (p1 >  p4)<<p1.Dump()<<" >  "<<p4.Dump();
            EXPECT_FALSE(p1 <  p4)<<p1.Dump()<<" <  "<<p4.Dump();
            EXPECT_TRUE (p1 >= p4)<<p1.Dump()<<" >= "<<p4.Dump();
            EXPECT_FALSE(p1 <= p4)<<p1.Dump()<<" <= "<<p4.Dump();
        }
        else if (p4[3] < point::precision)
        {
            EXPECT_TRUE (p1 >  p4)<<p1.Dump()<<" >  "<<p4.Dump();
            EXPECT_FALSE(p1 <  p4)<<p1.Dump()<<" <  "<<p4.Dump();
            EXPECT_TRUE (p1 >= p4)<<p1.Dump()<<" >= "<<p4.Dump();
            EXPECT_FALSE(p1 <= p4)<<p1.Dump()<<" <= "<<p4.Dump();
        }
        else 
        {
            EXPECT_FALSE(p1 >  p4)<<p1.Dump()<<" >  "<<p4.Dump();
            EXPECT_TRUE (p1 <  p4)<<p1.Dump()<<" <  "<<p4.Dump();
            EXPECT_FALSE(p1 >= p4)<<p1.Dump()<<" >= "<<p4.Dump();
            EXPECT_TRUE (p1 <= p4)<<p1.Dump()<<" <= "<<p4.Dump();
        }

        if(p2.R() > p4.R())
        {
            EXPECT_TRUE (p2 >  p4)<<p2.Dump()<<" >  "<<p4.Dump();
            EXPECT_FALSE(p2 <  p4)<<p2.Dump()<<" <  "<<p4.Dump();
            EXPECT_TRUE (p2 >= p4)<<p2.Dump()<<" >= "<<p4.Dump();
            EXPECT_FALSE(p2 <= p4)<<p2.Dump()<<" <= "<<p4.Dump();
        }
        else if (p2.Phi() > p4.Phi())
        {
            EXPECT_TRUE (p2 >  p4)<<p2.Dump()<<" >  "<<p4.Dump();
            EXPECT_FALSE(p2 <  p4)<<p2.Dump()<<" <  "<<p4.Dump();
            EXPECT_TRUE (p2 >= p4)<<p2.Dump()<<" >= "<<p4.Dump();
            EXPECT_FALSE(p2 <= p4)<<p2.Dump()<<" <= "<<p4.Dump();
        }
        else if (p2.Theta() > p4.Theta())
        {
            EXPECT_TRUE (p2 >  p4)<<p2.Dump()<<" >  "<<p4.Dump();
            EXPECT_FALSE(p2 <  p4)<<p2.Dump()<<" <  "<<p4.Dump();
            EXPECT_TRUE (p2 >= p4)<<p2.Dump()<<" >= "<<p4.Dump();
            EXPECT_FALSE(p2 <= p4)<<p2.Dump()<<" <= "<<p4.Dump();
        }
        else if (p4[3] < point::precision)
        {
            EXPECT_TRUE (p2 >  p4)<<p2.Dump()<<" >  "<<p4.Dump();
            EXPECT_FALSE(p2 <  p4)<<p2.Dump()<<" <  "<<p4.Dump();
            EXPECT_TRUE (p2 >= p4)<<p2.Dump()<<" >= "<<p4.Dump();
            EXPECT_FALSE(p2 <= p4)<<p2.Dump()<<" <= "<<p4.Dump();
        }
        else 
        {
            EXPECT_FALSE(p2 >  p4)<<p2.Dump()<<" >  "<<p4.Dump();
            EXPECT_TRUE (p2 <  p4)<<p2.Dump()<<" <  "<<p4.Dump();
            EXPECT_FALSE(p2 >= p4)<<p2.Dump()<<" >= "<<p4.Dump();
            EXPECT_TRUE (p2 <= p4)<<p2.Dump()<<" <= "<<p4.Dump();
        }

        EXPECT_TRUE (p1 == p1);
        EXPECT_FALSE(p1 != p1);
        EXPECT_TRUE (p1 >= p1);
        EXPECT_TRUE (p1 <= p1);
        EXPECT_FALSE(p1  < p1);
        EXPECT_FALSE(p1  > p1);

        EXPECT_TRUE (cp1 == p1);
        EXPECT_FALSE(cp1 != p1);
        EXPECT_TRUE (cp1 >= p1);
        EXPECT_TRUE (cp1 <= p1);
        EXPECT_FALSE(cp1  < p1);
        EXPECT_FALSE(cp1  > p1);
    }
}

TEST(point, arithmetic)
{

    point::precision = 1e-9;
    point p0 = point(3);
    point p1 = point(3);
    point p2 = point(3);
    point p3 = point(4);

    std::uniform_real_distribution<double> unid(-100.,100.);
    std::default_random_engine re;

    for(size_t k=0; k<NTEST; k++)
    {
        double x1 = unid(re);
        double y1 = unid(re);
        double z1 = unid(re);

        double x2 = unid(re);
        double y2 = unid(re);
        double z2 = unid(re);

        double x3 = unid(re);
        double y3 = unid(re);
        double z3 = unid(re);
        double t3 = unid(re);

        p0.SetPoints({x1,y2,z3});
        p1.SetPoints({x1,y1,z1});
        p2.SetPoints({x2,y2,z2});
        p3.SetPoints({x3,y3,z3,t3});

        point cp0= point(p0);
        point cp1= point(p1);
        point cp2= point(p2);
        point cp3= point(p3);

        p1 += p2;
        EXPECT_NEAR(p1[0], cp1[0]+cp2[0], point::precision);
        EXPECT_NEAR(p1[1], cp1[1]+cp2[1], point::precision);
        EXPECT_NEAR(p1[2], cp1[2]+cp2[2], point::precision);

        p1 = cp1;
        p1 -= p2;
        EXPECT_NEAR(p1[0], cp1[0]-cp2[0], point::precision);
        EXPECT_NEAR(p1[1], cp1[1]-cp2[1], point::precision);
        EXPECT_NEAR(p1[2], cp1[2]-cp2[2], point::precision);

        p1 = cp1;
        p1 *= p2;
        EXPECT_NEAR(p1[0], cp1[0]*cp2[0], point::precision);
        EXPECT_NEAR(p1[1], cp1[1]*cp2[1], point::precision);
        EXPECT_NEAR(p1[2], cp1[2]*cp2[2], point::precision);

        p1 = cp1;
        p1 /= p2;
        EXPECT_NEAR(p1[0], cp1[0]/cp2[0], point::precision);
        EXPECT_NEAR(p1[1], cp1[1]/cp2[1], point::precision);
        EXPECT_NEAR(p1[2], cp1[2]/cp2[2], point::precision);

        p1=cp1;
        p3+=p1;
        EXPECT_NEAR(p3[0], cp1[0]+cp3[0], point::precision);
        EXPECT_NEAR(p3[1], cp1[1]+cp3[1], point::precision);
        EXPECT_NEAR(p3[2], cp1[2]+cp3[2], point::precision);
        EXPECT_NEAR(p3[3], cp3[3], point::precision);

        p3=cp3;
        p3-=p1;
        EXPECT_NEAR(p3[0], cp3[0]-cp1[0], point::precision);
        EXPECT_NEAR(p3[1], cp3[1]-cp1[1], point::precision);
        EXPECT_NEAR(p3[2], cp3[2]-cp1[2], point::precision);
        EXPECT_NEAR(p3[3], cp3[3], point::precision);

        p3=cp3;
        p3*=p1;
        EXPECT_NEAR(p3[0], cp3[0]*cp1[0], point::precision);
        EXPECT_NEAR(p3[1], cp3[1]*cp1[1], point::precision);
        EXPECT_NEAR(p3[2], cp3[2]*cp1[2], point::precision);
        EXPECT_NEAR(p3[3], cp3[3], point::precision);

        p3=cp3;
        p3/=p1;
        EXPECT_NEAR(p3[0], cp3[0]/cp1[0], point::precision);
        EXPECT_NEAR(p3[1], cp3[1]/cp1[1], point::precision);
        EXPECT_NEAR(p3[2], cp3[2]/cp1[2], point::precision);
        EXPECT_NEAR(p3[3], cp3[3], point::precision);

        p1=cp1;
        p1 += t3;
        EXPECT_NEAR(p1[0], cp1[0]+t3, point::precision);
        EXPECT_NEAR(p1[1], cp1[1]+t3, point::precision);
        EXPECT_NEAR(p1[2], cp1[2]+t3, point::precision);

        p1=cp1;
        p1 -= t3;
        EXPECT_NEAR(p1[0], cp1[0]-t3, point::precision);
        EXPECT_NEAR(p1[1], cp1[1]-t3, point::precision);
        EXPECT_NEAR(p1[2], cp1[2]-t3, point::precision);

        p1=cp1;
        p1 *= t3;
        EXPECT_NEAR(p1[0], cp1[0]*t3, point::precision);
        EXPECT_NEAR(p1[1], cp1[1]*t3, point::precision);
        EXPECT_NEAR(p1[2], cp1[2]*t3, point::precision);

        p1=cp1;
        p1 /= t3;
        EXPECT_NEAR(p1[0], cp1[0]/t3, point::precision);
        EXPECT_NEAR(p1[1], cp1[1]/t3, point::precision);
        EXPECT_NEAR(p1[2], cp1[2]/t3, point::precision);

        point p4 = cp1 + cp2;
        EXPECT_NEAR(p4[0], cp1[0]+cp2[0], point::precision);
        EXPECT_NEAR(p4[1], cp1[1]+cp2[1], point::precision);
        EXPECT_NEAR(p4[2], cp1[2]+cp2[2], point::precision);

        p4 = cp1 - cp2;
        EXPECT_NEAR(p4[0], cp1[0]-cp2[0], point::precision);
        EXPECT_NEAR(p4[1], cp1[1]-cp2[1], point::precision);
        EXPECT_NEAR(p4[2], cp1[2]-cp2[2], point::precision);

        p4 = cp1 * cp2;
        EXPECT_NEAR(p4[0], cp1[0]*cp2[0], point::precision);
        EXPECT_NEAR(p4[1], cp1[1]*cp2[1], point::precision);
        EXPECT_NEAR(p4[2], cp1[2]*cp2[2], point::precision);

        p4 = cp1 / cp2;
        EXPECT_NEAR(p4[0], cp1[0]/cp2[0], point::precision);
        EXPECT_NEAR(p4[1], cp1[1]/cp2[1], point::precision);
        EXPECT_NEAR(p4[2], cp1[2]/cp2[2], point::precision);

        EXPECT_ANY_THROW(cp1 + cp3);
        EXPECT_ANY_THROW(cp1 - cp3);
        EXPECT_ANY_THROW(cp1 * cp3);
        EXPECT_ANY_THROW(cp1 / cp3);

        EXPECT_ANY_THROW(cp3 + cp1);
        EXPECT_ANY_THROW(cp3 - cp1);
        EXPECT_ANY_THROW(cp3 * cp1);
        EXPECT_ANY_THROW(cp3 / cp1);

        p4=cp1;
        p4 += t3;
        EXPECT_NEAR(p4[0], cp1[0]+t3, point::precision);
        EXPECT_NEAR(p4[1], cp1[1]+t3, point::precision);
        EXPECT_NEAR(p4[2], cp1[2]+t3, point::precision);

        p4=cp1;
        p4 -= t3;
        EXPECT_NEAR(p4[0], cp1[0]-t3, point::precision);
        EXPECT_NEAR(p4[1], cp1[1]-t3, point::precision);
        EXPECT_NEAR(p4[2], cp1[2]-t3, point::precision);

        p4=cp1;
        p4 *= t3;
        EXPECT_NEAR(p4[0], cp1[0]*t3, point::precision);
        EXPECT_NEAR(p4[1], cp1[1]*t3, point::precision);
        EXPECT_NEAR(p4[2], cp1[2]*t3, point::precision);

        p4=cp1;
        p4 /= t3;
        EXPECT_NEAR(p4[0], cp1[0]/t3, point::precision);
        EXPECT_NEAR(p4[1], cp1[1]/t3, point::precision);
        EXPECT_NEAR(p4[2], cp1[2]/t3, point::precision);

        p0 += t3;
        EXPECT_NEAR(p0[0], cp0[0]+t3, point::precision);
        EXPECT_NEAR(p0[1], cp0[1]+t3, point::precision);
        EXPECT_NEAR(p0[2], cp0[2]+t3, point::precision);

        p0=cp0;
        p0 -= t3;
        EXPECT_NEAR(p0[0], cp0[0]-t3, point::precision);
        EXPECT_NEAR(p0[1], cp0[1]-t3, point::precision);
        EXPECT_NEAR(p0[2], cp0[2]-t3, point::precision);

        p0=cp0;
        p0 *= t3;
        EXPECT_NEAR(p0[0], cp0[0]*t3, point::precision);
        EXPECT_NEAR(p0[1], cp0[1]*t3, point::precision);
        EXPECT_NEAR(p0[2], cp0[2]*t3, point::precision);

        p0=cp0;
        p0 /= t3;
        EXPECT_NEAR(p0[0], cp0[0]/t3, point::precision);
        EXPECT_NEAR(p0[1], cp0[1]/t3, point::precision);
        EXPECT_NEAR(p0[2], cp0[2]/t3, point::precision);

        p0 = cp0;
        point pp = cp0 + t3;
        EXPECT_NEAR(pp[0], cp0[0]+t3, point::precision);
        EXPECT_NEAR(pp[1], cp0[1]+t3, point::precision);
        EXPECT_NEAR(pp[2], cp0[2]+t3, point::precision);

        pp = cp0 - t3;
        EXPECT_NEAR(pp[0], cp0[0]-t3, point::precision);
        EXPECT_NEAR(pp[1], cp0[1]-t3, point::precision);
        EXPECT_NEAR(pp[2], cp0[2]-t3, point::precision);

        pp = cp0 * t3;
        EXPECT_NEAR(pp[0], cp0[0]*t3, point::precision);
        EXPECT_NEAR(pp[1], cp0[1]*t3, point::precision);
        EXPECT_NEAR(pp[2], cp0[2]*t3, point::precision);

        pp = cp0 / t3;
        EXPECT_NEAR(pp[0], cp0[0]/t3, point::precision);
        EXPECT_NEAR(pp[1], cp0[1]/t3, point::precision);
        EXPECT_NEAR(pp[2], cp0[2]/t3, point::precision);
    }
}

#pragma endregion
#pragma region - Test vector2D
TEST(vector2D, ctor)
{
    vector2D v0 = vector2D();
    EXPECT_NEAR(v0.Phi(), 0.  ,point::precision);
    EXPECT_NEAR(v0.Length(),1.,point::precision);

    std::uniform_real_distribution<double> unid(-100.,100.);
    std::default_random_engine re;

    for(size_t k=0; k<NTEST; k++)
    {
        double l = unid(re);
        double p = unid(re);
        
        vector2D vp = vector2D(p);
        
        double twopi = 2. * acos(-1.);
        double phi = std::fmod(p, twopi);
        if(phi < 0) phi += twopi;

        EXPECT_NEAR(vp.Phi(), phi, point::precision);
        EXPECT_NEAR(vp.Length(), 1., point::precision);
        EXPECT_NEAR(vp.X(), cos(phi), point::precision);
        EXPECT_NEAR(vp.Y(), sin(phi), point::precision);

        vector2D vl = vector2D(l, p);
        EXPECT_NEAR(vl.Phi(), phi, point::precision);
        EXPECT_NEAR(vl.Length(), std::abs(l), point::precision);
        EXPECT_NEAR(vl.X(), std::abs(l)*cos(phi), point::precision);
        EXPECT_NEAR(vl.Y(), std::abs(l)*sin(phi), point::precision);
        EXPECT_NEAR(vl.Z(), 0.                  , point::precision);

        vector2D cvl = vector2D(vl);
        EXPECT_NEAR(cvl.Phi(), vl.Phi(), point::precision);
        EXPECT_NEAR(cvl.Length(), vl.Length(), point::precision);
        EXPECT_NEAR(cvl.X(), vl.X(), point::precision);
        EXPECT_NEAR(cvl.Y(), vl.Y(), point::precision);
        EXPECT_NEAR(cvl.Z(), vl.Z(), point::precision);

        double x1 = unid(re);
        double y1 = unid(re);
        p = (atan2(y1,x1)>=0)?atan2(y1,x1):2.*acos(-1) + atan2(y1,x1);
        l = sqrt(x1*x1 + y1*y1);

        vector2D vxy1 = vector2D(point(2,x1,y1));
        EXPECT_NEAR(vxy1.Phi(),    p, point::precision);
        EXPECT_NEAR(vxy1.Length(), l, point::precision);
        EXPECT_NEAR(vxy1.X(),     x1, point::precision);
        EXPECT_NEAR(vxy1.Y(),     y1, point::precision);   
        EXPECT_NEAR(vxy1.Z(),     0., point::precision);

        double x2 = unid(re);
        double y2 = unid(re);
        p = (atan2(y2-y1,x2-x1)>=0)?atan2(y2-y1,x2-x1):2.*acos(-1) + atan2(y2-y1,x2-x1);
        l = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

        vector2D vxy2 = vector2D(point(2,x1,y1),point(2,x2,y2));
        EXPECT_NEAR(vxy2.Phi(), p, point::precision);
        EXPECT_NEAR(vxy2.Length(), l, point::precision);
        EXPECT_NEAR(vxy2.X(), x2-x1, point::precision);
        EXPECT_NEAR(vxy2.Y(), y2-y1, point::precision);   
        EXPECT_NEAR(vxy2.Z(), 0, point::precision);

        vector2D vxy3 = vector2D(point(2,x2,y2),point(2,x1,y1));
        p = (atan2(y1-y2,x1-x2)>=0)?atan2(y1-y2,x1-x2):2.*acos(-1) + atan2(y1-y2,x1-x2);
        l = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        EXPECT_NEAR(vxy3.Phi(), p, point::precision);
        EXPECT_NEAR(vxy3.Length(), l, point::precision);
        EXPECT_NEAR(vxy3.X(), x1-x2, point::precision);
        EXPECT_NEAR(vxy3.Y(), y1-y2, point::precision);   
        EXPECT_NEAR(vxy3.Z(), 0, point::precision);
    }
}

TEST(vector2D,exception)
{
    std::uniform_real_distribution<double> unid(-100.,100.);
    std::default_random_engine re;

    double x1= unid(re);
    double y1= unid(re);
    double x2= unid(re);
    double y2= unid(re);

    EXPECT_ANY_THROW(vector2D(point(1,x1)));
    EXPECT_ANY_THROW(vector2D(point(1,x1),point(1,x2)));
    EXPECT_ANY_THROW(vector2D(point(2,x1,y1),point(1,x2)));
    EXPECT_ANY_THROW(vector2D(point(1,x1),point(2,x2,y2)));
}

TEST(vector2D,modifier)
{
    double twopi = 2.*acos(-1.);

    std::uniform_real_distribution<double> unid(-100.,100.);
    std::default_random_engine re;

    double l = unid(re);
    double p = unid(re);
    
    double phi = std::fmod(p, twopi);
    if(phi < 0) phi += twopi;

    vector2D v = vector2D(l,p);
    EXPECT_NEAR(v.Length(), std::abs(l), point::precision);
    EXPECT_NEAR(v.Phi(), phi, point::precision);
    EXPECT_NEAR(v.X(), std::abs(l)*cos(phi), point::precision);
    EXPECT_NEAR(v.Y(), std::abs(l)*sin(phi), point::precision);

    for(size_t k=0; k<NTEST; k++)
    {
        double p0 = v.Phi();

        l = unid(re);
        p = unid(re);
        phi = std::fmod(p, twopi);
        if(phi < 0) phi += twopi;

        v.SetLength(l);
        EXPECT_NEAR(v.Length(), std::abs(l), point::precision);    
        EXPECT_NEAR(v.Phi(), p0, point::precision);
        EXPECT_NEAR(v.X(), std::abs(l)*cos(p0), point::precision);
        EXPECT_NEAR(v.Y(), std::abs(l)*sin(p0), point::precision);
        EXPECT_NEAR(v.Z(), 0, point::precision);

        v.SetPhi(p);
        EXPECT_NEAR(v.Length(), std::abs(l), point::precision);    
        EXPECT_NEAR(v.Phi(), phi, point::precision);
        EXPECT_NEAR(v.X(), std::abs(l)*cos(phi), point::precision);
        EXPECT_NEAR(v.Y(), std::abs(l)*sin(phi), point::precision);
        EXPECT_NEAR(v.Z(), 0, point::precision);

        point pt = v.EndPoint();
        EXPECT_EQ(pt.size(), 2);
        EXPECT_NEAR(pt[0], v.X(), point::precision);
        EXPECT_NEAR(pt[1], v.Y(), point::precision);        
    }
}

TEST(vector2D,operator)
{
    std::uniform_real_distribution<double> unidl(0,1000.);
    std::uniform_real_distribution<double> unidp(0,2.*acos(-1));
    std::default_random_engine re;

    for(size_t k=0; k<NTEST; k++)
    {
        double l1 = unidl(re);
        double p1 = unidp(re);

        double l2 = unidl(re);
        double p2 = unidp(re);
        
        vector2D v1 = vector2D(l1,p1);
        vector2D v2 = vector2D(l2,p2);
        vector2D v3 = vector2D(l1,p2);

        EXPECT_TRUE(v1 == v1);
        EXPECT_TRUE(v1 >= v1);
        EXPECT_TRUE(v1 <= v1);
        EXPECT_TRUE(v2 == v2);
        EXPECT_TRUE(v2 >= v2);
        EXPECT_TRUE(v2 <= v2);

        if((l1 != l2) || (p1 != p2))
        {
            EXPECT_FALSE(v1 == v2);
            EXPECT_FALSE(v1 == v2);
        }
        else
        {
            EXPECT_TRUE(v1 == v2);
            EXPECT_TRUE(v2 == v1);
        }

        if(std::abs(l1 - l2) < point::precision)
        {
            if(p1 - p2 > point::precision)
            {
                EXPECT_TRUE (v1 > v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_FALSE(v1 < v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_TRUE (v2 < v1)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_FALSE(v2 > v1)<<v1.Dump()<<" >  "<<v2.Dump();
            }
            else
            {
                EXPECT_TRUE (v1 < v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_FALSE(v1 > v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_TRUE (v2 > v1)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_FALSE(v2 < v1)<<v1.Dump()<<" >  "<<v2.Dump();
            }
        }
        else if(l1 - l2 > point::precision)
        {
            EXPECT_TRUE (v1 > v2)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_FALSE(v1 < v2)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_TRUE (v2 < v1)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_FALSE(v2 > v1)<<v1.Dump()<<" >  "<<v2.Dump();
        }
        else
        {
            EXPECT_TRUE (v1 < v2)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_FALSE(v1 > v2)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_TRUE (v2 > v1)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_FALSE(v2 < v1)<<v1.Dump()<<" >  "<<v2.Dump();
        }

        if(p1 - p2 > point::precision)
        {
            EXPECT_TRUE (v1 > v3)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_FALSE(v1 < v3)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_TRUE (v3 < v1)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_FALSE(v3 > v1)<<v1.Dump()<<" >  "<<v3.Dump();
        }
        else
        {
            EXPECT_TRUE (v1 < v3)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_FALSE(v1 > v3)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_TRUE (v3 > v1)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_FALSE(v3 < v1)<<v1.Dump()<<" >  "<<v3.Dump();
        }
    }
}

TEST(vector2D, arithmetic)
{
    point::precision = 1e-9;

    double twopi = 2.*acos(-1.);
    std::uniform_real_distribution<double> real(-10,10);
    std::uniform_real_distribution<double> unidl(0,1000.);
    std::uniform_real_distribution<double> unidp(0,2.*acos(-1));

    std::default_random_engine re;

    for(size_t k=0; k<NTEST; k++)
    {
        double l1 = unidl(re);
        double p1 = unidp(re);

        double l2 = unidl(re);
        double p2 = unidp(re);

        double s = real(re);
        
        vector2D v1 = vector2D(l1,p1);
        vector2D v2 = vector2D(l2,p2);

        vector2D cv1 = vector2D(v1);
        vector2D cv2 = vector2D(v2);

        v1 += v2;
        EXPECT_NEAR(v1.X(), cv1.X()+cv2.X(), point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Y()+cv2.Y(), point::precision);
        EXPECT_NEAR(v1.Length(), sqrt((cv1.X()+cv2.X())*(cv1.X()+cv2.X()) + (cv1.Y()+cv2.Y())*(cv1.Y()+cv2.Y())), point::precision);
        EXPECT_NEAR(v1.Phi(), (atan2(cv1.Y()+cv2.Y(),cv1.X()+cv2.X())>=0)?atan2(cv1.Y()+cv2.Y(),cv1.X()+cv2.X()):twopi+atan2(cv1.Y()+cv2.Y(),cv1.X()+cv2.X()), point::precision);

        v1 = cv1;
        v2 = cv2;
        v1 -= v2;
        EXPECT_NEAR(v1.X(), cv1.X()-cv2.X(), point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Y()-cv2.Y(), point::precision);
        EXPECT_NEAR(v1.Length(), sqrt((cv1.X()-cv2.X())*(cv1.X()-cv2.X()) + (cv1.Y()-cv2.Y())*(cv1.Y()-cv2.Y())), point::precision);
        EXPECT_NEAR(v1.Phi(), (atan2(cv1.Y()-cv2.Y(),cv1.X()-cv2.X())>=0)?atan2(cv1.Y()-cv2.Y(),cv1.X()-cv2.X()):twopi+atan2(cv1.Y()-cv2.Y(),cv1.X()-cv2.X()),point::precision);

        v1 = cv1;
        v1 *= s;
        EXPECT_NEAR(v1.X(), cv1.X()*s, point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Y()*s, point::precision);
        EXPECT_NEAR(v1.Length(), cv1.Length()*std::abs(s), point::precision);
        EXPECT_NEAR(v1.Phi(), (atan2(cv1.Y()*s,cv1.X()*s)>=0)?atan2(cv1.Y()*s,cv1.X()*s):twopi+atan2(cv1.Y()*s,cv1.X()*s), point::precision);

        v1 = cv1;
        v1 /= s;
        EXPECT_NEAR(v1.X(), cv1.X()/s, point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Y()/s, point::precision);
        EXPECT_NEAR(v1.Length(), cv1.Length()/std::abs(s), point::precision);
        EXPECT_NEAR(v1.Phi(), (atan2(cv1.Y()/s,cv1.X()/s)>=0)?atan2(cv1.Y()/s,cv1.X()/s):twopi+atan2(cv1.Y()/s,cv1.X()/s), point::precision);

        double v3 = cv1 * cv2;
        EXPECT_NEAR(v3, cv1.X()*cv2.X()+cv1.Y()*cv2.Y(), point::precision);

        v1 = cv1 * s;
        EXPECT_NEAR(v1.X(), cv1.X()*s, point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Y()*s, point::precision);
        EXPECT_NEAR(v1.Length(), cv1.Length()*std::abs(s), point::precision);
        EXPECT_NEAR(v1.Phi(), (atan2(cv1.Y()*s,cv1.X()*s)>=0)?atan2(cv1.Y()*s,cv1.X()*s):twopi+atan2(cv1.Y()*s,cv1.X()*s), point::precision);

        v2 = s *cv2 ;
        EXPECT_NEAR(v2.X(), cv2.X()*s, point::precision);
        EXPECT_NEAR(v2.Y(), cv2.Y()*s, point::precision);
        EXPECT_NEAR(v2.Length(), cv2.Length()*std::abs(s), point::precision);
        EXPECT_NEAR(v2.Phi(), (atan2(cv2.Y()/s,cv2.X()/s)>=0)?atan2(cv2.Y()/s,cv2.X()/s):twopi+atan2(cv2.Y()/s,cv2.X()/s), point::precision);

        v1 = cv1 / s;
        EXPECT_NEAR(v1.X(), cv1.X()/s, point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Y()/s, point::precision);
        EXPECT_NEAR(v1.Length(), cv1.Length()/std::abs(s), point::precision);
        EXPECT_NEAR(v1.Phi(), (atan2(cv1.Y()/s,cv1.X()/s)>=0)?atan2(cv1.Y()/s,cv1.X()/s):twopi+atan2(cv1.Y()/s,cv1.X()/s), point::precision);
    }
}

#pragma endregion
#pragma region - Test vector3D

TEST(vector3D, ctor)
{
    vector3D v0 = vector3D();
    EXPECT_NEAR(v0.Phi(), 0.,point::precision);
    EXPECT_NEAR(v0.Theta(), acos(1.),point::precision);
    EXPECT_NEAR(v0.Length(),1.,point::precision);

    std::uniform_real_distribution<double> unid(-100.,100.);
    std::default_random_engine re;

    double twopi = 2. * acos(-1.);
    double pi = acos(-1.);

    for(size_t k=0; k<NTEST; k++)
    {
        double l = unid(re);
        double p = unid(re);
        double t = unid(re);

        double phi = std::fmod(p, twopi);
        if(phi < 0)
            phi += twopi;

        double theta = std::fmod(t, pi);
        if (theta < 0)
            theta += pi;
        
        double x=std::abs(l)*sin(theta)*cos(phi);
        double y=std::abs(l)*sin(theta)*sin(phi);
        double z=std::abs(l)*cos(theta);
        
        vector3D v3 = vector3D(point(3,x,y,z));
        EXPECT_NEAR(v3.Length(), sqrt(x*x + y*y + z*z), point::precision);
        EXPECT_NEAR(v3.Phi(), (atan2(y,x)>=0)?atan2(y,x):2*acos(-1)+atan2(y,x), point::precision);
        EXPECT_NEAR(v3.Theta(), theta, point::precision);
        EXPECT_NEAR(v3.X(), x, point::precision);
        EXPECT_NEAR(v3.Y(), y, point::precision);
        EXPECT_NEAR(v3.Z(), z, point::precision);

        vector3D v4 = vector3D(l,phi,theta);
        EXPECT_NEAR(v4.Length(), std::abs(l), point::precision);
        EXPECT_NEAR(v4.Phi()  , phi, point::precision);
        EXPECT_NEAR(v4.Theta(), theta, point::precision);
        EXPECT_NEAR(v4.X(), x, point::precision);
        EXPECT_NEAR(v4.Y(), y, point::precision);
        EXPECT_NEAR(v4.Z(), z, point::precision);
    }
}

TEST(vector3D,exception)
{
    std::uniform_real_distribution<double> unid(-100.,100.);
    std::default_random_engine re;

    double a_x= unid(re);
    double a_y= unid(re);
    double a_z= unid(re);

    point pa_1(1,a_x);
    point pa_2(2,a_x,a_y);
    point pa_3(3,a_x,a_y,a_z);

    EXPECT_ANY_THROW(vector3D(point(1,a_x)));
    EXPECT_ANY_THROW(vector3D(point(2,a_x,a_y)));

    EXPECT_ANY_THROW(vector3D(pa_1,pa_1));
    EXPECT_ANY_THROW(vector3D(pa_1,pa_2));
    EXPECT_ANY_THROW(vector3D(pa_1,pa_3));
    EXPECT_ANY_THROW(vector3D(pa_2,pa_1));
    EXPECT_ANY_THROW(vector3D(pa_2,pa_2));
    EXPECT_ANY_THROW(vector3D(pa_2,pa_3));
    EXPECT_ANY_THROW(vector3D(pa_3,pa_1));
    EXPECT_ANY_THROW(vector3D(pa_3,pa_2));
}

TEST(vector3D,operator)
{
    std::uniform_real_distribution<double> unidl(0,1000.);
    std::uniform_real_distribution<double> unidp(0,2.*acos(-1));
    std::uniform_real_distribution<double> unidt(0,acos(-1));
    std::default_random_engine re;

    for(size_t k=0; k<NTEST; k++)
    {
        double l1 = unidl(re);
        double p1 = unidp(re);
        double t1 = unidt(re);

        double l2 = unidl(re);
        double p2 = unidp(re);
        double t2 = unidt(re);
        
        vector3D v1 = vector3D(l1,p1,t1);
        vector3D v2 = vector3D(l2,p2,t2);
        vector3D v3 = vector3D(l1,p1,t2);

        EXPECT_TRUE(v1 == v1);
        EXPECT_TRUE(v1 >= v1);
        EXPECT_TRUE(v1 <= v1);
        EXPECT_TRUE(v2 == v2);
        EXPECT_TRUE(v2 >= v2);
        EXPECT_TRUE(v2 <= v2);

        if((l1 != l2) || (p1 != p2) || (t1 != t2))
        {
            EXPECT_FALSE(v1 == v2);
            EXPECT_FALSE(v1 == v2);
        }
        else
        {
            EXPECT_TRUE(v1 == v2);
            EXPECT_TRUE(v2 == v1);
        }

        if(std::abs(l1 - l2) < point::precision)
        {
            if((std::abs(p1 - p2) < point::precision) && (t1 > t2))
            {
                EXPECT_TRUE (v1 > v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_FALSE(v1 < v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_TRUE (v2 < v1)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_FALSE(v2 > v1)<<v1.Dump()<<" >  "<<v2.Dump();
            }
            else if ((std::abs(p1 - p2) < point::precision) & (t1 < t2))
            {
                EXPECT_FALSE(v1 > v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_TRUE (v1 < v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_FALSE(v2 < v1)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_TRUE (v2 > v1)<<v1.Dump()<<" >  "<<v2.Dump();
            }
            else if (p1 > p2)
            {
                EXPECT_TRUE (v1 > v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_FALSE(v1 < v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_TRUE (v2 < v1)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_FALSE(v2 > v1)<<v1.Dump()<<" >  "<<v2.Dump();
            }
            else
            {
                EXPECT_FALSE(v1 > v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_TRUE (v1 < v2)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_FALSE(v2 < v1)<<v1.Dump()<<" >  "<<v2.Dump();
                EXPECT_TRUE (v2 > v1)<<v1.Dump()<<" >  "<<v2.Dump();
            }
        }
        else if (l1 > l2)
        {
            EXPECT_TRUE (v1 > v2)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_FALSE(v1 < v2)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_TRUE (v2 < v1)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_FALSE(v2 > v1)<<v1.Dump()<<" >  "<<v2.Dump();
        }
        else
        {
            EXPECT_FALSE(v1 > v2)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_TRUE (v1 < v2)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_FALSE(v2 < v1)<<v1.Dump()<<" >  "<<v2.Dump();
            EXPECT_TRUE (v2 > v1)<<v1.Dump()<<" >  "<<v2.Dump();
        }        
            

        if(t1 >= t2)
        {
            EXPECT_TRUE (v1 > v3)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_FALSE(v1 < v3)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_TRUE (v3 < v1)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_FALSE(v3 > v1)<<v1.Dump()<<" >  "<<v3.Dump();
        }
        else
        {
            EXPECT_FALSE(v1 > v3)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_TRUE (v1 < v3)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_FALSE(v3 < v1)<<v1.Dump()<<" >  "<<v3.Dump();
            EXPECT_TRUE (v3 > v1)<<v1.Dump()<<" >  "<<v3.Dump();
        }
    }
}

TEST(vector3D,aritmetic)
{
    std::uniform_real_distribution<double> unids(-1000.,1000.);
    std::uniform_real_distribution<double> unidl(0,1000.);
    std::uniform_real_distribution<double> unidp(0,2.*acos(-1));
    std::uniform_real_distribution<double> unidt(0,acos(-1));
    std::default_random_engine re;

    for(size_t k=0; k<NTEST; k++)
    {
        double s  = unids(re);

        double l1 = unidl(re);
        double p1 = unidp(re);
        double t1 = unidt(re);

        double l2 = unidl(re);
        double p2 = unidp(re);
        double t2 = unidt(re);

        double l3 = unidl(re);
        double p3 = unidp(re);

        vector3D v1 = vector3D(l1,p1,t1);
        vector3D v2 = vector3D(l2,p2,t2);
        vector3D v3 = vector2D(l3,p3);

        vector3D cv1 = vector3D(v1);
        vector3D cv2 = vector3D(v2);
        vector2D cv3 = vector2D(v3);

        v1 += v2;
        EXPECT_NEAR(v1.X(), cv1.X()+cv2.X(), point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Y()+cv2.Y(), point::precision);
        EXPECT_NEAR(v1.Z(), cv1.Z()+cv2.Z(), point::precision);
        EXPECT_NEAR(v1.Length(), sqrt((cv1.X()+cv2.X())*(cv1.X()+cv2.X()) + (cv1.Y()+cv2.Y())*(cv1.Y()+cv2.Y()) + (cv1.Z()+cv2.Z())*(cv1.Z()+cv2.Z())), point::precision);
        EXPECT_NEAR(v1.Phi(), (atan2(cv1.Y()+cv2.Y(),cv1.X()+cv2.X())>=0)?atan2(cv1.Y()+cv2.Y(),cv1.X()+cv2.X()):2*acos(-1)+atan2(cv1.Y()+cv2.Y(),cv1.X()+cv2.X()), point::precision);
        EXPECT_NEAR(v1.Theta(), (acos((cv1.Z()+cv2.Z())/v1.Length())>=0)?acos((cv1.Z()+cv2.Z())/v1.Length()):acos(-1)+acos((cv1.Z()+cv2.Z())/v1.Length()), point::precision);

        v1 = cv1;
        v2 = cv2;
        v1 -= v2;
        EXPECT_NEAR(v1.X(), cv1.X()-cv2.X(), point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Y()-cv2.Y(), point::precision);
        EXPECT_NEAR(v1.Z(), cv1.Z()-cv2.Z(), point::precision);
        EXPECT_NEAR(v1.Length(), sqrt((cv1.X()-cv2.X())*(cv1.X()-cv2.X()) + (cv1.Y()-cv2.Y())*(cv1.Y()-cv2.Y()) + (cv1.Z()-cv2.Z())*(cv1.Z()-cv2.Z())), point::precision);
        EXPECT_NEAR(v1.Phi(), (atan2(cv1.Y()-cv2.Y(),cv1.X()-cv2.X())>=0)?atan2(cv1.Y()-cv2.Y(),cv1.X()-cv2.X()):2*acos(-1)+atan2(cv1.Y()-cv2.Y(),cv1.X()-cv2.X()), point::precision);
        EXPECT_NEAR(v1.Theta(), (acos((cv1.Z()-cv2.Z())/v1.Length())>=0)?acos((cv1.Z()-cv2.Z())/v1.Length()):acos(-1)+acos((cv1.Z()-cv2.Z())/v1.Length()), point::precision);

        v1 = cv1;
        v1 += v3;
        EXPECT_NEAR(v1.X(), cv1.X()+v3.X(), point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Y()+v3.Y(), point::precision);
        EXPECT_NEAR(v1.Z(), cv1.Z(), point::precision);
        EXPECT_NEAR(v1.Length(), sqrt((cv1.X()+v3.X())*(cv1.X()+v3.X()) + (cv1.Y()+v3.Y())*(cv1.Y()+v3.Y())+cv1.Z()*cv1.Z()), point::precision);
        EXPECT_NEAR(v1.Phi(), (atan2(cv1.Y()+v3.Y(),cv1.X()+v3.X())>=0)?atan2(cv1.Y()+v3.Y(),cv1.X()+v3.X()):2*acos(-1)+atan2(cv1.Y()+v3.Y(),cv1.X()+v3.X()), point::precision);
        EXPECT_NEAR(v1.Theta(), (acos(cv1.Z()/v1.Length())>=0)?acos(cv1.Z()/v1.Length()):acos(-1.)+acos(cv1.Z()/v1.Length()), point::precision);

        v1 = cv1;
        v1 -= v3;
        EXPECT_NEAR(v1.X(), cv1.X()-v3.X(), point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Y()-v3.Y(), point::precision);
        EXPECT_NEAR(v1.Z(), cv1.Z(), point::precision);
        EXPECT_NEAR(v1.Length(), sqrt((cv1.X()-v3.X())*(cv1.X()-v3.X()) + (cv1.Y()-v3.Y())*(cv1.Y()-v3.Y())+cv1.Z()*cv1.Z()), point::precision);
        EXPECT_NEAR(v1.Phi(), (atan2(cv1.Y()-v3.Y(),cv1.X()-v3.X())>=0)?atan2(cv1.Y()-v3.Y(),cv1.X()-v3.X()):2*acos(-1)+atan2(cv1.Y()-v3.Y(),cv1.X()-v3.X()), point::precision);
        EXPECT_NEAR(v1.Theta(), (acos(cv1.Z()/v1.Length())>=0)?acos(cv1.Z()/v1.Length()):acos(-1.)+acos(cv1.Z()/v1.Length()), point::precision);

        v1 = cv1;
        v1 *= s;
        EXPECT_NEAR(v1.Length(), cv1.Length()*std::abs(s), point::precision);
        EXPECT_NEAR(v1.Phi()   , (atan2(cv1.Y()*s,cv1.X()*s)>=0)?atan2(cv1.Y()*s,cv1.X()*s):2*acos(-1)+atan2(cv1.Y()*s,cv1.X()*s), point::precision);
        EXPECT_NEAR(v1.Theta() , (acos (cv1.Z()*s/(cv1.Length()*std::abs(s)))>=0)?acos(cv1.Z()*s/(cv1.Length()*std::abs(s))):acos(-1)+acos(cv1.Z()*s/(cv1.Length()*std::abs(s))), point::precision);

        v1 = cv1;
        v1 /= s;
        EXPECT_NEAR(v1.Length(), cv1.Length()/std::abs(s), point::precision);
        EXPECT_NEAR(v1.Phi()   , (atan2(cv1.Y()/s,cv1.X()/s)>=0)?atan2(cv1.Y()/s,cv1.X()/s):2*acos(-1)+atan2(cv1.Y()/s,cv1.X()/s), point::precision);
        EXPECT_NEAR(v1.Theta() , (acos (cv1.Z()/s/(cv1.Length()/std::abs(s)))>=0)?acos(cv1.Z()/s/(cv1.Length()/std::abs(s))):acos(-1)+acos(cv1.Z()*s/(cv1.Length()/std::abs(s))), point::precision);

        point::precision = 5e-9;

        v1 = cv1;
        v1 ^= v2;
        EXPECT_NEAR(v1.X(), cv1.Y()*cv2.Z()-cv1.Z()*cv2.Y(), point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Z()*cv2.X()-cv1.X()*cv2.Z(), point::precision);
        EXPECT_NEAR(v1.Z(), cv1.X()*cv2.Y()-cv1.Y()*cv2.X(), point::precision);

        v1 = cv1;
        v1 ^= v3;
        EXPECT_NEAR(v1.X(), -cv1.Z()*v3.Y(), point::precision);
        EXPECT_NEAR(v1.Y(), cv1.Z()*v3.X() , point::precision);
        EXPECT_NEAR(v1.Z(), cv1.X()*v3.Y()-cv1.Y()*v3.X(), point::precision);

        double val1;
        double val2;
        val1 = cv1 * cv2;
        EXPECT_NEAR(val1, cv1.X()*cv2.X()+cv1.Y()*cv2.Y()+cv1.Z()*cv2.Z(), point::precision);

        val1 = cv1 * cv3;
        EXPECT_NEAR(val1, cv1.X()*cv3.X()+cv1.Y()*cv3.Y(), point::precision);

        val2 = cv3 * cv1;
        EXPECT_NEAR(val2, val1, point::precision);

        v1 = cv1 * s;
        EXPECT_NEAR(v1.Length(), cv1.Length()*std::abs(s), point::precision);
        EXPECT_NEAR(v1.Phi()   , (atan2(cv1.Y()*s,cv1.X()*s)>=0)?atan2(cv1.Y()*s,cv1.X()*s):2*acos(-1)+atan2(cv1.Y()*s,cv1.X()*s), point::precision);
        EXPECT_NEAR(v1.Theta() , (acos (cv1.Z()*s/(cv1.Length()*std::abs(s)))>=0)?acos(cv1.Z()*s/(cv1.Length()*std::abs(s))):acos(-1)+acos(cv1.Z()*s/(cv1.Length()*std::abs(s))), point::precision);

        v2 = s *cv1 ;
        EXPECT_NEAR(v2.Length(), v1.Length(), point::precision);
        EXPECT_NEAR(v2.Phi()   , v1.Phi()   , point::precision);
        EXPECT_NEAR(v2.Theta() , v1.Theta() , point::precision);

        v1 = cv1 / s;
        EXPECT_NEAR(v1.Length(), cv1.Length()/std::abs(s), point::precision);
        EXPECT_NEAR(v1.Phi()   , (atan2(cv1.Y()/s,cv1.X()/s)>=0)?atan2(cv1.Y()/s,cv1.X()/s):2*acos(-1)+atan2(cv1.Y()/s,cv1.X()/s), point::precision);
        EXPECT_NEAR(v1.Theta() , (acos (cv1.Z()/s/(cv1.Length()/std::abs(s)))>=0)?acos(cv1.Z()/s/(cv1.Length()/std::abs(s))):acos(-1)+acos(cv1.Z()*s/(cv1.Length()/std::abs(s))), point::precision);

        vector3D cvl = v1^v2;
        EXPECT_NEAR(cvl.X(), v1.Y()*v2.Z()-v1.Z()*v2.Y(), point::precision);
        EXPECT_NEAR(cvl.Y(), v1.Z()*v2.X()-v1.X()*v2.Z(), point::precision);
        EXPECT_NEAR(cvl.Z(), v1.X()*v2.Y()-v1.Y()*v2.X(), point::precision);

        vector3D clv = v1^v3;
        EXPECT_NEAR(clv.X(), v1.Y()*0.-v1.Z()*v3.Y(), point::precision);
        EXPECT_NEAR(clv.Y(), v1.Z()*v3.X()-v1.X()*0., point::precision);
        EXPECT_NEAR(clv.Z(), v1.X()*v3.Y()-v1.Y()*v3.X(), point::precision);

        vector3D cvv = v1^v3;
        EXPECT_NEAR(cvv.X(), v1.Y()*0.-v1.Z()*v3.Y(), point::precision);
        EXPECT_NEAR(cvv.Y(), v1.Z()*v3.X()-v1.X()*0., point::precision);
        EXPECT_NEAR(cvv.Z(), v1.X()*v3.Y()-v1.Y()*v3.X(), point::precision);

        point::precision = 1e-9;
    }

}

#pragma endregion

