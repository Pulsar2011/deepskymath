#include <gtest/gtest.h>
#include <DSTmath/DSTmath.h>
#include <DSTmath/DSTgeo.h>
#include <random>
#include <functional>

using namespace DST::Math;

TEST(geometry,point_ctor)
{
    point p0 = point();
    EXPECT_EQ(p0.size(),0);

    for(size_t n = 0; n<4; n++)
    {
        point pp = point(n);
        EXPECT_EQ(pp.size(),n);
        for(size_t i=0; i<pp.size(); i++)
             EXPECT_EQ(pp[i],0.);
    }

    point p1 = point(1,1.);          EXPECT_EQ(p1.size(),1); for(size_t i=0; i<p1.size(); i++) EXPECT_DOUBLE_EQ(p1[i],static_cast<double>(i+1)); 
    point p2 = point(2,1.,2.);       EXPECT_EQ(p2.size(),2); for(size_t i=0; i<p2.size(); i++) EXPECT_DOUBLE_EQ(p2[i],static_cast<double>(i+1)); 
    point p3 = point(3,1.,2.,3.);    EXPECT_EQ(p3.size(),3); for(size_t i=0; i<p3.size(); i++) EXPECT_DOUBLE_EQ(p3[i],static_cast<double>(i+1)); 
    point p4 = point(4,1.,2.,3.,4.); EXPECT_EQ(p4.size(),4); for(size_t i=0; i<p4.size(); i++) EXPECT_DOUBLE_EQ(p4[i],static_cast<double>(i+1)); 

    point cpy = point(p4); 
    EXPECT_EQ(cpy.size(),4);
    for(size_t i=0; i<cpy.size(); i++) EXPECT_DOUBLE_EQ(cpy[i],p4[i]);

    EXPECT_EQ(point::precision,std::numeric_limits<double>::min()*1000.);
    EXPECT_FALSE(point::debug);
}

TEST(geometry,point_modifier)
{
    point pp = point(4,0,0,0,0);
    point cpy= point(pp);
    EXPECT_EQ(pp.size(),4);
    
    for(size_t i=0; i<pp.size(); i++)
        EXPECT_DOUBLE_EQ(pp[i],0.); 

    double lower_bound = -10;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    for(size_t k = 0; k < 10000; k++)
    {
        double x,y,z,t;
        x=unif(re);
        y=unif(re);
        z=unif(re);
        t=unif(re);

        pp.SetCoordinates(x,y,z,t);

        EXPECT_DOUBLE_EQ(pp[0],x);
        EXPECT_DOUBLE_EQ(pp[1],y);
        EXPECT_DOUBLE_EQ(pp[2],z);
        EXPECT_DOUBLE_EQ(pp[3],t);

        pp.SetCoordinate(0,z);
        pp.SetCoordinate(1,t);
        pp.SetCoordinate(2,x);
        pp.SetCoordinate(3,y);

        EXPECT_DOUBLE_EQ(pp[0],z);
        EXPECT_DOUBLE_EQ(pp[1],t);
        EXPECT_DOUBLE_EQ(pp[2],x);
        EXPECT_DOUBLE_EQ(pp[3],y);

        cpy.SetCoordinates(pp);
        EXPECT_DOUBLE_EQ(cpy[0],pp[0]);
        EXPECT_DOUBLE_EQ(cpy[1],pp[1]);
        EXPECT_DOUBLE_EQ(cpy[2],pp[2]);
        EXPECT_DOUBLE_EQ(cpy[3],pp[3]);

        cpy.SetX(x);
        cpy.SetY(y);
        cpy.SetZ(z);
        EXPECT_DOUBLE_EQ(cpy[0],x);
        EXPECT_DOUBLE_EQ(cpy[1],y);
        EXPECT_DOUBLE_EQ(cpy[2],z);
        EXPECT_DOUBLE_EQ(cpy[3],pp[3]);
    }
}

TEST(geometry,point_assert)
{
    point pp = point(4,0,0,0,0);
    point p1 = point(1,0);

    EXPECT_THROW(pp[4],std::out_of_range);
    EXPECT_THROW(pp[5],std::out_of_range);

    EXPECT_THROW(p1.Y(),std::out_of_range);
    EXPECT_THROW(p1.Z(),std::out_of_range);

    EXPECT_THROW(pp.SetCoordinate(4,2.),std::out_of_range);
    EXPECT_THROW(pp.SetCoordinate(5,2.),std::out_of_range);

    EXPECT_THROW(pp.SetCoordinates(p1),std::runtime_error);
    EXPECT_THROW(p1.SetCoordinates(pp),std::runtime_error);

    pp.Dump();
}

TEST(geometry,point_accessor)
{
    point p3 = point(3,0,0,0);
    point p2 = point(2,0,0);
    point p1 = point(1,0);

    double lower_bound = -10;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    for(size_t k = 0; k < 10000; k++)
    {
        double x,y,z;
        x=unif(re);
        y=unif(re);
        z=unif(re);

        double Phi   = (atan2(y,x) > 0) ? atan2(y,x) : 2.*DST::Math::MathCore::Pi() - atan2(y,x);
        double Theta = acos(z/sqrt(x*x+y*y+z*z));

        p3.SetCoordinates(x,y,z);
        p2.SetCoordinates(x,y,z);
        p1.SetCoordinates(x,y,z);
        const std::vector<double>& v3 = p3.GetPoint();
        const std::vector<double>& v2 = p2.GetPoint();
        const std::vector<double>& v1 = p1.GetPoint();

        EXPECT_EQ(v3.size(),p3.size());
        EXPECT_EQ(v2.size(),p2.size());
        EXPECT_EQ(v1.size(),p1.size());

        for(size_t i=0; i<v3.size(); i++)
            EXPECT_DOUBLE_EQ(v3[i],p3[i]);

        for(size_t i=0; i<v2.size(); i++)
            EXPECT_DOUBLE_EQ(v2[i],p2[i]);

        for(size_t i=0; i<v1.size(); i++)
            EXPECT_DOUBLE_EQ(v1[i],p1[i]);

        EXPECT_DOUBLE_EQ(p3.X(),x);
        EXPECT_DOUBLE_EQ(p3.Y(),y);
        EXPECT_DOUBLE_EQ(p3.Z(),z);

        EXPECT_DOUBLE_EQ(p2.X(),x);
        EXPECT_DOUBLE_EQ(p2.Y(),y);
        EXPECT_THROW    (p2.Z(),std::out_of_range);

        EXPECT_DOUBLE_EQ(p1.X(),x);
        EXPECT_THROW    (p1.Y(),std::out_of_range);
        EXPECT_THROW    (p1.Z(),std::out_of_range);

        EXPECT_NEAR(p3.R(),sqrt(x*x+y*y+z*z),1e-10);
        EXPECT_NEAR(p2.R(),sqrt(x*x+y*y),1e-10);
        EXPECT_NEAR(p1.R(),std::abs(x),1e-10);

        EXPECT_NEAR(p3.Phi(),Phi,1e-10);
        EXPECT_NEAR(p2.Phi(),Phi,1e-10);
        EXPECT_NEAR(p1.Phi(),0.,1e-10);

        EXPECT_NEAR(p3.Theta(),Theta,1e-10);
        EXPECT_NEAR(p2.Theta(),0.,1e-10);
        EXPECT_NEAR(p1.Theta(),0.,1e-10);
    }
}

TEST(geometry,point_math_operator)
{
    point pref  = point(3,0,0,0);
    point ptest;

    double lower_bound = -10;
    double upper_bound =  10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    point::precision=1e-10;

    for(size_t k = 0; k < 10000; k++)
    {
        double x0,y0,z0,x1,y1,z1;
        x0=unif(re);
        y0=unif(re);
        z0=unif(re);

        pref.SetCoordinates(x0,y0,z0);
        ptest = pref;

        EXPECT_DOUBLE_EQ(pref[0],ptest[0]);
        EXPECT_DOUBLE_EQ(pref[1],ptest[1]);
        EXPECT_DOUBLE_EQ(pref[2],ptest[2]);

        x1=unif(re);
        y1=unif(re);
        z1=unif(re);

        ptest.SetCoordinates(x1,y1,z1);
        EXPECT_DOUBLE_EQ(ptest[0],x1);
        EXPECT_DOUBLE_EQ(ptest[1],y1);
        EXPECT_DOUBLE_EQ(ptest[2],z1);

        ptest += pref;
        EXPECT_NEAR(ptest[0],x1+x0,point::precision);
        EXPECT_NEAR(ptest[1],y1+y0,point::precision);
        EXPECT_NEAR(ptest[2],z1+z0,point::precision);

        ptest -= pref;
        EXPECT_NEAR(ptest[0],x1,point::precision);
        EXPECT_NEAR(ptest[1],y1,point::precision);
        EXPECT_NEAR(ptest[2],z1,point::precision);

        ptest *= pref;
        EXPECT_NEAR(ptest[0],x1*x0,point::precision);
        EXPECT_NEAR(ptest[1],y1*y0,point::precision);
        EXPECT_NEAR(ptest[2],z1*z0,point::precision);

        ptest /= pref;
        EXPECT_NEAR(ptest[0],x1,point::precision);
        EXPECT_NEAR(ptest[1],y1,point::precision);
        EXPECT_NEAR(ptest[2],z1,point::precision);

        double scale = unif(re);
        double offset= unif(re);

        ptest += offset;
        EXPECT_NEAR(ptest[0],x1+offset,point::precision);
        EXPECT_NEAR(ptest[1],y1+offset,point::precision);
        EXPECT_NEAR(ptest[2],z1+offset,point::precision);

        ptest *= scale;
        EXPECT_NEAR(ptest[0],(x1+offset)*scale,point::precision);
        EXPECT_NEAR(ptest[1],(y1+offset)*scale,point::precision);
        EXPECT_NEAR(ptest[2],(z1+offset)*scale,point::precision);

        ptest /= scale;
        EXPECT_NEAR(ptest[0],x1+offset,point::precision);
        EXPECT_NEAR(ptest[1],y1+offset,point::precision);
        EXPECT_NEAR(ptest[2],z1+offset,point::precision);

        ptest -= offset;
        EXPECT_NEAR(ptest[0],x1,point::precision);
        EXPECT_NEAR(ptest[1],y1,point::precision);
        EXPECT_NEAR(ptest[2],z1,point::precision);
    }
}