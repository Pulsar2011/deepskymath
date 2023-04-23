#include "gTest_version.h"
#include <gtest/gtest.h>
#include <DSTmath/DSM_version.h>

namespace DST
{
    namespace Math
    {
    class version_Tester: public gGIT
    {
        protected:
        public:
            static std::vector<int16_t> version_to_test() {return DST::Math::gGIT::this_ver->getVersion();}
    };

    TEST(version_Tester, test_getVersion)
    {
        gGIT::version();
        std::vector<int16_t> ver = DST::Math::gGIT::this_ver->getVersion();

        ASSERT_EQ(ver.size(), 3);
        ASSERT_EQ(DST::Math::gGIT::isEqual(std::string("v")+std::to_string(ver[0])+std::string(".")+std::to_string(ver[1]) +std::string(".")+std::to_string(ver[2])+std::string("-")+this_VERSION_SHA1), true)<<std::string("v")+std::to_string(ver[0])+std::string(".")+std::to_string(ver[1]) +std::string(".")+std::to_string(ver[2])+std::string("-")+this_VERSION_SHA1<<"%"<<DST::Math::gGIT::this_ver->version;
    }

    TEST(version_Tester, test_readVersion)
    {
        gGIT::version();
        std::vector<int16_t> ver = DST::Math::gGIT::this_ver->readVersion(this_VERSION_SHORT);

        ASSERT_EQ(ver.size(), 3);
        ASSERT_EQ(DST::Math::gGIT::isEqual(std::string("v")+std::to_string(ver[0])+std::string(".")+std::to_string(ver[1]) +std::string(".")+std::to_string(ver[2])+std::string("-")+this_VERSION_SHA1), true);
    }

    TEST(version_Tester, test_private)
    {
        DST::Math::gGIT::version();
        std::vector<int16_t> ver0 = DST::Math::gGIT::this_ver->getVersion();
        std::vector<int16_t> ver1 = DST::Math::gGIT::this_ver->readVersion(this_VERSION_SHORT);

        ASSERT_EQ(ver1.size(), ver0.size());
        ASSERT_EQ(ver0[0], ver1[0]);
        ASSERT_EQ(ver0[1], ver1[1]);
        ASSERT_EQ(ver0[2], ver1[2]);
    }

    TEST(version_Tester, test_operator_newer)
    {
        DST::Math::gGIT::version();
        std::vector<int16_t> ver = DST::Math::gGIT::this_ver->readVersion(this_VERSION_SHORT);
        std::vector<int16_t> to_test(ver.size());

        for(size_t k=0; k<ver.size()+2; k++)
        {
            for(size_t i =0; i < ver.size(); i++)
                to_test[i] = ver[i];

            if(k>0 && k<=to_test.size())
                to_test[to_test.size()-k] += 1;
            else if (k > to_test.size())
            {
                to_test[0] +=1;
                to_test[1] +=1;
                to_test[2] +=1;
            }

            std::string newer = std::string("v")+std::to_string(to_test[0])
                              + std::string(".")+std::to_string(to_test[1])
                              + std::string(".")+std::to_string(to_test[2]);            
            
            if(to_test[0] == ver[0] &&
               to_test[1] == ver[1] &&
               to_test[2] == ver[2] )
            {
                ASSERT_EQ((*DST::Math::gGIT::this_ver)< newer,false)<<DST::Math::gGIT::version_short()<<" < "<<newer<<" ?!" <<(ver[0]<to_test[0] )<<"||"<<(ver[1]<to_test[1])<<"||"<<(ver[2]<to_test[2]);
                ASSERT_EQ((*DST::Math::gGIT::this_ver)<=newer,true )<<DST::Math::gGIT::version_short()<<" <= "<<newer<<" ?!"<<(ver[0]<=to_test[0])<<"||"<<(ver[1]<=to_test[1])<<"||"<<(ver[2]<=to_test[2]);

                ASSERT_EQ((*DST::Math::gGIT::this_ver)> newer,false)<<DST::Math::gGIT::version_short()<<" > "<<newer<<" ?!" <<(ver[0]>to_test[0] )<<"||"<<(ver[1]>to_test[1] )<<"||"<<(ver[2]>to_test[2] );
                ASSERT_EQ((*DST::Math::gGIT::this_ver)>=newer,true )<<DST::Math::gGIT::version_short()<<" >= "<<newer<<" ?!"<<(ver[0]>=to_test[0])<<"||"<<(ver[1]>=to_test[1])<<"||"<<(ver[2]>=to_test[2]);
            }
            else
            {
                ASSERT_EQ((*DST::Math::gGIT::this_ver)< newer,true)<<DST::Math::gGIT::version_short()<<" <  "<<newer<<" ?!"<<(ver[0]< to_test[0])<<"||"<<(ver[1]< to_test[1])<<"||"<<(ver[2]< to_test[2]);
                ASSERT_EQ((*DST::Math::gGIT::this_ver)<=newer,true)<<DST::Math::gGIT::version_short()<<" <= "<<newer<<" ?!"<<(ver[0]<=to_test[0])<<"||"<<(ver[1]<=to_test[1])<<"||"<<(ver[2]<=to_test[2]);

                ASSERT_EQ((*DST::Math::gGIT::this_ver)> newer,false)<<DST::Math::gGIT::version_short()<<" > "<<newer<<" ?!"<<(ver[0]> to_test[0])<<"||"<<(ver[1]> to_test[1])<<"||"<<(ver[2]> to_test[2]);
                ASSERT_EQ((*DST::Math::gGIT::this_ver)>=newer,false)<<DST::Math::gGIT::version_short()<<" >="<<newer<<" ?!"<<(ver[0]>=to_test[0])<<"||"<<(ver[1]>=to_test[1])<<"||"<<(ver[2]>=to_test[2]);
            }
        }
    }
    TEST(version_Tester, test_operator_older)
    {
        DST::Math::gGIT::version();
        std::vector<int16_t> ver = DST::Math::gGIT::this_ver->readVersion(this_VERSION_SHORT);
        std::vector<int16_t> to_test(ver.size());

        for(size_t k=0; k<ver.size()+2; k++)
        {
            for(size_t i =0; i < ver.size(); i++)
                to_test[i] = ver[i];
            
            if( k > 0 && k <= to_test.size())
                to_test[to_test.size()-k] -= ((to_test[to_test.size()-k]>0)?1:0);
            else if (k > to_test.size())
            {
                to_test[0] -= ((to_test[0]>0)?1:0);
                to_test[1] -= ((to_test[1]>0)?1:0); if(to_test[0] < 1) to_test[1] ++;
                to_test[2] -= ((to_test[2]>0)?1:0); if(to_test[0] < 1) to_test[2] ++;
            }

            std::string older = std::string("v")+std::to_string(to_test[0])
                              + std::string(".")+std::to_string(to_test[1])
                              + std::string(".")+std::to_string(to_test[2]);
            
            if(to_test[0] == ver[0] &&
               to_test[1] == ver[1] &&
               to_test[2] == ver[2] )
            {
                ASSERT_EQ((*DST::Math::gGIT::this_ver)<older ,false) <<DST::Math::gGIT::version_short()<<" < "<<older<<" ?!" <<(ver[0]<to_test[0] )<<"||"<<(ver[1]<to_test[1])<<"||"<<(ver[2]<to_test[2]);
                ASSERT_EQ((*DST::Math::gGIT::this_ver)<=older,true ) <<DST::Math::gGIT::version_short()<<" <= "<<older<<" ?!"<<(ver[0]<=to_test[0])<<"||"<<(ver[1]<=to_test[1])<<"||"<<(ver[2]<=to_test[2]);

                ASSERT_EQ((*DST::Math::gGIT::this_ver)>older ,false) <<DST::Math::gGIT::version_short()<<" > "<<older<<" ?!" <<(ver[0]>to_test[0] )<<"||"<<(ver[1]>to_test[1] )<<"||"<<(ver[2]>to_test[2] );
                ASSERT_EQ((*DST::Math::gGIT::this_ver)>=older,true ) <<DST::Math::gGIT::version_short()<<" >= "<<older<<" ?!"<<(ver[0]>=to_test[0])<<"||"<<(ver[1]>=to_test[1])<<"||"<<(ver[2]>=to_test[2]);
            }
            else if(to_test[0] == ver[0])
            {
                if(ver[1]<to_test[1])
                    ASSERT_EQ((*DST::Math::gGIT::this_ver)<older ,true) <<DST::Math::gGIT::version_short()<<" < "<<older<<" ?!" <<(ver[0]<to_test[0] )<<"||"<<(ver[1]<to_test[1])<<"||"<<(ver[2]<to_test[2]);
                else
                    ASSERT_EQ((*DST::Math::gGIT::this_ver)<older ,false) <<DST::Math::gGIT::version_short()<<" < "<<older<<" ?!" <<(ver[0]<to_test[0] )<<"||"<<(ver[1]<to_test[1])<<"||"<<(ver[2]<to_test[2]);
                    
                if(ver[2]<to_test[2])
                    ASSERT_EQ((*DST::Math::gGIT::this_ver)<older ,true) <<DST::Math::gGIT::version_short()<<" < "<<older<<" ?!" <<(ver[0]<to_test[0] )<<"||"<<(ver[1]<to_test[1])<<"||"<<(ver[2]<to_test[2]);
                else
                    ASSERT_EQ((*DST::Math::gGIT::this_ver)<older ,false) <<DST::Math::gGIT::version_short()<<" < "<<older<<" ?!" <<(ver[0]<to_test[0] )<<"||"<<(ver[1]<to_test[1])<<"||"<<(ver[2]<to_test[2]);
            }
            else
            {
                ASSERT_EQ((*DST::Math::gGIT::this_ver)<older, false)<<DST::Math::gGIT::version_short()<<" < "<<older<<" ?!" <<(ver[0]<to_test[0] )<<"||"<<(ver[1]<to_test[1] )<<"||"<<(ver[2]<to_test[2])<<" -> "<<((*gGIT::this_ver)<older);
                ASSERT_EQ((*DST::Math::gGIT::this_ver)<=older,false)<<DST::Math::gGIT::version_short()<<" <= "<<older<<" ?!"<<(ver[0]<=to_test[0])<<"||"<<(ver[1]<=to_test[1])<<"||"<<(ver[2]<=to_test[2])<<" -> "<<((*DST::Math::gGIT::this_ver)<=older);

                ASSERT_EQ((*DST::Math::gGIT::this_ver)>older, true)<<DST::Math::gGIT::version_short()<<" > "<<older<<" ?!" <<(ver[0]>to_test[0] )<<"||"<<(ver[1]>to_test[1] )<<"||"<<(ver[2]>to_test[2])<<" -> "<<((*DST::Math::gGIT::this_ver)>older);
                ASSERT_EQ((*DST::Math::gGIT::this_ver)>=older,true)<<DST::Math::gGIT::version_short()<<" >= "<<older<<" ?!"<<(ver[0]>=to_test[0])<<"||"<<(ver[1]>=to_test[1])<<"||"<<(ver[2]>=to_test[2])<<" -> "<<((*DST::Math::gGIT::this_ver)>=older);
            }
        }
    }
}
}

using namespace DST::Math;

TEST(test_version, version_equal)
{

    std::string ver=DST::Math::gGIT::version();
    ASSERT_EQ(!ver.compare(this_VERSION), true);
    ASSERT_EQ(DST::Math::gGIT::isEqual(ver), true);
    ASSERT_EQ(DST::Math::gGIT::isEqual(this_VERSION), true);

    ASSERT_EQ(DST::Math::gGIT::isNotEqual(ver), false);
    ASSERT_EQ(DST::Math::gGIT::isNotEqual(this_VERSION), false);

    ASSERT_EQ(DST::Math::gGIT::isEqual(this_VERSION_SHORT), false);
    ASSERT_EQ(DST::Math::gGIT::isEqual(this_VERSION_SHA1), false);

    ASSERT_EQ(DST::Math::gGIT::isNotEqual(this_VERSION_SHORT), true);
    ASSERT_EQ(DST::Math::gGIT::isNotEqual(this_VERSION_SHA1), true);
}

TEST(test_version, acessor)
{
    std::string ver = DST::Math::gGIT::version();
    ASSERT_EQ(!ver.compare(this_VERSION), true);

    std::string ver_short = DST::Math::gGIT::version_short();
    ASSERT_EQ(!ver_short.compare(this_VERSION_SHORT), true);

    std::string sha1 = DST::Math::gGIT::sha1();
    ASSERT_EQ(!sha1.compare(this_VERSION_SHA1), true);

    ASSERT_EQ(DST::Math::gGIT::checkVersion(this_VERSION), true);
    ASSERT_EQ(DST::Math::gGIT::checkVersionShort(this_VERSION_SHORT), true);
    ASSERT_EQ(DST::Math::gGIT::checkSHA1(this_VERSION_SHA1), true);
}

TEST(test_version, version_newer)
{
    ASSERT_EQ(DST::Math::gGIT::isNewerThan(this_VERSION_SHORT), true) <<DST::Math::gGIT::version_short()  <<" newer than "<<this_VERSION_SHORT<<" ? "<<gGIT::isNewerThan(this_VERSION_SHORT);
    ASSERT_EQ(DST::Math::gGIT::isOlderThan(this_VERSION_SHORT), false)<<DST::Math::gGIT::version_short()  <<" older than "<<this_VERSION_SHORT<<" ? "<<gGIT::isOlderThan(this_VERSION_SHORT);

    std::vector<int16_t> ver = version_Tester::version_to_test();
    ASSERT_EQ(ver.size(),3);

    for(int k=0; k<3; k++)
    {
        std::vector<int16_t> to_test = ver;
        to_test[k] = to_test[k]+1;

        std::string newer = std::string("v")+std::to_string(to_test[0])
                          + std::string(".")+std::to_string(to_test[1])
                          + std::string(".")+std::to_string(to_test[2]);
        ASSERT_EQ(DST::Math::gGIT::isOlderThan(newer),true)  <<DST::Math::gGIT::version_short()  <<" older than "<<newer<<" ? "<<DST::Math::gGIT::isOlderThan(newer);
        ASSERT_EQ(DST::Math::gGIT::isNewerThan(newer),false) <<DST::Math::gGIT::version_short()<<" newer than "<<newer<<" ? "  <<DST::Math::gGIT::isNewerThan(newer);
    }

    for(int k=0; k<3; k++)
    {
        std::vector<int16_t> to_test = ver;
        to_test[k] = to_test[k]-((to_test[k]>0)?1:0);

        std::string older = std::string("v")+std::to_string(to_test[0])
                          + std::string(".")+std::to_string(to_test[1])
                          + std::string(".")+std::to_string(to_test[2]);
        ASSERT_EQ(DST::Math::gGIT::isOlderThan(older),false) <<DST::Math::gGIT::version_short()<<" older than "<<older<<" ? "<<DST::Math::gGIT::isOlderThan(older);;
        ASSERT_EQ(DST::Math::gGIT::isNewerThan(older),true)  <<DST::Math::gGIT::version_short()<<" newer than "<<older<<" ? "<<DST::Math::gGIT::isNewerThan(older);;
    }
}