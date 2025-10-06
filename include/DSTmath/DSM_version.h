#ifndef DST_VERSION_H
#define DST_VERSION_H

#include <string>
#include <vector>
#include <stdexcept>
#include <string>

namespace DSM
{

  /*! 
  @class gGIT
  @brief Class to handle the version of the library
  @author GILLARD Willam
  @date 2024
  @licence CC-BY-NC 4.0
  
  This class is a singleton object which is used to handle the version of the library. It is used to check if the version of the library is compatible with the version of the application or library it is link to.
  Versionning is done using the git tag system. The version is stored in the following format: "vX.Y.Z-<number of commits since last tag>-g<short sha1 of the commit>". For example: "v1.2.3-4-g1234567". The version can be compared using the operators <, <=, >, >=, ==, !=. The comparison is done using the semantic versioning system. The number of commits since last tag and the sha1 are ignored in the comparison.
  */
  class gGIT
  {

#ifdef _gTEST_version_
  friend class version_Tester;
  FRIEND_TEST(version_Tester, test_getVersion);
  FRIEND_TEST(version_Tester, test_readVersion);
  FRIEND_TEST(version_Tester, test_private);
  FRIEND_TEST(version_Tester, test_operator_newer);
  FRIEND_TEST(version_Tester, test_operator_older);
#endif

  private:
   const std::string _VERSION;      //Full version string
   const std::string _VERSION_SHORT;//Short version string (vX.Y.Z)
   const std::string _VERSION_SHA1; //<number of commits since last tag>-g<short sha1 of the commit>
  
   ~gGIT(); //Destructor

   /*!
    @brief Default Constructor
    Initializes the version strings with the values defined during the build process.
   */
   gGIT();  //Constructor
   
   static gGIT *this_ver;

   /*!
    @brief Get the version as a vector of integers
    @return A vector of integers representing the version [X, Y, Z]
   */
   std::vector<int16_t> getVersion() const;

   /*!
    @brief Read a version string and convert it to a vector of integers
    @param ver The version string to read
    @return A vector of integers representing the version [X, Y, Z]
    @throw std::runtime_error if the version string is not valid
   */
   std::vector<int16_t> readVersion(const std::string& ver) const;

  public:
    static const std::string version();
    static const std::string version_short();
    static const std::string sha1();

    static bool checkVersion(const std::string& ver) {if(this_ver == NULL)gGIT::this_ver = new gGIT();return (ver.compare(this_ver->_VERSION)==0);}
    static bool checkVersionShort(const std::string& ver) {if(this_ver == NULL)gGIT::this_ver = new gGIT();return (ver.compare(this_ver->_VERSION_SHORT)==0);}
    static bool checkSHA1(const std::string& ver) {if(this_ver == NULL)gGIT::this_ver = new gGIT();return (ver.compare(this_ver->_VERSION_SHA1)==0);}

    bool operator<(const std::string& ver);
    bool operator<=(const std::string& ver);
    bool operator>(const std::string& ver);
    bool operator>=(const std::string& ver);
    bool operator==(const std::string& ver);
    bool operator!=(const std::string& ver);

    static bool isEqual(const std::string& ver)     {if(this_ver == NULL) gGIT::this_ver  = new gGIT(); return (*gGIT::this_ver)==ver;};
    static bool isNotEqual(const std::string& ver)  {if(this_ver == NULL) gGIT::this_ver  = new gGIT(); return (*gGIT::this_ver)!=ver;};
    static bool isOlderThan(const std::string& ver) {if(this_ver == NULL) gGIT::this_ver  = new gGIT(); return (*gGIT::this_ver)< ver;};
    static bool isNewerThan(const std::string& ver) {if(this_ver == NULL) gGIT::this_ver  = new gGIT(); return (*gGIT::this_ver)>=ver;};  
  
  };

}

#endif //DST_VERSION_H
