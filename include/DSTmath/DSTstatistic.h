//
//  pixtemplate.h
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under CC BY-NC 4.0
//  You may share and adapt this code with attribution, 
//  but not for commercial purposes.
//  Licence text: https://creativecommons.org/licenses/by-nc/4.0/

#ifndef __DST_MATH_template__
#define __DST_MATH_template__

#include <cstdio>
#include <map>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <stdexcept>
#include <string>
#include <random>

typedef std::map<int,double> value; ///< ! Probability distribution function : a set of variable x with their associated probability
typedef std::map<int,double> weight;
typedef std::vector<double> farray;

typedef std::map<long,double> lvalue; ///< ! Probability distribution function : a set of variable x with their associated probability
typedef std::map<long,double> lweight;

namespace DST
{
    namespace Math
    {
#pragma region - template helper

        template<typename U>
        constexpr bool is_allowed_math_type =
            std::is_same_v<U,uint8_t> || std::is_same_v<U,int8_t> ||
            std::is_same_v<U,uint16_t> || std::is_same_v<U,int16_t> ||
            std::is_same_v<U,uint32_t> || std::is_same_v<U,int32_t> ||
            std::is_same_v<U,uint64_t> || std::is_same_v<U,int64_t> ||
            std::is_same_v<U,size_t>  ||
            std::is_same_v<U,float>   || std::is_same_v<U,double> ;

#pragma endregion

#pragma region - defined pdf
        /**
         * @brief Sample a normale distrbution of mean \c m_ and standard deviation \c s_
         * 
         * @param m_ Mean value of the normale distribution
         * @param s_ Standard deviation of the normale distribution
         * @return double value sampled from the normale distribution
         */
        double sample_normal_distribution ( double m_, double s_ );

        template<typename T>
        T sample_poisson_distribution(double l_)
        {
            if(l_ < 0.0)
                throw std::invalid_argument("[DST::Math::sample_poisson_distribution] l_ must be >= 0\033[0m");

            if constexpr (std::is_floating_point_v<T>)
                throw std::invalid_argument("[DST::Math::sample_poisson_distribution] Floating point datatype are not supported for poisson distribution\033[0m");

            static thread_local std::mt19937_64 rng(static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count()));

            std::poisson_distribution<int> dist(l_);
            return dist(rng);
        }

        /**
         * @brief Sample a poisson distrbution of mean \c l_
         * 
         * @param l_ Mean value of the poisson distribution
         * @return double value sampled from the poisson distribution
         */
        double sample_poisson_distribution ( double l_ );

        
        /**
         * @brief Sample a uniform distrbution in range [ \c i_ , \c f_ ]
         * 
         * @param i_ lower bound of the uniform distribution
         * @param f_ upper bound of the uniform distribution
         * @return a random number sampled from the uniform distribution
         */
        template<typename T>
        T sample_uniform_distribution ( const T& i_, const T& f_ )
        {
            if(!is_allowed_math_type<T>)
                throw std::invalid_argument("\033[1;35;47m[DST::Math::get_uniform_distribution] Unsuported datatype or data structure\033[0m");

            if( i_ >= f_ )
                throw std::invalid_argument("\033[1;35;47m[DST::Math::get_uniform_distribution] Invalid range for uniform distribution\033[0m");

            T rndNum;
            static int seed = 0;
        
            static std::mt19937_64 rndgen;
        
            if(seed == 0)
            {
                seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                rndgen = std::mt19937_64(seed);
            }
            
            if constexpr (std::is_floating_point_v<T>)
            {
                std::uniform_real_distribution<T> dis = std::uniform_real_distribution<T>(i_, f_);
                rndNum = dis(rndgen);
                return rndNum;
            }
            else if constexpr (std::is_integral_v<T>)
            {
                std::uniform_int_distribution<T> dis = std::uniform_int_distribution<T>(i_, f_);
                rndNum = dis(rndgen);
                return rndNum;
            }
            else
                throw std::invalid_argument("\033[1;35;47m[DST::Math::get_uniform_distribution] Unrecognized numerical datatype\033[0m");
        
            return 0;
        }
        
        /**
         * @brief Sample a uniform distrbution between 0 and \c f_
         * 
         * @param f_ upper bound of the uniform distribution
         * @return a random number sampled from the uniform distribution
         */
        template<typename T>
        T sample_uniform_distribution ( const T& f_ )
        {
            return sample_uniform_distribution<T>( 0, f_ );
        }

#pragma endregion

#pragma region - binned_pdf class definition
    /**
     *  @class binned_pdf DSTstatistic.h "DSTmath/DSTstatistic.h"
     *  @brief Store statistical probability distribution function
     *  @details The class binned_pdf store a statistical probability distribution function (PDF) defined by a set of values and their associated probability. The PDF is used to generate random numbers according to the PDF using the MCMC slice sampling method.
     */
    class binned_pdf
    {
    private:

#pragma region -- private member
#pragma region - variable
        value  fval;
        weight fpdf;
        
        double xSlice, ySlice;
        double fmax;                ///< The maximum of the pdf
        double lowEdge, upEdge;     ///< Low edge and upper edge of the random variable
        double lowSlice, upSlice;   ///< Low edge and upper edge of the current slice
        
        bool isNormalized;
        
#pragma endregion

#pragma region - random generator
        int64_t seed;
        std::mt19937_64 rndgen;
#pragma endregion

#pragma region - methods
        void setup();                ///< Optimize the probability distribution function for the random generator.
        void init();                ///< Initialize random generator state.
        double normalize_pdf();
#pragma endregion
#pragma endregion
        
#pragma region -- protected member
    protected:
        virtual void Uinit(){};
        virtual void Usetup(){};
#pragma endregion

    public:

#pragma region -- public member
#pragma region  - random generator
        double gen();                                        ///< Use of the MCMC slice sampling methods to generate a random number according to the PDF.
        double operator()(void);
        
        inline double normalize(){return normalize_pdf();}
        
        /**
         *  @brief Change the generator seed
         *  @details Modify the seed used by the random generator.
         *
         *  @param s: new seed
         */
        void setSeed(int64_t s)
        {
            seed = s;
            rndgen.seed(seed);
        }
        
        /**
         *  @brief Retrive the random generator seed.
         *  @details Retrive the random generator seed.
         *
         *  @return this\f$\rightarrow\f$seed
         */
        int64_t getSeed(){return seed;}

#pragma endregion
#pragma region - constructor/destructor
        binned_pdf();
        binned_pdf(farray, farray);
        binned_pdf(value , weight);
        binned_pdf(const binned_pdf &);
        
        virtual ~binned_pdf();

#pragma endregion
#pragma region - operator

        void operator+=(const binned_pdf&   );
        void operator*=(const double   );
        void operator/=(const double   );
        
        binned_pdf  operator+(const binned_pdf& );
        binned_pdf  operator*(const double );
        binned_pdf  operator/(const double );
        
        binned_pdf& operator=(const binned_pdf& );

#pragma endregion
#pragma region - PDF accessor
        
        inline size_t size() const {return fval.size();}
        
#pragma endregion
#pragma region - PDF modifier
        
        void scale(double);
        void shift(double);
        void fill(double, double w = 1.0);
        
#pragma endregion
#pragma region - PDF Evaluation
        double eval(float ) const;    //!<  Estimatition of the probability.
        double eval(int   ) const;    //!<  Estimatition of the probability.
        double eval(double) const;    //!<  Estimatition of the probability.
        
        double getMean() const;
        double getRMS()  const;
        
        double getCDF(double) const;
        double getCDF() const;
        
        int getBinWithMostProbableValue( double &) const;
        
        inline double valueAt(int i) const {return (fval.find(i)!=fval.end())?fval.find(i)->second:0;}
        inline double probabilityAt(int i) const {return (fpdf.find(i)!=fpdf.end())?fpdf.find(i)->second:0;}
        
#pragma endregion
#pragma region - Debug tools
        double LowEdge(){return lowEdge;}
        double UpEdge (){return upEdge ;}
        
        void dump() const;

#pragma endregion
#pragma region - Predifined PDF
        static binned_pdf& gaus(double, double, unsigned int nbin = 5000);
        static binned_pdf& logNormal(double, double, unsigned int nbin = 5000);
#pragma endregion
#pragma endregion
    };
        
#pragma endregion

#pragma region - binned_pdf2 class definition
        /**
         *  @class binned_pdf2 DSTstatistic.h "DSTmath/DSTstatistic.h"
         *  @brief Store 2D statistical probability distribution function
         *  @details The class binned_pdf2 store a 2D statistical probability distribution function (PDF) defined by a set of values and their associated probability. The PDF is used to generate random numbers according to the PDF using the MCMC slice sampling method.
         */
        class binned_pdf2
        {
        protected:
#pragma region -- protected member
#pragma region - variable
            lvalue  fxval;
            lvalue  fyval;
            lweight fpdf;
            
            long fnx;
            
            double fmax;                ///< The maximum of the pdf
            double x_lowEdge, x_upEdge;    ///< Low edge and upper edge of the random variable along the x axis
            double y_lowEdge, y_upEdge; ///< Low edge and upper edge of the random variable along the y axis
#pragma endregion
#pragma region - methods
            double normalize_pdf();
            virtual void Uinit(){};
            virtual void Usetup(){};

#pragma endregion
#pragma endregion
            
        private:
#pragma region -- private member     
            double xSlice, ySlice, zSlice;
            bool isNormalized;
            
            std::vector<double> rand_gen;
            
#pragma endregion
#pragma region - random generator
            int64_t seed;
            std::mt19937_64 rndgen;

#pragma endregion

#pragma region - methods
            void setup();                ///< Optimize the probability distribution function for the random generator.
            void init();                ///< Initialize random generator state.
#pragma endregion

        public:

#pragma region -- public member
#pragma region  - random generator
            inline double normalize(){return normalize_pdf();}
            
            std::vector<double> gen();               ///< Use of the MCMC slice sampling methods to generate a random number according to the PDF.
            std::vector<double> operator()(void);
            
            /**
             *  @brief Change the generator seed
             *  @details Modify the seed used by the random generator.
             *
             *  @param s: new seed
             */
            void setSeed(int64_t s)
            {
                seed = s;
                rndgen.seed(seed);
            }
            
            /**
             *  @brief Retrive the random generator seed.
             *  @details Retrive the random generator seed.
             *
             *  @return this\f$\rightarrow\f$seed
             */
            int64_t getSeed(){return seed;}
            
#pragma endregion
#pragma region - constructor/destructor
            binned_pdf2();
            binned_pdf2(farray,farray, farray);
            binned_pdf2(lvalue ,lvalue , lweight);
            binned_pdf2(const binned_pdf2 &);
            
            virtual ~binned_pdf2();
            
#pragma endregion
#pragma region - operator
            
            void operator+=(const binned_pdf2&   );
            void operator*=(const double   );
            void operator/=(const double   );
            
            binned_pdf2  operator+(const binned_pdf2& );
            binned_pdf2  operator*(const double );
            binned_pdf2  operator/(const double );
            
            binned_pdf2& operator=(const binned_pdf2& );
            
#pragma endregion
#pragma region - PDF modifier
            
            void scale(double, double);
            void shift(double, double);
            
            inline size_t nBins_x() const {return fxval.size()-1;}
            inline size_t nBins_y() const {return fyval.size()-1;}
            
            double x_width(int i) const;
            double y_width(int i) const;
            
            double x_center(int i) const;
            double y_center(int i) const;
            
            double bin_content(int i, int j) const;
            
#pragma endregion
#pragma region - PDF Evaluation
            double eval(float , float ) const;    //!<  Estimatition of the probability.
            double eval(int   , int   ) const;    //!<  Estimatition of the probability.
            double eval(double, double) const;    //!<  Estimatition of the probability.
            
            double getXMean() const;
            double getXRMS()  const;
            
            double getYMean() const;
            double getYRMS()  const;
            
            double getCDF(double, double) const;
            
#pragma endregion
#pragma region - Debug tools
            inline double LowEdge(unsigned short axis = 0) const {return (axis == 0) ? x_lowEdge : y_lowEdge;}
            inline double UpEdge (unsigned short axis = 0) const {return (axis == 0) ? x_upEdge  : y_upEdge ;}
            
            void dump() const;
            
#pragma endregion
#pragma region - Predifined PDF
            static binned_pdf2& gaus(double *, double *, unsigned int nbin = 5000);
            static binned_pdf2 gaus(double, double, double, double, unsigned int nbin = 5000);
            static binned_pdf2& logNormal(double *, double *, unsigned int nbin = 5000);
            static binned_pdf2  logNormal(double, double, double, double, unsigned int nbin = 5000);
#pragma endregion
#pragma endregion
        };
#pragma endregion

#pragma region - important_sampling class definition
    /**
     *  @class important_sampling DSTstatistic.h "DSTmath/DSTstatistic.h"
     *  @brief Random sampling of PDF using the importantsampling methods to probe tails of PDF.
     *  @details The class \c important_sampling offers a random sampling of probability distribution function (PDF). It uses the important-sampling statistical method to generate random number with an emphasis on the rare-envent's probabilities. This offers the possibitlity to probe the tails of PDF and studdy the effects and importances of rare events.
     */
    class important_sampling
    {
    private:
        binned_pdf * normal_pdf;
        binned_pdf * biased_pdf;
        
        const float fscale;
        
        void biasPDF();
        
        void reset();
        
    public:
        important_sampling(float scale= 2.5);                    //!<  Default constructor
        important_sampling(const binned_pdf &, float scale= 2.5);        //!<  Constructor
        important_sampling(const important_sampling&);            //!<  Copy constructor
        
        virtual ~important_sampling()                            //!<  Destructor
        {
            if(normal_pdf != NULL)
                delete normal_pdf;
            
            if(biased_pdf != NULL)
                delete biased_pdf;
        }
        
        double operator()(double &);
        
        void setSeed(int64_t);
        int64_t getSeed();
        
        void setPDF(const binned_pdf &);
        
        
    };
#pragma endregion

#pragma region - pdf class definition
    class pdf;
    
    typedef std::vector<double> pdf_param;
    typedef double (pdf::*pdf_function)(double, pdf_param) const;
    typedef std::pair<pdf_param,pdf_function> pdf_distribution;
    typedef std::vector<pdf_distribution> pdf_list;
    
    /**
     *  @class pdf DSTstatistic.h "DSTmath/DSTstatistic.h"
     *  @brief Analytical probability distribution function
     *  @details The class \c pdf is used to perform operation on unbinned probability distribution function of statiscal quantities.
     */
    class pdf
    {

#pragma region -- protected memeber
    protected:
        double fmax;
        double pdf_norme;
#pragma endregion 

#pragma region -- private memeber
    private:

        //pdf_function_list fpdf;
        pdf_list fpdf;
        
        double pdf_lowEdge;
        double pdf_upEdge;

        double xSlice, ySlice;

#pragma region - random generator
        int64_t seed;
        std::mt19937_64 rndgen;
#pragma endregion
#pragma endregion

#pragma region -- public memeber
    public:
#pragma region - constructor/destructor
        pdf ();
        ~pdf(){fpdf.clear();}
        
#pragma endregion

#pragma region - Probability function managment
        //void add_pdf(pdf_function);
        void add_pdf(pdf_function, pdf_param);
        void add_pdf(pdf_distribution);
        
#pragma endregion

#pragma region - PDF estimator
        double operator()(double);
        inline double get_maximum() const {return fmax;}
#pragma endregion

#pragma region - Configuration
        void setRange(double, double);
#pragma endregion

#pragma region - Random kernel
        double operator()(void);
        double gen(void);
#pragma endregion

#pragma region - operator
        inline void operator*=(double scale){pdf_norme *= scale;}
        inline void operator/=(double scale){pdf_norme /= scale;}
#pragma endregion
#pragma endregion
    };
        
#pragma endregion


#pragma region - normal_distribution class definition
    class normal_distribution: public pdf
    {
    private:
        double mean, rms;
        
    public:
        normal_distribution(double m, double s):pdf(),mean(m),rms(s)
        {
            pdf_param p={mean, rms};
            
            fmax = gauss(mean,p);
            pdf_norme = 1./sqrt(rms*rms*2.*acos(-1.));
            setRange(mean-10.*rms, mean+10.*rms);
            
            add_pdf(static_cast<pdf_function>(&normal_distribution::gauss), p);
        }
        
        normal_distribution(const normal_distribution &N)
        {
            mean = N.mean;
            rms  = N.rms;
            fmax = N.fmax;
            pdf_norme = N.pdf_norme;
            
            setRange(mean-10.*rms, mean+10.*rms);
            
            pdf_param p={mean, rms};
            
            add_pdf(static_cast<pdf_function>(&normal_distribution::gauss), p);
        }
        
        inline void set_sigma(double s){rms=s;}
        inline void set_mean (double m){mean=m;}
        
    protected:
        double gauss(double x, pdf_param p) const
        {
            return exp(-0.5*(x-p[0])*(x-p[0])/(p[1]*p[1]));
        }
    };
    

#pragma endregion

#pragma region - log_normal_distribution class definition
    class log_normal_distribution: public pdf
    {
    private:
        double mean, rms;
        
    public:
        log_normal_distribution(double m, double s):pdf(),mean(m),rms(s)
        {
            pdf_param p={mean, rms};
            
            fmax = log_normal(mean, p);
            pdf_norme = 1./sqrt(rms*rms*2.*acos(-1));
            setRange(mean-10.*rms, mean+10.*rms);
            
            add_pdf(static_cast<pdf_function>(&log_normal_distribution::log_normal), p);
        }
        
        log_normal_distribution(const log_normal_distribution &N)
        {
            mean = N.mean;
            rms  = N.rms;
            fmax = N.fmax;
            pdf_norme = N.pdf_norme;
            
            pdf_param p={mean, rms};
            
            setRange(0, mean*(1+3*rms));
            
            add_pdf(static_cast<pdf_function>(&log_normal_distribution::log_normal), p);
        }
        
    protected:
        double log_normal(double x, pdf_param p) const
        {
            if(x < 0)
                return 0;
            
            return 1./x*exp(-0.5*(log(x/p[0])*log(x/p[0]))/(p[1]*p[1]));
        }
        
    };
#pragma endregion
 
    }
}

#endif /* defined(_statistic_) */
