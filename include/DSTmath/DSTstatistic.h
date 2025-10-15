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
      
#if __cplusplus < 201103L
      int MOD_RAND = RAND_MAX/1000;
      double DMOD_RAND = static_cast<double>( MOD_RAND); 
#endif
        
	double get_normale_distribution ( double m_, double s_ );
	double get_poisson_distribution ( double l_ );
    
    template<typename T>
    T  get_uniform_distribution ( T, T );
    
    template<typename T>
    T  get_uniform_distribution ( T );

    
    template<typename T>
    T get_uniform_distribution ( T i_, T f_ )
    {
        T rndNum;
        static int seed = 0;
        
#if __cplusplus >= 201103L
        
        static std::mt19937_64 rndgen;
        
        if(seed == 0)
        {
            seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
            rndgen = std::mt19937_64(seed);
        }
        
        std::uniform_int_distribution<T> dis = std::uniform_int_distribution<T>(i_, f_);
        rndNum = dis(rndgen);
#else
        if(seed == 0 )
        {
            seed = time(NULL);
            srand(seed);
        }
        rndNum = static_cast<T>(static_cast<double>( rand() % MOD_RAND )/DMOD_RAND ) * ( f_ - i_ ) + i_;
#endif
            
        return rndNum;
        }
        
        template<typename T>
        T get_uniform_distribution ( T f_ )
        {
            return get_uniform_distribution<T>( 0, f_ );
        }
        
#if __cplusplus >= 201103L
        template<>
        double  get_uniform_distribution<double> ( double, double );
        
        template<>
        float  get_uniform_distribution<float> ( float, float );
        
        template<>
        long double  get_uniform_distribution<long double> ( long double, long double );
#endif
        
        


        
#pragma mark - pixmath::binned_pdf class definition
	/**
	 *  @class binned_pdf pixtemplate.h "pixsim/pixtemplate.h"
	 *  @brief Store statistical probability distribution function
	 *  @details To do
	 */
	class binned_pdf
	{
	private:
#pragma mark -- private member
#pragma mark - variable
		value  fval;
		weight fpdf;
		
		double xSlice, ySlice;
		double fmax;				///< The maximum of the pdf
		double lowEdge, upEdge;		///< Low edge and upper edge of the random variable
		double lowSlice, upSlice;	///< Low edge and upper edge of the current slice
		
		bool isNormalized;
		
#pragma mark - random generator
		int64_t seed;
        
#if __cplusplus >= 201103L
        std::mt19937_64 rndgen;
#endif


#pragma mark - methods
		void setup();				///< Optimize the probability distribution function for the random generator.
		void init();				///< Initialize random generator state.
		double normalize_pdf();
		
	protected:
		virtual void Uinit(){};
		virtual void Usetup(){};

	public:
#pragma mark -- public member
#pragma mark  - random generator
		double gen();										///< Use of the MCMC slice sampling methods to generate a random number according to the PDF.
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
#if __cplusplus >= 201103L
			rndgen.seed(seed);
#endif
		}
		
		/**
		 *  @brief Retrive the random generator seed.
		 *  @details Retrive the random generator seed.
		 *
		 *  @return this\f$\rightarrow\f$seed
		 */
		int64_t getSeed(){return seed;}

#pragma mark - constructor/destructor
		binned_pdf();
		binned_pdf(farray, farray);
		binned_pdf(value , weight);
		binned_pdf(const binned_pdf &);
		
		virtual ~binned_pdf();

#pragma mark - operator

		void operator+=(const binned_pdf&   );
		void operator*=(const double   );
		void operator/=(const double   );
		
		binned_pdf  operator+(const binned_pdf& );
		binned_pdf  operator*(const double );
		binned_pdf  operator/(const double );
		
		binned_pdf& operator=(const binned_pdf& );

#pragma mark - PDF accessor
        
        inline size_t size() const {return fval.size();}
        
#pragma mark - PDF modifier
		
		void scale(double);
		void shift(double);
        void fill(double, double w = 1.0);
		
#pragma mark - PDF Evaluation
		double eval(float ) const;	//!<  Estimatition of the probability.
		double eval(int   ) const;	//!<  Estimatition of the probability.
		double eval(double) const;	//!<  Estimatition of the probability.
		
		double getMean() const;
		double getRMS()  const;
		
		double getCDF(double) const;
        double getCDF() const;
        
        int getBinWithMostProbableValue( double &) const;
        
        inline double valueAt(int i) const {return (fval.find(i)!=fval.end())?fval.find(i)->second:0;}
        inline double probabilityAt(int i) const {return (fpdf.find(i)!=fpdf.end())?fpdf.find(i)->second:0;}
		
#pragma mark - Debug tools
		double LowEdge(){return lowEdge;}
		double UpEdge (){return upEdge ;}
		
		void dump() const;

#pragma mark - Predifined PDF
		static binned_pdf& gaus(double, double, unsigned int nbin = 5000);
		static binned_pdf& logNormal(double, double, unsigned int nbin = 5000);
	};
        
#pragma mark - pixmath::binned_pdf2 class definition
        /**
         *  @class binned_pdf2 pixtemplate.h "pixsim/pixtemplate.h"
         *  @brief Store 2D statistical probability distribution function
         *  @details To do
         */
        class binned_pdf2
        {
        protected:
#pragma mark -- private member
#pragma mark - variable
            lvalue  fxval;
            lvalue  fyval;
            lweight fpdf;
            
            long fnx;
            
            double fmax;				///< The maximum of the pdf
            double x_lowEdge, x_upEdge;	///< Low edge and upper edge of the random variable along the x axis
            double y_lowEdge, y_upEdge; ///< Low edge and upper edge of the random variable along the y axis
            
            double normalize_pdf();

            
        private:
            
            double xSlice, ySlice, zSlice;
            bool isNormalized;
            
            std::vector<double> rand_gen;
            
#pragma mark - random generator
            int64_t seed;
            
#if __cplusplus >= 201103L
            std::mt19937_64 rndgen;
#endif
            
            
#pragma mark - methods
            void setup();				///< Optimize the probability distribution function for the random generator.
            void init();				///< Initialize random generator state.
            
        protected:
            virtual void Uinit(){};
            virtual void Usetup(){};
            
        public:
#pragma mark -- public member
#pragma mark  - random generator
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
#if __cplusplus >= 201103L
                rndgen.seed(seed);
#endif
            }
            
            /**
             *  @brief Retrive the random generator seed.
             *  @details Retrive the random generator seed.
             *
             *  @return this\f$\rightarrow\f$seed
             */
            int64_t getSeed(){return seed;}
            
#pragma mark - constructor/destructor
            binned_pdf2();
            binned_pdf2(farray,farray, farray);
            binned_pdf2(lvalue ,lvalue , lweight);
            binned_pdf2(const binned_pdf2 &);
            
            virtual ~binned_pdf2();
            
#pragma mark - operator
            
            void operator+=(const binned_pdf2&   );
            void operator*=(const double   );
            void operator/=(const double   );
            
            binned_pdf2  operator+(const binned_pdf2& );
            binned_pdf2  operator*(const double );
            binned_pdf2  operator/(const double );
            
            binned_pdf2& operator=(const binned_pdf2& );
            
#pragma mark - PDF modifier
            
            void scale(double, double);
            void shift(double, double);
            
            inline size_t nBins_x() const {return fxval.size()-1;}
            inline size_t nBins_y() const {return fyval.size()-1;}
            
            double x_width(int i) const;
            double y_width(int i) const;
            
            double x_center(int i) const;
            double y_center(int i) const;
            
            double bin_content(int i, int j) const;
            
#pragma mark - PDF Evaluation
            double eval(float , float ) const;	//!<  Estimatition of the probability.
            double eval(int   , int   ) const;	//!<  Estimatition of the probability.
            double eval(double, double) const;	//!<  Estimatition of the probability.
            
            double getXMean() const;
            double getXRMS()  const;
            
            double getYMean() const;
            double getYRMS()  const;
            
            double getCDF(double, double) const;
            
#pragma mark - Debug tools
            inline double LowEdge(unsigned short axis = 0) const {return (axis == 0) ? x_lowEdge : y_lowEdge;}
            inline double UpEdge (unsigned short axis = 0) const {return (axis == 0) ? x_upEdge  : y_upEdge ;}
            
            void dump() const;
            
#pragma mark - Predifined PDF
            static binned_pdf2& gaus(double *, double *, unsigned int nbin = 5000);
            static binned_pdf2 gaus(double, double, double, double, unsigned int nbin = 5000);
            static binned_pdf2& logNormal(double *, double *, unsigned int nbin = 5000);
            static binned_pdf2  logNormal(double, double, double, double, unsigned int nbin = 5000);
        };
	
#pragma mark - pixmath::important_sampling class definition
	/**
	 *  @class pixmath::important_sampling pixtemplate.h "pixsim/pixtemplate.h"
	 *  @brief Random sampling of PDF using the importantsampling methods to probe tails of PDF.
	 *  @details The class pixmath::important_sampling offers a random sampling of probability distribution function (PDF). It uses the important-sampling statistical method to generate random number with an emphasis on the rare-envent's probabilities. This offers the possibitlity to probe the tails of PDF and studdy the effects and importances of rare events.
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
		important_sampling(float scale= 2.5);					//!<  Default constructor
		important_sampling(const binned_pdf &, float scale= 2.5);		//!<  Constructor
		important_sampling(const important_sampling&);			//!<  Copy constructor
		
		virtual ~important_sampling()							//!<  Destructor
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
	
#pragma mark - pixmath::pdf class definition
	class pdf;
	
    typedef std::vector<double> pdf_param;
    typedef double (pdf::*pdf_function)(double, pdf_param) const;
    typedef std::pair<pdf_param,pdf_function> pdf_distribution;
    typedef std::vector<pdf_distribution> pdf_list;
	
	/**
	 *  @class binned_pdf pixtemplate.h "pixsim/pixtemplate.h"
	 *  @brief Analytical probability function
	 *  @details The class pixmath::pdf is used to perform operation on unbinned probability distribution function of statiscal quantities.
	 */
	class pdf
	{
#pragma mark -- protected memeber
	protected:
		double fmax;
		double pdf_norme;
		
	private:
#pragma mark -- private memeber

		//pdf_function_list fpdf;
		pdf_list fpdf;
		
		double pdf_lowEdge;
		double pdf_upEdge;

		double xSlice, ySlice;
		
#pragma mark - random generator
		int64_t seed;

#if __cplusplus >= 201103L
		std::mt19937_64 rndgen;
#endif
		
	public:

#pragma mark -- public memeber
#pragma mark - constructor/destructor
		pdf ();
		~pdf(){fpdf.clear();}
		
#pragma mark - Probability function managment
		//void add_pdf(pdf_function);
        void add_pdf(pdf_function, pdf_param);
        void add_pdf(pdf_distribution);
		
#pragma mark - PDF estimator
		double operator()(double);
		inline double get_maximum() const {return fmax;}
		
#pragma mark - Configuration
		void setRange(double, double);
				
#pragma mark - Random kernel
		double operator()(void);
        double gen(void);
		
#pragma mark - operator

		inline void operator*=(double scale){pdf_norme *= scale;}
		inline void operator/=(double scale){pdf_norme /= scale;}
		
	};
		
#pragma mark - pixmath::normal_distribution class definition
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
	

#pragma mark - pixmath::log_normal_distribution class definition
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
 
    }
}

#endif /* defined(_statistic_) */
