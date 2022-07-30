//
//  template_utils.cc
//  
//
//  Created by GILLARD William on 03/02/15.
//  __MyCompanyName__
//	Copyright (c) 2015, All rights reserved
//

#include <iostream>
#include <limits>
#if __cplusplus >= 201103L
#include <chrono>
#else
#include <ctime>
#endif
#include <DSTmath/template_utils.h>

namespace DST
{
    namespace Math
    {
	double get_normale_distribution ( double m_, double s_ )
	{
		double x1, x2, y;
		x1 = drand48();
		x2 = drand48();
		
		// Boc-Muller method
		y = sqrt ( -2 * log ( x1 ) ) * cos ( 2. * M_PI * x2 );
		return m_ + s_ * y;
	}
	
	double get_poisson_distribution ( double l_ )
	{
		int k = 10. * drand48();
		int f = 1;
		for ( int i = 1; i <= k; i++ )
			f *= i;
		
		return exp ( - l_ ) * pow ( l_, k ) / f;
	}
        
#if __cplusplus >= 201103L
        
        template<>
        double get_uniform_distribution<double>( double i_, double f_ )
        {
            double rndNum;
            static int seed = 0;
            
            static std::mt19937_64 rndgen;
            
            if(seed == 0)
            {
                seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                rndgen = std::mt19937_64(seed);
            }
            
            std::uniform_real_distribution<double> dis = std::uniform_real_distribution<double>(i_, f_);
            rndNum = dis(rndgen);
            
            return rndNum;
        }
        
        template<>
        float get_uniform_distribution<float>( float i_, float f_ )
        {
            float rndNum;
            static int seed = 0;
            
            static std::mt19937_64 rndgen;
            
            if(seed == 0)
            {
                seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                rndgen = std::mt19937_64(seed);
            }
            
            std::uniform_real_distribution<float> dis = std::uniform_real_distribution<float>(i_, f_);
            rndNum = dis(rndgen);
            
            return rndNum;
        }
        
        template<>
        long double get_uniform_distribution<long double>( long double i_, long double f_ )
        {
            long double rndNum;
            static int seed = 0;
            
            static std::mt19937_64 rndgen;
            
            if(seed == 0)
            {
                seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                rndgen = std::mt19937_64(seed);
            }
            
            std::uniform_real_distribution<long double> dis = std::uniform_real_distribution<long double>(i_, f_);
            rndNum = dis(rndgen);
            
            return rndNum;
        }
#endif

	
#pragma mark -
#pragma mark - statistic::binned_pdf implementation
	
#pragma mark • constructor/destructor
 
	/**
	 *  @brief Default constructor
	 *  @details The default constructor initialize the probability distribution function to a null value. Then, the random generator will return 0 as long as the function distribution isn't defined.
	 */
	binned_pdf::binned_pdf()
	{
		init();
	}
	
	/**
	 *  @brief Constructor
	 *  @details Construc probability density function from a combination of value and their probability. The arrays of value and probability must have the same dimension. The array of value will be sorted in incrasing order and the array of probability will be sorted accordingly. Thereafter, the constructed probability density function is nomalized to 1.
	 *  @param V : array of value.
	 *  @param P : array of probabilities.
	 */
	binned_pdf::binned_pdf(farray V, farray P)
	{
		init();
		
		if(V.size() != P.size())
		{
			std::cerr<<"\033[43;31m *** ERROR *** \033[0m"<<std::endl<<"\033[31m    [statistic::binned_pdf::binned_pdf] Variable and their associated probability \033[43;31m MUST \033[31m have same legth.\033[31m"<<std::endl;
			return;
		}
		
		// BUILD A TEMPORARY ARRAY TO SORT INPUT VALUE IN INCREASING ORDER
		bool sorted = true;
		std::map<double, unsigned int> sortArray = std::map<double, unsigned int>();
		for(unsigned int k = 0; k < V.size(); k++)
		{
			if(sortArray.size() > 0)
			{
				if(V[k] <= (--sortArray.end())->first && sorted == true)
				{
					std::cerr<<"\033[43;31m !!! WARNING !!! \033[0m"<<std::endl<<"\033[31m    [statistic::binned_pdf::binned_pdf] Variable array isn't in increasing order. Array will be sorted and its associated probability array will be sorted accordingly.\033[31m"<<std::endl;
					sorted = false;
				}
			}
			
			sortArray.insert(std::pair<double, unsigned int>(V[k],k));
		}
		
		unsigned int n = 0;
		for(std::map<double,unsigned int>::iterator it = sortArray.begin(); it != sortArray.end(); it++)
		{
			fval[n] = it->first;
			fpdf[n] = P[it->second];
		}
		
		setup();
	}
	
	/**
	 *  @brief Constructor
	 *  @details Construct probability density function from a combination of value and their probability. Thereafter, the constructed probability density function is nomalized to 1.
	 *  @param V : array of value.
	 *  @param P : array of probabilities.
	 */
	binned_pdf::binned_pdf(value V, weight P)
	{
		init();
		
		if(V.size() != P.size())
		{
			std::cerr<<"\033[43;31m *** ERROR *** \033[0m"<<std::endl<<"\033[31m    [statistic::binned_pdf::binned_pdf] Variable and their associated probability \033[43;31m MUST \033[31m have same legth.\033[31m"<<std::endl;
			return;
		}
		
		weight::iterator ft = P.end();
		for(value::iterator it = V.begin(); it != V.end(); it++)
		{
			ft = P.find(it->first);
			if( ft == P.end() )
				continue;
			
			fval[it->first] = it->second;
			fpdf[ft->first] = ft->second;
		}
		
		setup();
	}
	
	/**
	 *  @brief Copy constructor.
	 *  @details A copy of the probability distribution fuction is stored in memory. The two PDF (this one and the original) aren't sharing sharing the same random generator. Thus, they could be both used to generate intependent random number based on the same distribution.
	 *  @param PDF : Original probability distribution fucntion to be copied into this.
	 */
	binned_pdf::binned_pdf(const binned_pdf& PDF)
	{
		init();
		
		fval.insert(PDF.fval.begin(),PDF.fval.end());
		fpdf.insert(PDF.fpdf.begin(),PDF.fpdf.end());
		
		setup();
	}
	
	/**
	 *  @brief Destructor.
	 *  @details Clear probability distribution function from the memory
	 */
	binned_pdf::~binned_pdf()
	{
		fval.clear();
		fpdf.clear();
		
		init();
	}
	
#pragma mark • Initialization
	/**
	 * @details Initialize the random generator state to default value. The random seed is generated based on the clock time.
	 */
	void binned_pdf::init()
	{
	
		fval   = value();
		fpdf   = weight();
		
		xSlice = std::numeric_limits<double>::min();
		ySlice = std::numeric_limits<double>::min();
		fmax   = std::numeric_limits<double>::min();
		upEdge  = xSlice;
		lowEdge = xSlice;
		
		isNormalized = false;
		
		Uinit();

#if __cplusplus >= 199711L
		// obtain a seed from the timer
		seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		
        rndgen = std::mt19937_64(seed);
#else
		seed = time(NULL);
		srand(seed);
#endif

	}
	
	/**
	 *  @details Private method used to optimize the PDF that will be passed to the random generator. PDF ranges are extracted and the PDF is normalized to 1.\n The virtual methods statistic::binned_pdf::Usetup() is called upon return.
	 */
	void binned_pdf::setup()
	{
		lowEdge = fval.begin()->second;
		upEdge  = (--fval.end())->second;
		
		//--- NORMALIZE THE PDF DISTRIBUTION
		fmax = normalize_pdf();
		
		//--- CHOOSE INITIAL VALUE FOR THE RANDOM GENERATOR.
 
#if __cplusplus >= 199711L
		std::uniform_real_distribution<double>dis(upEdge, lowEdge);
		xSlice = dis(rndgen);
#else
        xSlice = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND *( upEdge - lowEdge ) + lowEdge;
#endif
		
		Usetup();
	}
	
	/**
	 *  @brief PDF nomalization
	 *  @details Compute the integral of the PDF and normalize the PDF such that its integral is equal to 1. The integration is done using Sipson's rule of integration.
	 *  @return The maximum probability that the PDF can take.
	 */
	double binned_pdf::normalize_pdf()
	{
		if(isNormalized)
			return fmax;
		
		double norm = getCDF((--fval.end())->second);
		double max = 0;
		
		if(norm < 0)
			std::cerr<<"\033[43;31m *** ERROR *** \033[0m"<<std::endl<<"\033[31m     [statistic::binned_pdf::binned_pdf] Negative nomalization : \033[34m"<<norm<<"\033[0m"<<std::endl;
		
		for(weight::iterator ft = fpdf.begin(); ft != fpdf.end(); ft++)
		{
			ft->second /= norm;
			
			if(ft->second >= max)
				max = ft->second;
		}
		
		isNormalized = true;
		
		return max;
	}
	
#pragma mark • random generator
	/**
	 *  @brief Generate random number based on the PDF.
	 *  @details Markov-Chain Monte-Carlo slice sampling algorithm is used to sample the pdf and randomly generate a number according to the probability. Futhermore, in order to allow rare event to occurs, an important sampling is combined to the slice sampling algortim to artificially increase probability of the tails of the distribution.
	 *  @return A random number.
	 */
	double binned_pdf::gen()
	{
		if(fval.size() < 1)
			return lowEdge-1;
		
		if(fabs(xSlice) <= std::numeric_limits<double>::min() || xSlice <= fval.begin()->second || xSlice >= (--fval.end())->second )
		  setup();

		if(!isNormalized)
			fmax = normalize_pdf();
		
		if(fmax <= 0)
			return lowEdge-1;
        
        //- RANDOMLY PICK A SLICE IN Y
#if __cplusplus >= 199711L
        std::uniform_real_distribution<double> ydis(0, eval(xSlice));
        ySlice = ydis(rndgen);
        while(ySlice >= fmax) ySlice = ydis(rndgen);
#else
        ySlice = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND * eval( xSlice );
        while(ySlice >= fmax)
	  ySlice = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*eval( xSlice );
#endif
    
		//- FIND X BOUNDARY FOR THE Y SLICE
	double min, max;
	weight::const_iterator ik = fpdf.begin();
        while(ik->second < ySlice && ik != fpdf.end()) ik++;
	min = fval.find(ik->first)->second;
		
	ik = fpdf.end();
	ik--;
	while(ik->second < ySlice && ik != fpdf.begin()) ik--;
	max = fval.find(ik->first)->second;
		
	lowSlice = (min <= max)?min:max;
	upSlice  = (max >= min)?max:min;

	if(upSlice == lowSlice)
	  {
	    double step = (upSlice - lowSlice)/100.;
	    upSlice     += step;
	    lowSlice    -= step;
	  }
		
		//- RANDOMLY PICK A X WITHIN THE SLICE BOUNDARY
        double xgen;
        
#if __cplusplus >= 199711L
        std::uniform_real_distribution<double> xdis(lowSlice, upSlice);
        xgen = xdis(rndgen);
        while(eval(xgen) < ySlice)	xgen = xdis(rndgen);
#else
        xgen = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*( upSlice - lowSlice ) + lowSlice;
        while(eval(xgen) < ySlice)
	  xgen = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*( upSlice - lowSlice ) + lowSlice;
#endif
		
		xSlice = xgen;
		
		return xSlice;
	}
	
	/**
	 *  @brief Generate random number based on the PDF.
	 *  @details Markov-Chain Monte-Carlo slice sampling algorithm is used to sample the pdf and randomly generate a number according to the probability. Futhermore, in order to allow rare event to occurs, an important sampling is combined to the slice sampling algortim to artificially increase probability of the tails of the distribution.
	 *  @return A random number.
	 */
	double binned_pdf::operator()(void)
	{
		return gen();
	}

#pragma mark • PDF Evaluation
	/**
	 *  @details Estimate the probability of \f$x\f$ to be selected. If \f$x\f$ is outside the PDF boundary, their isn't any interpolation done and the returned probability will be 0.
	 *  @param x: a value to test.
	 *  @return the probability associated to the value of \f$x\f$
	 */
	double binned_pdf::eval(int x) const
	{
		return eval(static_cast<double>( x ));
	}
	
	/**
	 *  @details Estimate the probability of \f$x\f$ to be selected. If \f$x\f$ is outside the PDF boundary, their isn't any interpolation done and the returned probability will be 0.
	 *  @param x: a value to test.
	 *  @return the probability associated to the value of \f$x\f$
	 */
	double binned_pdf::eval(float x) const
	{
		return eval(static_cast<double>( x ));
	}
	
	/**
	 *  @details Estimate the probability of \f$x\f$ to be selected. If \f$x\f$ is outside the PDF boundary, their isn't any interpolation done and the returned probability will be 0.
	 *  @param x: a value to test.
	 *  @return the probability associated to the value of \f$x\f$
	 */
	double binned_pdf::eval(double x) const
	{
		if(x < lowEdge || x > upEdge)
			return 0;
		
		if(fval.size() < 1)
			return 0;
		
		value::const_iterator it;
		
		double a=0, b=0, fa=0, fb=0, slope=0, cte=0;
		
		if((x - lowEdge) <= (upEdge - x))
		{
			it = fval.begin();
			
			while (it->second < x && it != fval.end())
			{
				a = it->second;
				fa= fpdf.find(it->first)->second;
				it++;
			}
			
			if(it == fval.end())
				return fa;
			
			b = it->second;
			fb= fpdf.find(it->first)->second;
		}
		else
		{
			it = fval.end();
			it --;
			while (it->second > x && it != fval.begin())
			{
				b = it->second;
				fb= fpdf.find(it->first)->second;
				
				it--;
			}
			
			a = it->second;
			fa= fpdf.find(it->first)->second;

		}
		
		slope = (fb-fa)/(b-a);
		cte   = fb - slope*b;
		
		return slope*x + cte;
	}
	
	/**
	 *  @details Print out onto stderr the PDF distribution
	 */
	void binned_pdf::dump() const
	{
		weight::const_iterator ft = fpdf.end();
		value ::const_iterator ik = fval.begin();
		
        std::cout<<"\033[31m--------------------------------------------------------\033[0m"<<std::endl
                 <<" Mean : "<<getMean()<<std::endl
                 <<" RMS : "<<getRMS()<<std::endl<<std::endl;
		std::cout<<"#	x	P(x)	CDF(x)"<<std::endl;
		while(ik != fval.end())
		{
			ft = fpdf.find(ik->first);
			if(ft == fpdf.end())
				continue;
			
			std::cout<<ik->first<<"	"<<ik->second<<"	"<<ft->second<<"	"<<getCDF(ik->second)<<std::endl;
			
			ik++;
		}
		
	}
	
   /**
    *  @details Use of the Simpson's rule to compute the integral from the PDF lowEdge to  \f$x\f$ and return the cumulated distribution at \f$x\f$
    *
    *  @param x: upper limit to compute the cumulatied distribution function
    *
    *  @return cumulated distribution at \f$x\f$.
    */
   double binned_pdf::getCDF() const
   {
       return getCDF((--fval.end())->second);
   }
        
	/**
	 *  @details Use of the Simpson's rule to compute the integral from the PDF lowEdge to  \f$x\f$ and return the cumulated distribution at \f$x\f$
	 *
	 *  @param x: upper limit to compute the cumulatied distribution function
	 *
	 *  @return cumulated distribution at \f$x\f$.
	 */
	double binned_pdf::getCDF(double x) const
	{
		double norm = 0;
		
		value::const_iterator it = fval.begin();
		while( it->second < x && it != fval.end() )
		{
			double a = it->second;
			double fa= fpdf.find(it->first)->second;
			
			it++; if(it == fval.end()) continue;
			double f0 = fpdf.find(it->first)->second;
			
			it++; if(it == fval.end()) continue;
			double b = it->second;
			double fb= fpdf.find(it->first)->second;
			
			if(b < a)
				std::cerr<<"\033[43;31m *** ERROR *** \033[0m"<<std::endl<<"\033[0m     [statistic::binned_pdf::binned_pdf] Array value not srted in increasing order.\033[0m"<<std::endl;
			
			if((fa + 4.*f0 + fb) < 0)
				std::cerr<<"\033[31m[statistic::binned_pdf::binned_pdf]\033[43;31m !!! WARNING !!!\033[0m"<<std::endl<<"\033[31m     Will give a negative norm.\033[0m"<<std::endl
				<<"     |- b - a :"<<b<<" - "<<a<<std::endl
				<<"     |- fa    :"<<fa<<std::endl
				<<"     |- f0    :"<<f0<<std::endl
				<<"     `- fb    :"<<fb<<std::endl;
			
				
			norm += (b-a)/6.*(fa + 4.*f0 + fb);
		}
		
		return norm;
	}
	
	/**
	 *  @brief Get the mean of the PDF
	 *  @details Evaluate the mean value of the PDF
	 *
	 *  @return \f$\langle x\rangle = \frac{\sum\left( x_i\times g(x_i) \right)}{\sum g(x_i)}\f$
	 */
	double binned_pdf::getMean() const
	{
		double mean = 0;
		double count= 0;
		
		value::const_iterator  it = fval.begin();
		weight::const_iterator ik = fpdf.begin();
		
		while( it != fval.end() && ik != fpdf.end())
		{
			mean += it->second*ik->second;
			count += ik->second;
			
			it++;
			ik++;
		}
		
		if(count == 0)
			std::cerr<<"\033[43;31m *** ERROR *** \033[0m"<<std::endl<<"\033[31m     [statistic::binned_pdf::binned_pdf] Can't estimate mean of a NULL PDF\033[0m"<<std::endl;
		
		mean/=count;
		
		return mean;
	}
	
	/**
	 *  @brief Get the root mean square of the PDF
	 *  @details Evaluate the root mean square of the PDF
	 *
	 *  @return \f$\langle x^2\rangle = \sqrt{\frac{\sum x_i^2\times g(x_i)}{\sum g(x_i)}}\f$
	 */
	double binned_pdf::getRMS() const
	{
		double RMS = 0;
		double count= 0;
		
		value::const_iterator  it = fval.begin();
		weight::const_iterator ik = fpdf.begin();
		
		while( it != fval.end() && ik != fpdf.end())
		{
            if( ik->second <= std::numeric_limits<double>::min() )
                continue;
            
			RMS += ik->second*(it->second)*(it->second);
			count += ik->second;
			
			it++;
			ik++;
		}
		
        RMS /= count;
		RMS = sqrt(RMS);
		
		if(count == 0)
			std::cerr<<"\033[43;31m *** ERROR *** \033[0m"<<std::endl<<"\033[31m     [statistic::binned_pdf::binned_pdf] Can't estimate RMS of a NULL PDF\033[0m"<<std::endl;
		
		return RMS;
	}
        
    /**
     *  @details Find the bin that contains the modt probable value
     *  @details MPV upon return the most probable value
     *  @return Bin with the most probable value
     */
    int binned_pdf::getBinWithMostProbableValue(double& MPV) const
    {
        double max = 0;
        int iBin = -1;
        for(weight::const_iterator ik = fpdf.begin(); ik != fpdf.end(); ik++)
        {
            if(ik->second >= max)
            {
                max = ik->second;
                iBin=ik->first;
            }
        }
        
        MPV = fval.find(iBin)->second;
        return iBin;
    }

#pragma mark • Operator
#pragma mark – Self-modifier opperator

	/**
	 *  @brief Mathematical operator
	 *  @details Sum probability of two PDF
	 *  @param P PDF to add to this.
	 *  @return this\f$\rightarrow\f$fpdf = this\f$\rightarrow\f$fpdf + P
	 */
	void binned_pdf::operator+=(const binned_pdf &P)
	{
		value::const_iterator it = P.fval.begin();
		value::iterator it_this  = fval.end();
		
		//-- IDENTIFY INDEX OF FIRST AND LAST ELEMENT
		int it_first = fval.begin()->first;
		
		double dx = (1-(fval.begin()->second - (--fval.end())->second)/fval.size()/3.);
		
		if( P.lowEdge > upEdge  || P.upEdge < lowEdge )
			std::cerr<<"\033[31m[statistic::binned_pdf::binned_pdf]\033[43;31m !!! WARNING !!! \033[0m"<<std::endl<<"\033[33m     Their isn't any overlap between the two distribution. Their will be discontinuities.\033[0m"<<std::endl;
		else if( P.lowEdge != lowEdge  && P.upEdge != upEdge )
			std::cerr<<"\033[31m[statistic::binned_pdf::binned_pdf]\033[43;31m !!! WARNING !!! \033[0m"<<std::endl<<"\033[33m     The two distributions have different bundaries. Their might be discontinuities.\033[0m"<<std::endl;
		
		if(lowEdge*dx > P.lowEdge)
		{
			// COUNT NUMBER OF ELEMENT BELLOW LOW EDGE
			while(it->second < lowEdge*dx && it != P.fval.end()){it_first--; it++;}
		}
		
		it = P.fval.begin();
		double weight = 0;
		while(it != P.fval.end())
		{
			if(it == P.fval.end())
				break;
			
			//---- 3 posible case :
			if(it->second < lowEdge*dx)
			{
				// A) P.fval[it] bellow the lower limit of this PDF : this PDF range is extended with negativ indexes.
				
				weight = (P.fpdf.find(it->first))->second;
				
				fval.insert(std::pair<int,double>(it_first, it->second));
				fpdf.insert(std::pair<int,double>(it_first, weight    ));
				it_first++;
			}
			else if (it->second > lowEdge && it->second < upEdge)
			{
				// B) P.fval[it] is in the range of this pdf. Weight are incremented by the weight of P.
				
				if(fabs(it->second - fval.begin()->second) < fabs(it->second - (--fval.end())->second))
				{
					it_this = fval.begin();
					while(it_this->second < it->second && it_this != fval.end()) it_this++;
					
					weight = P.eval(it_this->second);
					if(weight > 0)
						fpdf.find(it_this->first)->second += weight;
				}
				else
				{
					it_this = fval.end();
					it_this--;
					while(it_this->second > it->second && it_this != fval.begin()) it_this--;
					
					weight = P.eval(it_this->second);
					if(weight > 0)
						fpdf.find(it_this->first)->second += weight;
				}
			}
			else if(it->second > upEdge)
			{
				// C) P.fval[it] is above the upper limite of this pdf. this pdf is extended.
				weight = (P.fpdf.find(it->first))->second;
				
				fval.insert(std::pair<int,double>((--fval.end())->first+1, it->second));
				fpdf.insert(std::pair<int,double>((--fpdf.end())->first+1, weight    ));
			}
			
			it++;
			
		}
		
		isNormalized = false;
		
		lowEdge = (lowEdge > P.lowEdge) ? P.lowEdge : lowEdge;
		upEdge  = (upEdge  < P.upEdge ) ? P.upEdge  : upEdge ;

	}
	
	/**
	 *  @brief Mathematical operator
	 *  @details Scale probability.
	 *  @param scaler Multiplication factor used to scale the probabilities.
	 *  @return this\f$\rightarrow\f$fpdf = this\f$\rightarrow\f$fpdf * scaler
	 */
	void binned_pdf::operator*=(const double scaler)
	{
		if(scaler <= 0)
			std::cerr<<"\033[31m[statistic::binned_pdf::binned_pdf]\033[43;31m !!! WARNING !!! \033[0m"<<std::endl<<"\033[33m     Scaling a PDF with negative number isn't physical.\033[0m"<<std::endl;
		
		for( weight::iterator it = fpdf.begin(); it != fpdf.end(); it++)
		{
			(it->second) *= scaler;
		}
		
		isNormalized = false;
		
	}
	
	/**
	 *  @brief Mathematical operator
	 *  @details Scale probability.
	 *  @param scaler Multiplication factor used to scale the probabilities.
	 *  @return this\f$\rightarrow\f$fpdf = this\f$\rightarrow\f$fpdf / scaler
	 */
	void binned_pdf::operator/=(const double scaler)
	{
		operator*=(1./scaler);
		
		isNormalized = false;
	}
	
#pragma mark – modifier operator
	/**
	 *  @brief Mathematical operator
	 *  @details Sum up probability between two pdf.
	 *  @param scaler Multiplication factor used to scale the probabilities.
	 *  @return a new PDF distribtuion B such that B.fpdf = this\f$\rightarrow\f$fpdf + P.fpdf
	 */
	binned_pdf binned_pdf::operator+(const binned_pdf& P)
	{
		binned_pdf P2 = *this;
		
		P2 += P;
		
		return P2;
	}
	
	/**
	 *  @brief Mathematical operator
	 *  @details Scale the probability of PDF.
	 *  @param scaler Multiplication factor used to scale the probabilities.
	 *  @return a new PDF distribtuion B such that B.fpdf = this\f$\rightarrow\f$fpdf * scaler
	 */
	binned_pdf binned_pdf::operator*(const double scaler)
	{
		binned_pdf P2 = *this;
		
		P2*=scaler;
		
		return P2;
	}
	
	/**
	 *  @brief Mathematical operator
	 *  @details Scale the probability of PDF.
	 *  @param scaler Multiplication factor used to scale the probabilities.
	 *  @return a new PDF distribtuion B such that B.fpdf = this\f$\rightarrow\f$fpdf / scaler
	 */
	binned_pdf binned_pdf::operator/(const double scaler)
	{
		binned_pdf P2 = *this;
		
		P2/=scaler;
		
		return P2;
	}
	
	/**
	 *  @details Shift the mean of the pdf distribution by a given offset.
	 *  @param offset The offset used to shift the distribution
	 *  @note After shifting the PDF, the distribution is re-normalized to 1.
	 */
	void binned_pdf::shift(double offset)
	{
		double mean = getMean();
		
		for(value::iterator it = fval.begin(); it != fval.end(); it++)
		{
			double oldVal = it->second;
			it->second = ((oldVal - mean) + offset) + mean;
		}
		
		lowEdge =     fval.begin()->second;
		upEdge  = (-- fval.end() )->second;
		
		isNormalized = false;
		normalize_pdf();
	}
	
	/**
	 *  @details Scale the RMS of the pdf distribution.
	 *  @param scaleFactor The scaling factor used to scale the RMS of the distribution
	 *  @note After scalling the distribution, the scaled PDF is re-normalized to 1.
	 */
	void binned_pdf::scale(double scaleFactor)
	{
		double mean = getMean();
		for(value::iterator it = fval.begin(); it != fval.end(); it++)
		{
			double oldVal = it->second;
			it->second = (oldVal - mean)*scaleFactor + mean;
		}
		
		lowEdge =     fval.begin()->second;
		upEdge  = (-- fval.end() )->second;
		
		isNormalized = false;
		normalize_pdf();
	}
        
   /**
    *  @details Increment probability of the given bin.
    *  @param P the bin to increment probability
    *  @param W the weight to increment
    */
    void binned_pdf::fill(double P, double W)
    {
        bool match = false;
        value::iterator it = fval.begin();
        
        while(!match || it != fval.end())
        {
            if(P <= it->second)
            {
                match=true;
                continue;
            }
        }
        
        if(it == fval.end())
            return;
        
        if(it == fval.begin() && P <= it->second)
            return;
        
        weight::iterator iw = fpdf.find(it->first);
        if(iw == fpdf.end())
            return;
        
        iw->second += W;
        
    }
	
#pragma mark – assignement operator
	/**
	 *  @brief Assignement operator
	 *  @details Copy this to a new statistic::binned_pdf PDF distribution. The random generator seed ins,t copyed during the copy to insure independance between the two random number generators
	 *  @param P A PDF to copy.
	 *  @return a new PDF distribtuion which is a copy of P.
	 */
	binned_pdf & binned_pdf::operator=(const binned_pdf &P)
	{
		if(this != &P)
		{
			binned_pdf tmp(P);
			
			std::swap(fval, tmp.fval);
			std::swap(fpdf, tmp.fpdf);
			std::swap(lowEdge, tmp.lowEdge);
			std::swap(upEdge, tmp.upEdge);
			
			setup();
		}
		
		return *this;
	}
	
#pragma mark • Predifined PDF
	/**
	 *  @details Generate a binned gaussian distribution
	 *
	 *  @param mu    Mean of the gaussian
	 *  @param sigma standard deviation of the gaussian
	 *  @param nbin  Number of bin use to sample the gaussian distribution
	 *
	 *  @return A binned probability distribution with a gaussian shape.
	 */
	binned_pdf & binned_pdf::gaus(double mu, double sigma, unsigned int nbin)
	{
		value val;
		weight w ;
		
		unsigned int k = 0;
		for(double x0 = mu-10.*fabs(sigma); x0 < mu+10.*fabs(sigma); x0 += 20*fabs(sigma)/static_cast<double>(nbin))
		{
			val[k] = x0;
			w  [k] = exp(-0.5*(x0 - mu)*(x0 - mu)/sigma/sigma);
			
			k++;
		}
		
		binned_pdf *P = new binned_pdf(val, w);
		
		return *P;
	}
	
	/**
	 *  @details Generate a binned log-normal distribution
	 *
	 *  @param mu    Mean of the log-normal
	 *  @param sigma standard deviation of the log-normal
	 *  @param nbin  Number of bin use to sample the log-normal distribution
	 *
	 *  @return A binned probability distribution with a gaussian shape.
	 */
	binned_pdf & binned_pdf::logNormal(double mu, double sigma, unsigned int nbin)
	{
		value val;
		weight w ;
		
		unsigned int k = 0;
		double low = (mu-10.*fabs(sigma) < 0)? 0 : mu-10.*fabs(sigma);
		double up  = mu+10.*fabs(sigma);
		
		for(double x0 = low; x0 < up; x0 += fabs(up-low)/static_cast<double>(nbin))
		{
			if(x0 == 0)
				continue;
			
			val[k] = x0;
			w  [k] = 1./x0*exp(-0.5*(log(x0) - mu)*(log(x0) - mu)/sigma/sigma);
			
			k++;
		}
		
		binned_pdf *P = new binned_pdf(val, w);
		
		return *P;
	}
        
#pragma mark - statistic::binned_pdf2 implementation
        
#pragma mark • constructor/destructor
        
        /**
         *  @brief Default constructor
         *  @details The default constructor initialize the probability distribution function to a null value. Then, the random generator will return 0 as long as the function distribution isn't defined.
         */
        binned_pdf2::binned_pdf2()
        {
            init();
        }
        
        /**
         *  @brief Constructor
         *  @details Construc probability density function from a combination of value and their probability. The arrays of value and probability must have the same dimension. The array of value will be sorted in incrasing order and the array of probability will be sorted accordingly. Thereafter, the constructed probability density function is nomalized to 1.
         *  @param V : array of value.
         *  @param P : array of probabilities.
         */
        binned_pdf2::binned_pdf2(farray Vx, farray Vy, farray P)
        {
            init();
            
            if(Vx.size()*Vy.size() != P.size())
            {
                std::cerr<<"\033[43;31m *** ERROR *** \033[0m"<<std::endl<<"\033[31m    [statistic::binned_pdf2::binned_pdf2] Variables and their associated probability \033[43;31m MUST \033[31m have same legth.\033[31m"<<std::endl;
                return;
            }
            
            // BUILD A TEMPORARY ARRAY TO SORT INPUT VALUE IN INCREASING ORDER
            bool sorted_x = true;
            std::map<double, unsigned int> sortArray_x = std::map<double, unsigned int>();
            for(unsigned int k = 0; k < Vx.size(); k++)
            {
                if(sortArray_x.size() > 0)
                {
                    if(Vx[k] <= (--sortArray_x.end())->first && sorted_x == true)
                    {
                        std::cerr<<"\033[43;31m !!! WARNING !!! \033[0m"<<std::endl<<"\033[31m    [statistic::binned_pdf2::binned_pdf2] X variable array isn't in increasing order. Array will be sorted and its associated probability array will be sorted accordingly.\033[31m"<<std::endl;
                        sorted_x = false;
                    }
                }
                
                sortArray_x.insert(std::pair<double, unsigned int>(Vx[k],k));
            }
            
            // BUILD A TEMPORARY ARRAY TO SORT INPUT VALUE IN INCREASING ORDER
            bool sorted_y = true;
            std::map<double, unsigned int> sortArray_y = std::map<double, unsigned int>();
            for(unsigned int k = 0; k < Vy.size(); k++)
            {
                if(sortArray_y.size() > 0)
                {
                    if(Vy[k] <= (--sortArray_y.end())->first && sorted_y == true)
                    {
                        std::cerr<<"\033[43;31m !!! WARNING !!! \033[0m"<<std::endl<<"\033[31m    [statistic::binned_pdf2::binned_pdf2] Y variable array isn't in increasing order. Array will be sorted and its associated probability array will be sorted accordingly.\033[31m"<<std::endl;
                        sorted_y = false;
                    }
                }
                
                sortArray_y.insert(std::pair<double, unsigned int>(Vy[k],k));
            }
            
            while (Vy.size()/static_cast<unsigned int>(fnx) > 1)
                fnx *= 10;
            
            unsigned int nx = 0;
            unsigned int ny = 0;
            for(std::map<double,unsigned int>::iterator ix = sortArray_x.begin(); ix != sortArray_x.end(); ix++)
            {
                if(nx >= Vx.size())
                    nx = 0;
                
                fxval[nx] = ix->first;
                
                ny = 0;
                for(std::map<double,unsigned int>::iterator iy = sortArray_y.begin(); iy != sortArray_y.end(); iy++)
                {
                    fyval[ny] = iy->first;
                    fpdf[nx*fnx + ny]  = P[ix->second*fnx + iy->second];
                    ny++;
                }
                nx++;
            }
            
            setup();
        }
        
        /**
         *  @brief Constructor
         *  @details Construct probability density function from a combination of value and their probability. Thereafter, the constructed probability density function is nomalized to 1.
         *  @param V : array of value.
         *  @param P : array of probabilities.
         */
        binned_pdf2::binned_pdf2(lvalue Vx, lvalue Vy, lweight P)
        {
            init();
            
            if(Vx.size()*Vy.size() != P.size())
            {
                std::cerr<<"\033[43;31m *** ERROR *** \033[0m"<<std::endl<<"\033[31m    [statistic::binned_pdf2::binned_pdf2] Variables and their associated probability \033[43;31m MUST \033[31m have same legth.\033[31m"<<std::endl;
                std::cout<<"   - x axis dimension : "<<Vx.size()<<std::endl
                         <<"   - y axis dimension : "<<Vy.size()<<std::endl
                         <<"   - array dimension  : "<<P.size ()<<std::endl;
                return;
            }
            
            while (Vy.size()/static_cast<unsigned int>(fnx) > 1)
                fnx *= 10;
            
            lweight::iterator ft = P.end();
            for(lvalue::const_iterator ix = Vx.begin(); ix != Vx.end(); ix++)
            for(lvalue::const_iterator iy = Vy.begin(); iy != Vy.end(); iy++)
            {
                ft = P.find(ix->first*fnx + iy->first);
                
                if( ft == P.end() )
                {
                    std::cout<<"NOT FOUND"<<std::endl;
                    continue;
                }
                
                fxval[ix->first] = ix->second;
                fyval[iy->first] = iy->second;
                fpdf [ft->first] = ft->second;
            }
            
            setup();

        }
        
        /**
         *  @brief Copy constructor.
         *  @details A copy of the probability distribution fuction is stored in memory. The two PDF (this one and the original) aren't sharing sharing the same random generator. Thus, they could be both used to generate intependent random number based on the same distribution.
         *  @param PDF : Original probability distribution fucntion to be copied into this.
         */
        binned_pdf2::binned_pdf2(const binned_pdf2& PDF)
        {
            init();
            
            fxval.insert(PDF.fxval.begin(),PDF.fxval.end());
            fyval.insert(PDF.fyval.begin(),PDF.fyval.end());
            fpdf .insert(PDF.fpdf.begin() ,PDF.fpdf.end());
            fnx  = PDF.fnx;
            x_upEdge  = PDF.x_upEdge;
            x_lowEdge = PDF.x_lowEdge;
            
            y_upEdge  = PDF.y_upEdge;
            y_lowEdge = PDF.y_lowEdge;
            
            setup();
        }
        
        /**
         *  @brief Destructor.
         *  @details Clear probability distribution function from the memory
         */
        binned_pdf2::~binned_pdf2()
        {
            fxval.clear();
            fyval.clear();
            fpdf.clear();
            
            init();
        }
        
#pragma mark • Initialization
        /**
         * @details Initialize the random generator state to default value. The random seed is generated based on the clock time.
         */
        void binned_pdf2::init()
        {
            
            fxval  = lvalue();
            fyval  = lvalue();
            fpdf   = lweight();
            
            xSlice = std::numeric_limits<double>::min();
            ySlice = std::numeric_limits<double>::min();
            zSlice = std::numeric_limits<double>::min();
            fmax   = std::numeric_limits<double>::min();
            
            x_upEdge  = xSlice;
            x_lowEdge = xSlice;
            
            y_upEdge  = ySlice;
            y_lowEdge = ySlice;
            
            isNormalized = false;
            fnx = 10;
            
            Uinit();
            
#if __cplusplus >= 199711L
            // obtain a seed from the timer
            seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
            
            rndgen = std::mt19937_64(seed);
#else
            seed = time(NULL);
            srand(seed);
#endif
            
        }
        
        /**
         *  @details Private method used to optimize the PDF that will be passed to the random generator. PDF ranges are extracted and the PDF is normalized to 1.\n The virtual methods statistic::binned_pdf::Usetup() is called upon return.
         */
        void binned_pdf2::setup()
        {
            rand_gen = std::vector<double>();
            
            x_lowEdge = fxval.begin()->second;
            x_upEdge  = (--fxval.end())->second;
            
            y_lowEdge = fyval.begin()->second;
            y_upEdge  = (--fyval.end())->second;
            
            if(fnx == 10)
            {
                while (fyval.size()/static_cast<unsigned int>(fnx) > 1)
                    fnx *= 10;
            }
            
            //--- NORMALIZE THE PDF DISTRIBUTION
            fmax = normalize_pdf();
            
            //--- SORT X and Y axis
            /*std::map<double, int> xsort = std::map<double, int>();
            std::map<double, int> ysort = std::map<double, int>();
            weight tmp = weight();
            
            for(value::const_iterator it = fxval.begin(); it != fxval.end(); it++)
                xsort[it->second] = it->first;

            for(value::const_iterator it = fyval.begin(); it != fyval.end(); it++)
                ysort[it->second] = it->first;
            
            for(value::const_iterator ix = fxval.begin(); ix != fxval.end(); ix++)
            for(value::const_iterator iy = fyval.begin(); iy != fyval.end(); iy++)
                tmp[ix->first * fnx + iy->first] = fpdf[ix->first * fnx + iy->first];
            
            int ik = 0;
            int jk = 0;
            for(std::map<double, int>::const_iterator jx = xsort.begin(); jx != xsort.end(); jx++)
            {
                fxval[ik] = jx->first;
                jk = 0;
                for(std::map<double, int>::const_iterator jy = ysort.begin(); jy != ysort.end(); jy++)
                {
                    fyval[jk] = jy->first;
                    fpdf[ik * fnx + jk] = fpdf[jx->second * fnx + jy->second];
                    
                    jk++;
                }
                
                ik++;
            }
            
            xsort.clear();
            ysort.clear();
            tmp.clear();
             */
            
            
            //--- CHOOSE INITIAL VALUE FOR THE RANDOM GENERATOR.
            
#if __cplusplus >= 199711L
            std::uniform_real_distribution<double>xdis(x_upEdge, x_lowEdge);
            xSlice = xdis(rndgen);
            
            std::uniform_real_distribution<double>ydis(y_upEdge, y_lowEdge);
            ySlice = ydis(rndgen);
#else
            xSlice = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND *( x_upEdge - x_lowEdge ) + x_lowEdge;
            ySlice = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND *( y_upEdge - y_lowEdge ) + y_lowEdge;
#endif
            
            Usetup();
        }
        
        /**
         *  @brief PDF nomalization
         *  @details Compute the integral of the PDF and normalize the PDF such that its integral is equal to 1. The integration is done using Sipson's rule of integration.
         *  @return The maximum probability that the PDF can take.
         */
        double binned_pdf2::normalize_pdf()
        {
            if(isNormalized)
                return fmax;
            
            double norm = getCDF((--fxval.end())->second, (--fyval.end())->second);
            double max = 0;
            
            if(norm < 0)
                std::cerr<<"\033[43;31m *** ERROR *** \033[0m"<<std::endl<<"\033[31m     [statistic::binned_pdf2::binned_pdf2] Negative nomalization : \033[34m"<<norm<<"\033[0m"<<std::endl;
            
            for(lweight::iterator ft = fpdf.begin(); ft != fpdf.end(); ft++)
            {
                ft->second /= norm;
                
                if(ft->second >= max)
                    max = ft->second;
            }
            
            isNormalized = true;
            
            return max;
        }
        
#pragma mark • random generator
        /**
         *  @brief Generate random number based on the PDF.
         *  @details Markov-Chain Monte-Carlo slice sampling algorithm is used to sample the pdf and randomly generate a number according to the probability. Futhermore, in order to allow rare event to occurs, an important sampling is combined to the slice sampling algortim to artificially increase probability of the tails of the distribution.
         *  @return A random number.
         */
        std::vector<double> binned_pdf2::gen()
        {
            
            
            if(fxval.size() < 1)
            {
                rand_gen.clear();
                rand_gen.push_back(x_lowEdge-1);
                rand_gen.push_back(y_lowEdge-1);
                return rand_gen;
            }
            
            if(fyval.size() < 1)
            {
                rand_gen.clear();
                rand_gen.push_back(x_lowEdge-1);
                rand_gen.push_back(y_lowEdge-1);
                return rand_gen;
            }
            
            //if(fabs(xSlice) <= std::numeric_limits<double>::min() || xSlice <= fxval.begin()->second || xSlice >= (--fxval.end())->second )
            //    setup();
            
            //if(fabs(ySlice) <= std::numeric_limits<double>::min() || ySlice <= fyval.begin()->second || ySlice >= (--fyval.end())->second )
            //    setup();
            
            if(!isNormalized)
                fmax = normalize_pdf();
            
            if(fmax <= 0)
            {
                rand_gen.clear();
                rand_gen.push_back(x_lowEdge-1);
                rand_gen.push_back(y_lowEdge-1);
                return rand_gen;
            }
            
            rand_gen.clear();
            
            //- RANDOMLY PICK A SLICE IN Y
#if __cplusplus >= 199711L
            std::uniform_real_distribution<double> zdis(0, eval(xSlice, ySlice));
            zSlice = zdis(rndgen);
            while(zSlice >= fmax) zSlice = zdis(rndgen);
#else
            zSlice = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND * eval( xSlice, ySlice );
            while(zSlice >= fmax)
                zSlice = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*eval( xSlice, ySlice );
#endif
            
            
            //- RANDOMLY PICK A X WITHIN THE SLICE BOUNDARY
            double xgen, ygen;
            
#if __cplusplus >= 199711L
            std::uniform_real_distribution<double> xdis(x_lowEdge, x_upEdge);
            std::uniform_real_distribution<double> ydis(y_lowEdge, y_upEdge);
            xgen = xdis(rndgen);
            ygen = ydis(rndgen);
            while(eval(xgen, ygen) < zSlice)	{xgen = xdis(rndgen); ygen = ydis(rndgen);}
#else
            xgen = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*( x_upEdge - x_lowEdge ) + x_lowEdge;
            ygen = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*( y_upEdge - y_lowEdge ) + y_lowEdge;
            while(eval(xgen, ygen) < zSlice)
            {
                xgen = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*( x_upEdge - x_lowEdge ) + x_lowEdge;
                ygen = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*( y_upEdge - y_lowEdge ) + y_lowEdge;
            }
#endif
            
            xSlice = xgen;
            ySlice = ygen;
            
            rand_gen.push_back(xSlice);
            rand_gen.push_back(ySlice);
            
            return rand_gen;
        }
        
        /**
         *  @brief Generate random number based on the PDF.
         *  @details Markov-Chain Monte-Carlo slice sampling algorithm is used to sample the pdf and randomly generate a number according to the probability. Futhermore, in order to allow rare event to occurs, an important sampling is combined to the slice sampling algortim to artificially increase probability of the tails of the distribution.
         *  @return A random number.
         */
        std::vector<double> binned_pdf2::operator()(void)
        {
            return gen();
        }
        
#pragma mark • PDF Evaluation
        
        double binned_pdf2::x_width(int i) const
        {
            lvalue::const_iterator it = fxval.find(i);
            if(it == fxval.end() || it == --fxval.end())
                return 0;
            
            double x0 = it->second;
            
            it++;
            
            return fabs(it->second-x0);
        }
        
        double binned_pdf2::y_width(int i) const
        {
            lvalue::const_iterator it = fyval.find(i);
            if(it == fyval.end() || it == --fyval.end())
                return 0;
            
            double x0 = it->second;
            
            it++;
            
            return fabs(it->second-x0);
        }
        
        double binned_pdf2::x_center(int i) const
        {
            lvalue::const_iterator it = fxval.find(i);
            if(it == fxval.end())
                return std::numeric_limits<double>::max();
            
            return it->second+x_width(i)/2;
        }
        
        double binned_pdf2::y_center(int i) const
        {
            lvalue::const_iterator it = fyval.find(i);
            if(it == fyval.end())
                return std::numeric_limits<double>::max();
            
            return it->second+y_width(i)/2;
        }
        
        double binned_pdf2::bin_content(int i, int j) const
        {
            lweight::const_iterator itw = fpdf.find(i*fnx+j);
            if(itw == fpdf.end())
                return 0;
            
            return itw->second;
        }
        
        /**
         *  @details Estimate the probability of \f$x\f$ to be selected. If \f$x\f$ is outside the PDF boundary, their isn't any interpolation done and the returned probability will be 0.
         *  @param x: a value to test.
         *  @return the probability associated to the value of \f$x\f$
         */
        double binned_pdf2::eval(int x, int y) const
        {
            return eval(static_cast<double>( x ), static_cast<double>( y ));
        }
        
        /**
         *  @details Estimate the probability of \f$x\f$ to be selected. If \f$x\f$ is outside the PDF boundary, their isn't any interpolation done and the returned probability will be 0.
         *  @param x: a value to test.
         *  @return the probability associated to the value of \f$x\f$
         */
        double binned_pdf2::eval(float x, float y) const
        {
            return eval(static_cast<double>( x ), static_cast<double>( y ));
        }
        
        /**
         *  @details Estimate the probability of \f$x\f$ to be selected. If \f$x\f$ or \f$y\f$ is outside the PDF boundary, their isn't any interpolation done and the returned probability will be 0. For any combination of \f$x\f$ and \f$y\f$ coordinates in allowed range, the value of the closest neighboorhood is returnd.
         *  @param x: a value to test.
         *  @param y: a value to test.
         *  @return the probability associated to the coordinates of \f$x\f$ and \f$y\f$
         */
        double binned_pdf2::eval(double x, double y) const
        {
            if(x <= x_lowEdge || x >= x_upEdge || y <= y_lowEdge || y >= y_upEdge)
                return 0;
            
            if(fxval.size() < 1 || fyval.size() < 1)
                return 0;
            
            lvalue::const_iterator itx1 = fxval.begin();
            lvalue::const_iterator ity1 = fyval.begin();
            
            //-- Find x1 the closest to x
            while( itx1->second < x && itx1 != fxval.end() )
                itx1++;
            
            if(itx1 == fxval.end())
                return 0;
            
            //-- Find y1 the closest to y
            while( ity1->second < y && ity1 != fyval.end() )
                ity1++;
            
            if(ity1 == fyval.end())
                return 0;
         
            long ix = (itx1 != fxval.begin()) ? itx1->first - 1 : itx1->first;
            long iy = (ity1 != fyval.begin()) ? ity1->first - 1 : ity1->first;
            double D = sqrt( pow(x - fxval.find(ix)->second,2.) + pow(y - fyval.find(iy)->second, 2.));
            
            for(long i = ix; i < ix+1; i++)
            for(long j = iy; j < iy+1; j++)
            {
                if( sqrt( pow(x - fxval.find(i)->second,2.) + pow(y - fyval.find(j)->second, 2.)) < D )
                {
                    D = sqrt( pow(x - fxval.find(i)->second,2.) + pow(y - fyval.find(j)->second, 2.));
                    itx1 = fxval.find(i);
                    ity1 = fyval.find(j);
                }
            }
            
            
            return fpdf.find( itx1->first * fnx + ity1->first )->second;
        }
        
        
        /**
         *  @details Print out onto stderr the PDF distribution
         */
        void binned_pdf2::dump() const
        {
            lweight::const_iterator ft = fpdf.begin();
            lvalue ::const_iterator ix = fxval.begin();
            lvalue ::const_iterator iy = fyval.begin();
            
            std::cout<<"\033[31m--------------------------------------------------------\033[0m"<<std::endl
                     <<" Mean X : "<<getXMean()<<"   -    Mean Y : "<<getYMean()<<std::endl
                     <<" RMS  X : "<<getXRMS()<<"    -    RMS  Y : "<<getYRMS()<<std::endl<<std::endl;
            std::cout<<"#   x   y   P(x)   CDF(x)"<<std::endl;
            while(ft != fpdf.end())
            {
                int xbin = ft->first / fnx;
                int ybin = ft->first - xbin*fnx;
                
                ix = fxval.find(xbin);
                iy = fyval.find(ybin);
                if(ix == fxval.end() || iy == fyval.end())
                {
                    ft++;
                    continue;
                }
                
                std::cout<<ft->first<<"    "<<iy->second<<"    "<<ix->second<<"    "<<ft->second<<"    "<<getCDF(ix->second,iy->second)<<std::endl;
                
                ft++;
            }
            
        }
        
        /**
         *  @details Use of the Simpson's rule to compute the integral from the PDF lowEdge to  \f$(x,y)\f$ and return the cumulated distribution at \f$(x,y)\f$. The interval \f$[x_0,x]\f$ is subdivided in \f$m\f$ subinterval of equal size \f$h=\frac{x-x_0}{m}\f$ and the interval \f$[y_0,y]\f$ is subdivided in \f$n\f$ subinterval of equal size \f$k=\frac{y-y_0}{n}\f$ where \f$n\f$ and \f$m\f$ correspond to the number of bin in range \f$[x_0,x]\f$ and \f$[y_0,y]\f$, respectively. Then, the 2D trapezoidal interpolation together with the simpson rules are used to compute the integral :
         \f[
         \int_{x_0}^x\int_{y_0}^y f(x,y)dxdy = \frac{1}{4}hk\left(f(x_0,y_0)+f(x_0,y)+f(x,y_0)+f(x,y)+2\sum_{i=0}^{m-1}f(x_i,y_0)+2\sum_{i=0}^{m-1}f(x_i,y)+2\sum_{j=0}^{n-1}f(x_0,y_j)+2\sum_{j=0}^{n-1}f(x,y_j)+4\sum_{j=0}^{n-1}\sum_{i=0}^{m-1}f(x_i,y_j)
         \f]
         *
         *  @param x: upper limit along the x axis to compute the cumulatied distribution function
         *  @param y: upper limit along the y to compute the cumulatied distribution function
         *
         *  @return cumulated distribution at \f$\int_{x_0}^x\int_{y_0}^y f(x,y)dxdy\f$.
         */
        double binned_pdf2::getCDF(double x, double y) const
        {
            if( (x < x_lowEdge) || (x > x_upEdge) || (y < y_lowEdge) || (y > y_upEdge) )
                return 0;
            
            //FIND fxval nearest to x
            unsigned int m = 1;
            lvalue::const_iterator ix = fxval.begin();
            while (ix->second < x && ix != fxval.end())
            {
                ix++;
                m++;
            }
            if(m == fxval.size())
                m--;
            
            //FIND fyval nearest to x
            unsigned int n = 1;
            lvalue::const_iterator iy = fyval.begin();
            while (iy->second < y  && iy != fyval.end() )
            {
                iy++;
                n++;
            }
            if(n == fyval.size())
                n--;

            
            double h = (x - x_lowEdge)/static_cast<double>(m);
            double k = (y - y_lowEdge)/static_cast<double>(n);
            
            double norm = 0;
            
            norm = eval(x_lowEdge, y_lowEdge) + eval(x_lowEdge, y) + eval(x, y_lowEdge) + eval(x,y);
            for(unsigned int i = 0; i < m; i++)
                norm += 2.*eval( x_lowEdge + h*static_cast<double>(i), y_lowEdge ) + 2.*eval( x_lowEdge + h*static_cast<double>(i), y );
            
            for(unsigned int j = 0; j < n; j++)
                norm += 2.*eval( x_lowEdge, y_lowEdge + k*static_cast<double>(j) ) + 2.*eval( x, y_lowEdge + k*static_cast<double>(j) );
            
            for(unsigned int j = 0; j < n; j++)
            for(unsigned int i = 0; i < m; i++)
                norm += 4.*eval( x_lowEdge + h*static_cast<double>(i), y_lowEdge + k*static_cast<double>(j) );
            
            norm *= 1./4.*h*k;
            
            
            return norm;
        }
        
        
        
        /**
         *  @brief Get the mean \f$x\f$ value
         *  @details Evaluate the mean \f$x\f$ value
         *
         *  @return \f$\langle x\rangle = \frac{\sum\left( x_i\times g(x_i) \right)}{\sum g(x_i)}\f$ where \f$g(x_i)=\sum_{j=0}^{n-1}f(x_i,y_j)\f$
         */
        double binned_pdf2::getXMean() const
        {
            double mean = 0;
            
            weight w = weight();
            value  x = value() ;
            
            for(lvalue::const_iterator ix = fxval.begin(); ix != fxval.end(); ix++)
            {
                w[static_cast<int>(ix->first)]=0;
                x[static_cast<int>(ix->first)]=ix->second;
            
                for(lvalue::const_iterator iy = fyval.begin(); iy != fyval.end(); iy++)
                    w.find(static_cast<int>(ix->first))->second += fpdf.find(ix->first*fnx + iy->first)->second;
            }
            
            binned_pdf tmp = binned_pdf(x,w);
            
            mean = tmp.getMean();
            
            x.clear();
            w.clear();
            tmp.~binned_pdf();
            
            return mean;
        }
        
        /**
         *  @brief Get the mean \f$y\f$ value
         *  @details Evaluate the mean \f$y\f$ value
         *
         *  @return \f$\langle y\rangle = \frac{\sum\left( y_i\times g(y_i) \right)}{\sum g(y_i)}\f$ where \f$g(y_i)=\sum_{j=0}^{m-1}f(x_j,y_i)\f$
         */
        double binned_pdf2::getYMean() const
        {
            double mean = 0;
            
            weight w;
            value  y = value() ;
            
            for(lvalue::const_iterator iy = fyval.begin(); iy != fyval.end(); iy++)
            {
                w[static_cast<int>(iy->first)]=0;
                y[static_cast<int>(iy->first)]=iy->second;
                
                for(lvalue::const_iterator ix = fxval.begin(); ix != fxval.end(); ix++)
                    w[static_cast<int>(iy->first)] += fpdf.find(ix->first*fnx + iy->first)->second;
            }
            
            binned_pdf tmp = binned_pdf(y,w);
            
            mean = tmp.getMean();
            
            y.clear();
            w.clear();
            tmp.~binned_pdf();
            
            return mean;
        }
        
        /**
         *  @brief Get the standard deviation of the PDF along the x axis
         *  @details Evaluate the standard defiation of the PDF along the x axis
         *
         *  @return \f$\langle x^2\rangle = \sqrt{\frac{\sum x_i^2\times g(x_i)^2 - \langle x\rangle^2}{\sum g(x_i)^2}}\f$ where \f$g(x_i)=\sum_{j=0}^{n-1}f(x_i,y_j)\f$
         */
        double binned_pdf2::getXRMS() const
        {
            double RMS = 0;
            
            weight w;
            value  x = value() ;
            
            for(lvalue::const_iterator ix = fxval.begin(); ix != fxval.end(); ix++)
            {
                w[static_cast<int>(ix->first)]=0;
                x[static_cast<int>(ix->first)]=ix->second;
                
                for(lvalue::const_iterator iy = fyval.begin(); iy != fyval.end(); iy++)
                    w[ix->first] += fpdf.find(ix->first*fnx + iy->first)->second;
            }
            
            binned_pdf tmp = binned_pdf(x,w);
            
            RMS = tmp.getRMS();
            
            x.clear();
            w.clear();
            tmp.~binned_pdf();
            
            return RMS;
        }
        
        /**
         *  @brief Get the standard deviation of the PDF along the y axis
         *  @details Evaluate the standard defiation of the PDF along the y axis
         *
         *  @return \f$\langle y^2\rangle = \sqrt{\frac{\sum y_i^2\times g(y_i)^2 - \langle y\rangle^2}{\sum g(y_i)^2}}\f$ where \f$g(y_i)=\sum_{j=0}^{n-1}f(x_j,y_i)\f$
         */
        double binned_pdf2::getYRMS() const
        {
            double RMS = 0;
            
            weight w;
            value  y = value() ;

            for(lvalue::const_iterator iy = fyval.begin(); iy != fyval.end(); iy++)
            {
                w[static_cast<int>(iy->first)]=0;
                y[static_cast<int>(iy->first)]=iy->second;
                
                for(lvalue::const_iterator ix = fxval.begin(); ix != fxval.end(); ix++)
                    w[iy->first] += fpdf.find(ix->first*fnx + iy->first)->second;
            }
            
            binned_pdf tmp = binned_pdf(y,w);
            
            RMS = tmp.getRMS();
            
            y.clear();
            w.clear();
            tmp.~binned_pdf();
            
            return RMS;
        }
        
#pragma mark • Operator
#pragma mark – Self-modifier opperator
        
        /**
         *  @brief Mathematical operator
         *  @details Sum probability of two PDF
         *  @param P PDF to add to this.
         *  @return this\f$\rightarrow\f$fpdf = this\f$\rightarrow\f$fpdf + P
         */
        void binned_pdf2::operator+=(const binned_pdf2 &P)
        {
            bool re_map = false;
            double dx = 0;
            double dy = 0;
            
            if( P.x_lowEdge > x_upEdge  || P.x_upEdge < x_lowEdge || P.y_lowEdge > y_upEdge  || P.y_upEdge < y_lowEdge )
            {
                std::cerr<<"\033[31m[statistic::binned_pdf2]\033[43;31m !!! WARNING !!! \033[0m"<<std::endl<<"\033[33m     Their isn't any overlap between the two distribution. I remap the current pdf to account for new boundarys.\033[0m"<<std::endl;
                re_map = true;
            }
            else if( (P.x_lowEdge < x_lowEdge  || P.x_upEdge > x_upEdge) || (P.y_lowEdge < y_lowEdge || P.y_upEdge > y_upEdge) )
            {
                std::cerr<<"\033[31m[statistic::binned_pdf2]\033[43;31m !!! WARNING !!! \033[0m"<<std::endl<<"\033[33m     The two distributions have different bundaries. I remap the current pdf to account for new boundarys.\033[0m"<<std::endl;
                re_map = true;
            }
            
            if(re_map)
            {
                //GET AVERAGE dx and dy of this
                dx = 0;
                dy = 0;
                
                for(lvalue::const_iterator ix = fxval.begin(); ix != (--fxval.end()); ix++)
                    dx+=fxval.find(ix->first+1)->second - ix->second;
                dx/=static_cast<double>(fxval.size()-1);
                
                for(lvalue::const_iterator iy = fyval.begin(); iy != (--fyval.end()); iy++)
                    dy+=fyval.find(iy->first+1)->second - iy->second;
                dy/=static_cast<double>(fyval.size()-1);
                
                double x_min = std::min(x_lowEdge, P.x_lowEdge);
                double x_max = std::max(x_upEdge , P.x_upEdge);
                double y_min = std::min(y_lowEdge, P.y_lowEdge);
                double y_max = std::max(y_upEdge , P.y_upEdge);
                
                lvalue tmpx;
                lvalue tmpy;
                lvalue tmp_pdf;
                
                long nelem = static_cast<long>((y_max-y_min)/dy);
                long f     = 10;
                while(nelem/f > 1)
                    f*=10;
                
                long ix=0, iy=0;
                
                for(long i = 0; i < static_cast<long>((x_max-x_min)/dx); i++)
                {
                    double x = dx*static_cast<double>( i ) + x_min;
                    if( (x >   x_upEdge && x < P.x_lowEdge) ||
                        (x > P.x_upEdge && x <   x_lowEdge) )
                        continue;
                    
                    lvalue::const_iterator ii = (--fxval.end());
                    while(x > ii->second && ii != fxval.begin())
                        ii--;
                    
                    if(x_width(ii->first) > 10*dx && (x > P.x_upEdge || x < P.x_lowEdge))
                        continue;
                    
                    tmpx[ix]=x;
                    
                    iy=0;
                    for(long j = 0; j < static_cast<long>((y_max-y_min)/dy); j++)
                    {
                        double y = dy*static_cast<double>( j ) + y_min;
                        if( (y >   y_upEdge && x < P.y_lowEdge) ||
                            (y > P.y_upEdge && x <   y_lowEdge) )
                            continue;
                        
                        lvalue::const_iterator ij = (--fyval.end());
                        while(y > ij->second && ij != fyval.begin())
                            ij--;
                        
                        if(y_width(ij->first) > 10*dy && (y > P.y_upEdge || y < P.y_lowEdge))
                            continue;
                        
                        tmpy[iy]=y;
                        tmp_pdf[ix*f+iy]=eval(tmpx[ix], tmpy[iy]) + P.eval(tmpx[ix], tmpy[iy]);
                        
                        iy++;
                    }
                    
                    ix++;
                }
                
                std::swap(fxval, tmpx);
                std::swap(fyval, tmpy);
                std::swap(fnx, f);
                std::swap(fpdf, tmp_pdf);
                
                x_lowEdge = x_min;
                x_upEdge  = x_max;
                
                y_lowEdge = y_min;
                y_upEdge  = y_max;
                
                tmpx.clear();
                tmpy.clear();
                tmp_pdf.clear();
            }
            else
            {
                for(lvalue::const_iterator ix = fxval.begin(); ix != fxval.end(); ix++)
                for(lvalue::const_iterator iy = fyval.begin(); iy != fyval.end(); iy++)
                {
                    fpdf.find(ix->first*fnx + iy->first)->second += P.eval(ix->second, iy->second);
                }
            }
            
            isNormalized = false;
            setup();
        }
        
        /**
         *  @brief Mathematical operator
         *  @details Scale probability.
         *  @param scaler Multiplication factor used to scale the probabilities.
         *  @return this\f$\rightarrow\f$fpdf = this\f$\rightarrow\f$fpdf * scaler
         */
        void binned_pdf2::operator*=(const double scaler)
        {
            if(scaler <= 0)
                std::cerr<<"\033[31m[statistic::binned_pdf::binned_pdf]\033[43;31m !!! WARNING !!! \033[0m"<<std::endl<<"\033[33m     Scaling a PDF with negative number isn't physical.\033[0m"<<std::endl;
            
            for( lweight::iterator it = fpdf.begin(); it != fpdf.end(); it++)
            {
                (it->second) *= scaler;
            }
            
            isNormalized = false;
            
        }
        
        /**
         *  @brief Mathematical operator
         *  @details Scale probability.
         *  @param scaler Multiplication factor used to scale the probabilities.
         *  @return this\f$\rightarrow\f$fpdf = this\f$\rightarrow\f$fpdf / scaler
         */
        void binned_pdf2::operator/=(const double scaler)
        {
            operator*=(1./scaler);
            
            isNormalized = false;
        }
        
#pragma mark – modifier operator
        /**
         *  @brief Mathematical operator
         *  @details Sum up probability between two pdf.
         *  @param scaler Multiplication factor used to scale the probabilities.
         *  @return a new PDF distribtuion B such that B.fpdf = this\f$\rightarrow\f$fpdf + P.fpdf
         */
        binned_pdf2 binned_pdf2::operator+(const binned_pdf2& P)
        {
            binned_pdf2 P2 = *this;
            
            P2 += P;
            
            return P2;
        }
        
        /**
         *  @brief Mathematical operator
         *  @details Scale the probability of PDF.
         *  @param scaler Multiplication factor used to scale the probabilities.
         *  @return a new PDF distribtuion B such that B.fpdf = this\f$\rightarrow\f$fpdf * scaler
         */
        binned_pdf2 binned_pdf2::operator*(const double scaler)
        {
            binned_pdf2 P2 = *this;
            
            P2*=scaler;
            
            return P2;
        }
        
        /**
         *  @brief Mathematical operator
         *  @details Scale the probability of PDF.
         *  @param scaler Multiplication factor used to scale the probabilities.
         *  @return a new PDF distribtuion B such that B.fpdf = this\f$\rightarrow\f$fpdf / scaler
         */
        binned_pdf2 binned_pdf2::operator/(const double scaler)
        {
            binned_pdf2 P2 = *this;
            
            P2/=scaler;
            
            return P2;
        }
        
        /**
         *  @details Shift the mean of the pdf distribution by a given offset.
         *  @param offset The offset used to shift the distribution
         *  @note After shifting the PDF, the distribution is re-normalized to 1.
         */
        void binned_pdf2::shift(double xoffset, double yoffset)
        {
            double meanX = getXMean();
            double meanY = getYMean();
            
            for(lvalue::iterator itx = fxval.begin(); itx != fxval.end(); itx++)
            {
                double oldVal = itx->second;
                itx->second = ((oldVal - meanX) + xoffset) + meanX;
            }
            
            x_lowEdge =     fxval.begin()->second;
            x_upEdge  = (-- fxval.end() )->second;
            
            for(lvalue::iterator ity = fyval.begin(); ity != fyval.end(); ity++)
            {
                double oldVal = ity->second;
                ity->second = ((oldVal - meanY) + yoffset) + meanY;
            }
            
            y_lowEdge =     fyval.begin()->second;
            y_upEdge  = (-- fyval.end() )->second;
            
            isNormalized = false;
            normalize_pdf();
        }
        
        /**
         *  @details Scale the RMS of the pdf distribution.
         *  @param scaleFactor The scaling factor used to scale the RMS of the distribution
         *  @note After scalling the distribution, the scaled PDF is re-normalized to 1.
         */
        void binned_pdf2::scale(double XscaleFactor, double YscaleFactor)
        {
            double meanX = getXMean();
            double meanY = getYMean();
            for(lvalue::iterator itx = fxval.begin(); itx != fxval.end(); itx++)
            {
                double oldVal = itx->second;
                itx->second = (oldVal - meanX)*XscaleFactor + meanX;
            }
            
            for(lvalue::iterator ity = fyval.begin(); ity != fyval.end(); ity++)
            {
                double oldVal = ity->second;
                ity->second = (oldVal - meanY)*YscaleFactor + meanY;
            }
            
            x_lowEdge =     fxval.begin()->second;
            x_upEdge  = (-- fxval.end() )->second;
            
            y_lowEdge =     fyval.begin()->second;
            y_upEdge  = (-- fyval.end() )->second;
            
            isNormalized = false;
            normalize_pdf();
        }
        
#pragma mark – assignement operator
        /**
         *  @brief Assignement operator
         *  @details Copy this to a new statistic::binned_pdf PDF distribution. The random generator seed ins,t copyed during the copy to insure independance between the two random number generators
         *  @param P A PDF to copy.
         *  @return a new PDF distribtuion which is a copy of P.
         */
        binned_pdf2 & binned_pdf2::operator=(const binned_pdf2 &P)
        {
            if(this != &P)
            {
                binned_pdf2 tmp(P);
                
                std::swap(fxval, tmp.fxval);
                std::swap(fyval, tmp.fyval);
                std::swap(fnx, tmp.fnx);
                std::swap(fpdf, tmp.fpdf);
                std::swap(x_lowEdge, tmp.x_lowEdge);
                std::swap(x_upEdge, tmp.x_upEdge);
                std::swap(y_lowEdge, tmp.y_lowEdge);
                std::swap(y_upEdge, tmp.y_upEdge);
                
                setup();
            }
            
            return *this;
        }
        
#pragma mark • Predifined PDF
        binned_pdf2 binned_pdf2::gaus(double mux, double muy, double sigmax, double sigmay, unsigned int nbin)
        {
            double *mu = new double[2]; mu[0]=mux; mu[1]=muy;
            double *sigma = new double[2]; sigma[0]=sigmax; sigma[1]=sigmay;
            
            binned_pdf2 G = gaus(mu, sigma, nbin);
            delete [] mu;
            delete [] sigma;
            
            return G;
        }
        
        /**
         *  @details Generate a binned gaussian distribution
         *
         *  @param mu    Mean of the gaussian
         *  @param sigma standard deviation of the gaussian
         *  @param nbin  Number of bin use to sample the gaussian distribution
         *
         *  @return A binned probability distribution with a gaussian shape.
         */
        binned_pdf2 & binned_pdf2::gaus(double *mu, double *sigma, unsigned int nbin)
        {
            std::cout<<"Generate 2D gaus distribution with :"<<std::endl
                     <<"sigma = ("<<sigma[0]<<" ; "<<sigma[1]<<")"<<std::endl
                     <<"mu    = ("<<mu[0]<<" ; "<<mu[1]<<")"<<std::endl
                     <<"nbins = "<<nbin<<std::endl;
            
            lvalue valx;
            lvalue valy;
            lweight w ;
            
            long fnmax = 10;
            while(nbin/static_cast<unsigned int>(fnmax) > 1)
                fnmax *= 10;

            double lowx = mu[0]-5.*fabs(sigma[0]);
            double upx  = mu[0]+5.*fabs(sigma[0]);
            
            double lowy = mu[1]-5.*fabs(sigma[1]);
            double upy  = mu[1]+5.*fabs(sigma[1]);
            
            for(long kx = 0; kx < static_cast<long>(nbin); kx++)
            {
                double x0 = lowx + (upx-lowx)/static_cast<double>(nbin-1)*static_cast<double>(kx);
                valx[kx] = x0;
                
                for(long ky = 0; ky < static_cast<long>(nbin); ky++)
                {
                    double y0 = lowy + (upy-lowy)/static_cast<double>(nbin-1)*static_cast<double>(ky);
                    
                    valy[ky] = y0;
                    w  [kx*fnmax+ky] = exp(-0.5*pow((x0 - mu[0])/sigma[0],2.))*exp(-0.5*pow((y0 - mu[1])/sigma[1],2.));
                }
            }

            binned_pdf2 *P = new binned_pdf2(valx, valy, w);
            
            valx.clear();
            valy.clear();
            w.clear();
            
            return *P;
        }
        
        binned_pdf2 binned_pdf2::logNormal(double mux, double muy, double sigmax, double sigmay, unsigned int nbin)
        {
            double *mu = new double[2]; mu[0]=mux; mu[1]=muy;
            double *sigma = new double[2]; sigma[0]=sigmax; sigma[1]=sigmay;
            
            binned_pdf2 G = logNormal(mu, sigma, nbin);
            delete [] mu;
            delete [] sigma;
            
            return G;
        }
        
        /**
         *  @details Generate a binned log-normal distribution
         *
         *  @param mu    Mean of the log-normal
         *  @param sigma standard deviation of the log-normal
         *  @param nbin  Number of bin use to sample the log-normal distribution
         *
         *  @return A binned probability distribution with a gaussian shape.
         */
        binned_pdf2 & binned_pdf2::logNormal(double *mu, double *sigma, unsigned int nbin)
        {
            lvalue valx;
            lvalue valy;
            lweight w ;
            
            long fnmax = 10;
            while(nbin/static_cast<unsigned int>(fnmax) > 1)
                fnmax *= 10;
            
            double lowx = (mu[0]-10.*fabs(sigma[0]) < 0)? 0 : mu[0]-10.*fabs(sigma[0]);
            double upx  =  mu[0]+10.*fabs(sigma[0]);
            
            double lowy = (mu[1]-10.*fabs(sigma[1]) < 0)? 0 : mu[1]-10.*fabs(sigma[1]);
            double upy  =  mu[1]+10.*fabs(sigma[1]);
            
            for(int kx = 0; kx < static_cast<long>(nbin); kx++)
            {
                double x0 = lowx + (upx-lowx)/static_cast<double>(nbin-1)*static_cast<double>(kx);
                
                if(x0 == 0)
                    continue;
                
                valx[kx] = x0;

                for(int ky = 0; ky < static_cast<long>(nbin); ky++)
                {
                    double y0 = lowy + (upy -lowy)/static_cast<double>(nbin-1)*static_cast<double>(ky);
                    
                    valy[ky] = y0;
                    
                    w  [kx*fnmax+ky] = 1./x0*exp(-0.5*(log(x0/mu[0])*log(x0/mu[0]))
                                         /sigma[0]/sigma[0])
                                     * 1./y0*exp(-0.5*(log(y0/mu[1])*log(y0/mu[1]))/sigma[1]/sigma[1]);
                }
            }
            
            binned_pdf2 *P = new binned_pdf2(valx, valy, w);
            
            
            valx.clear();
            valy.clear();
            w.clear();
            
            return *P;
        }
	
#pragma mark - statistic::important_sampling implementation
#pragma mark • Contructor/Destriuctor
	
	/**
	 *  @details Default constructor with an empty PDF distribution.
	 *  @param scale PDF scaling factor
	 */
	important_sampling::important_sampling(float scale): fscale(scale)
	{
		normal_pdf = NULL;
		biased_pdf = NULL;
	}
	
	/**
	 *  @details Constructs an statistic::important_sampling object based on an input statistic::binned_pdf.
	 *
	 *  @param P     Input PDF
	 *  @param scale scaling factor applied by the important sampling algorithm to artificially bias the distribution
	 */
	important_sampling::important_sampling(const binned_pdf &P, float scale):fscale(scale)
	{
		normal_pdf = new binned_pdf(P);
		biased_pdf = NULL;
	}
	
	/**
	 *  @details Copy constructor
	 *  @param P statistic::important_sampling object to copy.
	 */
	important_sampling::important_sampling(const important_sampling &P): fscale(P.fscale)
	{
		normal_pdf = new binned_pdf(*(P.normal_pdf));
		biased_pdf = new binned_pdf(*(P.biased_pdf));
	}
	
	/**
	 *  Reset content of the PDF from memory.
	 */
	void important_sampling::reset()
	{
		if(biased_pdf != NULL)
			delete biased_pdf;
		
		biased_pdf = NULL;
	}
	
#pragma mark • Randomness
	/**
	 *  @details Set the random seed
	 *  @param seed Random seed
	 */
	void important_sampling::setSeed(int64_t seed)
	{
		if(normal_pdf != NULL)
		{
			normal_pdf->setSeed(seed);
		}
	}
	
	/**
	 *  @details Retrive the random seed used by the random generator	 *
	 *  @return the random seed
	 */
	int64_t important_sampling::getSeed()
	{
		int64_t seed = 0;
		
		if(normal_pdf != NULL)
		{
			seed = normal_pdf->getSeed();
		}
		
		return seed;
	}
	
#pragma mark • PDF configuration
	/**
	 *  Modify PDF stored by in the pimath::ramdom_sampling memory with an other PDF.
	 *  @param P statistic::binned_pdf to be used to generate random number
	 */
	void important_sampling::setPDF(const binned_pdf &P)
	{
		reset();
		
		if(normal_pdf != NULL)
			delete normal_pdf;
		
		normal_pdf =  new binned_pdf(P);
	}
	
	/**
	 *  This is a private methods which biased the PDF function to artificially increases the probability of rare event to occurs. The biased correspond to an enlargment of the standard deviation of the optimal PDF :
	 \f[ 
		P'(x,\mu,\sigma') = P(x,\mu,\sigma\times f)
	 \f]
	 with \f$\mu\f$ the mean and \f$\sigma\f$ the 1 standard deviation of the PDF \f$P(x)\f$ and \f$f\f$ the PDF important-sampling scaling factor.
	 */
	void important_sampling::biasPDF()
	{
		if(biased_pdf != NULL)
		{
			delete biased_pdf;
			biased_pdf = NULL;
		}
		
		biased_pdf = new binned_pdf(*(normal_pdf));
		biased_pdf->scale(fscale);
	}
	
#pragma mark • Random generator
	/**
	 *  @details Generates random number from a biased PDF and calculate the true prtobability of the generated event.
	 *
	 *  @param w upon return true probability of the genrated event.
	 *  @return A generated number
	 *  @note Bicause the PDF is artificially biased to increase the probability of rare event, the true probability of the generated number is estimated as :
	 \f[
		w = \frac{P(x,\mu,\sigma)}{P(x,\mu,\sigma\times f)}
	 \f]
	 */
	double important_sampling::operator()(double &w)
	{
		w = 0;
		if(normal_pdf == NULL)
			return 0;
		
		if(biased_pdf == NULL)
			biasPDF();
			
		double scaledP = normal_pdf->gen()*fscale;
		while(biased_pdf->eval(scaledP) <= 0) scaledP = normal_pdf->gen()*fscale;
		
		w = normal_pdf->eval(scaledP)/biased_pdf->eval(scaledP);
		
		return scaledP;
	}
	
#pragma mark - statistic::pdf class implementation
	pdf::pdf ()
	{
		fmax = 1;
		
		pdf_norme = 1;
		
		pdf_lowEdge = 0;
		pdf_upEdge  = 0;
		
		xSlice = std::numeric_limits<double>::min();
		ySlice = std::numeric_limits<double>::min();
		
#if __cplusplus >= 199711L
		// obtain a seed from the timer
		seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	        rndgen = std::mt19937_64(seed);
#else
		seed = time(NULL);
		srand(seed);
#endif
	}
	
    void pdf::add_pdf(pdf_function ptr, pdf_param p)
	{
		fpdf.push_back(pdf_distribution(p,ptr));
	}
	
#pragma mark – PDF estimator
        
	double pdf::operator()(double x)
	{
		double w =  0;
		pdf_list::iterator it = fpdf.begin();
		
		while(it != fpdf.end())
		{
			pdf_function f = (*it).second;
			w += (this->*f)(x,(*it).first);
			
			it++;
		}
		
		return w*pdf_norme;
	}
	
#pragma mark – pdf
	void pdf::setRange(double min, double max)
	{
		if(min == max)
		{
			double step = (max-min)/100.;
			min -= step;
			max += step;
		}
		pdf_lowEdge = (min < max)?min:max;
		pdf_upEdge  = (max > min)?max:min;
		
	}
	
#pragma mark – Random kernel
	
    double pdf::gen(void)
    {
        return this->operator()();
    }
        
	double pdf::operator()(void)
	{
		
		if(pdf_upEdge <= pdf_lowEdge)
			setRange(pdf_upEdge, pdf_lowEdge);
        
#if __cplusplus >= 199711L
        std::uniform_real_distribution<double> xdis(pdf_lowEdge, pdf_upEdge);
#endif
		
		if( fabs(xSlice) <= std::numeric_limits<double>::min() )
		{
#if __cplusplus >= 199711L
			xSlice = xdis(rndgen);
#else
			xSlice = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*(pdf_upEdge - pdf_lowEdge) + pdf_lowEdge;
#endif
		}
		
#if __cplusplus >= 199711L
		std::uniform_real_distribution<double> ydis(0, (*this)(xSlice));
		ySlice = ydis(rndgen);
#else
		ySlice = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*((*this)(xSlice));
#endif

#if __cplusplus >= 199711L
		double xgen = xdis(rndgen);
		while((*this)(xgen) < ySlice || (*this)(xgen) <= 0)
			xgen = xdis(rndgen);
#else
		double xgen = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*(pdf_upEdge - pdf_lowEdge) + pdf_lowEdge;
        while((*this)(xgen) < ySlice)
	  xgen = static_cast<double>( rand() % MOD_RAND )/DMOD_RAND*(pdf_upEdge - pdf_lowEdge) + pdf_lowEdge;
#endif
		
		xSlice = xgen;
		
		fmax = (ySlice > fmax) ? ySlice : fmax;
		
		return xSlice;
	}
    }
		
}




