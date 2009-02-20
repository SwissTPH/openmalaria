#include "infection.h"

class EmpiricalWithinHostModel : public Infection {
public:
  EmpiricalWithinHostModel ();
  double getNewDensity(double * transformedLaggedDensities, int ageOfInfection);
  void initialiseInfection(double * transfomedLaggedDensities);
  void setInflationFactors(double inflationMean, double inflationVariance);
private:
  double getInflatedDensity(double nonInflatedDensity);
  double sigma_noise(int ageOfInfection);
  double samplePatentValue(double mu, double sigma, double lowerBound);
  double sampleSubPatentValue(double mu, double sigma, double upperBound);
  double inverseBoxCoxTransform(double transformedValue);
  double boxCoxTransform(double untransformedValue);
  static const int _maximumDurationInDays=418; 
  double _maximumPermittedAmplificationPerCycle;
  double _subPatentLimit;
  double _inflationVariance;
  double _inflationMean;
  double _lambda;
  double _alpha1;
	double _alpha2;	
  double _alpha3;
  double _sigma_alpha1;	
	double _sigma_alpha2;	
	double _sigma_alpha3;
  double _sigma0_res;	
  double _sigmat_res;
  double _mu_beta1[_maximumDurationInDays];
  double _sigma_beta1[_maximumDurationInDays];
  double _mu_beta2[_maximumDurationInDays];
  double _sigma_beta2[_maximumDurationInDays];
  double _mu_beta3[_maximumDurationInDays];
  double _sigma_beta3[_maximumDurationInDays];  
};