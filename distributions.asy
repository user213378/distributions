/*
// distributions.asy
*** Version 0.2a  ... ***

- Add several new distributions (include PDF and CDF)

// Supported on the whole real line

Will be updated later ...

*** Version 0.2  8/26/2020 ***

- Add several new distributions (include PDF and CDF)

// Supported on a bounded interval

PERT distribution
Raised cosine distribution
Reciprocal distribution
Triangular distribution
Trapezoidal distribution
Truncated normal distribution
U-quadratic distribution
Wigner semicircle distribution
Continuous Bernoulli distribution

// Supported on intervals of length 2π – directional distributions
von Mises distribution
Wrapped normal distribution
Wrapped exponential distribution
Wrapped Lévy distribution
Wrapped Cauchy distribution
Wrapped asymmetric Laplace distribution

// Supported on semi-infinite intervals, usually [0,+infinity)
Inverse Gaussian distribution
Lévy distribution
Log-Cauchy distribution
Log-Laplace distribution
Log-logistic distribution
Lomax distribution
Mittag-Leffler distribution
Nakagami distribution
Rice distribution
Shifted Gompertz distribution

*** Version 0.1 8/23-24/2020 ***

- Rearrange the old distributions (NEW)
- Add several new distributions (include PDF and CDF)

The contents:

// Continuous distributions

// Supported on a bounded interval

- Arcsine distribution
- The Beta distribution
- The Logit normal Distribution
- The Flat (Uniform) Distribution
- Irwin–Hall distribution
- The Bates distribution
- Kumaraswamy distribution

// Supported on semi-infinite intervals, usually [0,+infinity)

- The Beta prime distribution
- Birnbaum–Saunders distribution (Fatigue Life Distribution)
- Chi distribution
- The Chi-squared distribution
- Dagum distribution
- The Exponential distribution
- Exponential-logarithmic distribution
- The F-distribution
- Folded normal distribution
- Fréchet distribution
- The Gamma distribution
- Generalized gamma distribution
- Generalized Pareto distribution
- Gamma/Gompertz distribution
- Gompertz distribution
- Half-normal distribution
- The Lognormal Distribution
- The Pareto (Type I) Distribution
- The Rayleigh distribution
- Type-2 Gumbel distribution
- The Weibull Distribution

// Supported on the whole real line

- The Cauchy distribution
- The Generalized Normal (version 1)
- The Generalized Normal (version 2)
- Gumbel distribution
- The Landau Distribution
- The Laplace distribution
- The Logistic distribution
- The Normal distribution (Gauss)
- The Student's t-distribution
- Type-1 Gumbel distribution

*** Version 0.0 8/21-22/2020 ***

This package has the following distributions (in gsl.cc): (include PDF and CDF)
- The Normal distribution (Gauss)
- The Cauchy distribution
- The Exponential distribution
- The Laplace distribution
- The Generalized Normal (version 1)
- The Generalized Normal (version 2)
- The Gamma distribution
- The Rayleigh distribution
- The Chi-squared distribution
- The Logistic distribution
- The Student's t-distribution
- The Beta distribution
- The Weibull Distribution
- The Pareto (Type I) Distribution
- The Flat (Uniform) Distribution
- The Lognormal Distribution
- The F-distribution
- Type-1 Gumbel distribution
- Type-2 Gumbel distribution
- Gumbel distribution
- Fréchet distribution
- The Landau Distribution

*/


private import gsl;

typedef real REAL(real);
typedef pair PAIR(real);

// https://en.wikipedia.org/wiki/List_of_probability_distributions

/* Continuous distributions */

// Supported on a bounded interval

/* Arcsine distribution */

REAL ArcsinedistributionPDF(real a=0, real b=1){
  return new real(real x){
    return 1/(pi*sqrt((x-a)*(b-x)));
  };
};
REAL ArcsinedistributionCDF(real a=0, real b=1){
  return new real(real x){
    return 2/pi*asin(sqrt((x-a)/(b-a)));
  };
};

/* The Beta distribution */

REAL BetadistributionPDF(real alpha=1, real beta=1){
  return new real(real x){
    return pdf_beta(x,a=alpha,b=beta);
  };
};
REAL BetadistributionCDF(real alpha=1, real beta=1){
  return new real(real x){
    return cdf_beta_P(x,a=alpha,b=beta);
  };
};

/* The Logit normal Distribution */
// https://en.wikipedia.org/wiki/Logit-normal_distribution

REAL LogitnormaldistributionPDF(real mu=0, real sigma=1){
  return new real(real x){
    real a=1/(sqrt(2*pi)*sigma*x*(1-x));
	real b=log(x/(1-x));
    return a*exp(-0.5*((b-mu)/sigma)^2);
  };
};
REAL LogitnormaldistributionCDF(real mu=0, real sigma=1){
  return new real(real x){
    real b=log(x/(1-x));
    return 0.5*(1+erf((b-mu)/sqrt(2*sigma*sigma)));
  };
};

/* The Flat (Uniform) Distribution */

REAL UniformdistributionPDF(real a=0, real b=1){
  return new real(real x){
    return pdf_flat(x,a,b);
  };
};
REAL UniformdistributionCDF(real a=0, real b=1){
  return new real(real x){
    return cdf_flat_P(x,a,b);
  };
};

/* Irwin–Hall distribution */
// https://en.wikipedia.org/wiki/Irwin%E2%80%93Hall_distribution

REAL IrwinHalldistributionPDF(int n=1){
  return new real(real x){
    real a=1/factorial(n-1);
    real s=0;
	for (int k=0; k <= x; ++k){
	  s=s+(-1)^k*choose(n,k)*(x-k)^(n-1);
	}	
    return a*s;
  };
};

REAL IrwinHalldistCDF(int n=1){
  return new real(real x){
    real a=1/factorial(n);
    real s=0;
	for (int k=0; k <= x; ++k){
	  s=s+(-1)^k*choose(n,k)*(x-k)^n;
	}	
    return a*s;
  };
};

PAIR IrwinHalldistributionCDF(int n=1){
  return new pair(real x){
    return (x <= 0) ? (x,0) : (x >= n) ? (x,1) : (x,IrwinHalldistCDF(n)(x));
  };
};

/* The Bates distribution */

REAL BatesdistributionPDF(int n=1, real a=0, real b=1){
  return new real(real x){
	real s=0;
	for (int k=0; k <= n*x ; ++k){
	  // https://commons.wikimedia.org/wiki/File:BatesPDF.svg
	  real X=n*x-k;
	  real m = (X < 0) ? -1 : ( X == 0) ? 0 : 1;
	  s=s+(-1)^k*choose(n,k)*X^(n-1)*m;
	}	
    return s*n/factorial(n-1);
  };
};

REAL BatesdistributionCDF(int n=1, real a=0, real b=1){
  return new real(real x){
	real s=0;
	for (int k=0; k <= n*x ; ++k){
	  // https://commons.wikimedia.org/wiki/File:BatesCDF.svg
	  real X=n*x-k;
	  real m = (X < 0) ? -1 : ( X == 0) ? 0 : 1;
	  s=s+(-1)^k*choose(n,k)*X^(n)*m;
	}	
    return s/factorial(n);
  };
};

/* Kumaraswamy distribution */

REAL KumaraswamydistributionPDF(real a=0.5, real b=0.5){
  return new real(real x){
    return a*b*x^(a-1)*(1-x^a)^(b-1);
  };
};

REAL KumaraswamydistributionCDF(real a=0.5, real b=0.5){
  return new real(real x){
    return 1-(1-x^a)^b;
  };
};

/* PERT distribution */

REAL PERTdistributionPDF(real a=0, real b=10, real c=100){
  return new real(real x){
    real alpha=1+4*(b-a)/(c-a);
	real beta=1+4*(c-b)/(c-a);
    real A=(x-a)^(alpha-1)*(c-x)^(beta-1);
	real B=beta(alpha,beta)*(c-a)^(alpha+beta-1);
    return A/B;
  };
};

REAL PERTdistributionCDF(real a=0, real b=10, real c=100){
  return new real(real x){
    real alpha=1+4*(b-a)/(c-a);
	real beta=1+4*(c-b)/(c-a);  
	real z=(x-a)/(c-a);
    return beta(alpha,beta,z)/beta(alpha,beta);
  };
};

/* Raised cosine distribution */

REAL RaisedcosinedistributionPDF(real mu=0, real s=0.5){
  return new real(real x){
    real A=pi*(x-mu)/s;
    return (1+cos(A))/(2*s);
  };
};

REAL RaisedcosinedistributionCDF(real mu=0, real s=0.5){
  return new real(real x){
    real A=(x-mu)/s;
    return 0.5*(1+A+sin(A*pi)/pi);
  };
};

/* Reciprocal distribution */

REAL ReciprocaldistributionPDF(real a=0.5, real b=1){
  return new real(real x){
    return 1/(x*log(b/a));
  };
};

REAL ReciprocaldistributionCDF(real a=0.5, real b=1){
  return new real(real x){
    return log(x/a)/log(b/a);
  };
};

/* Triangular distribution */
// https://en.wikipedia.org/wiki/Triangular_distribution

REAL TriangulardistPDF(real a=0.25, real c=0.5, real b=1){
  return new real(real x){
    if (x < a) { return 0;}
	  else if ( x >= a && x < c){ return 2*(x-a)/((b-a)*(c-a));}
	    else if (x == c) { return 2/(b-a);}
		  else if ( x > c && x <= b) { return 2*(b-x)/((b-a)*(b-c));}
            else { return 0;}
  };
};

REAL TriangulardistCDF(real a=0.25, real c=0.5, real b=1){
  return new real(real x){
    if ( x <= a) { return 0;}
	  else if ( x > a && x <= c) { return (x-a)^2/((b-a)*(c-a));}
	    else if ( x > c && x < b) { return 1-(b-x)^2/((b-a)*(b-c));}
		  else { return 1;}
  };
};

PAIR TriangulardistributionPDF(real a=0.25, real c=0.5, real b=1){
  return new pair(real x){
    if (x < a || x > b) { return (x,0);}
	  else { return (x,TriangulardistPDF(a,c,b)(x));}
  };
};

PAIR TriangulardistributionCDF(real a=0.25, real c=0.5, real b=1){
  return new pair(real x){
    if ( x <= a) { return (x,0);}
	  else if ( x >= b) { return (x,1);}
	    else { return (x,TriangulardistCDF(a,c,b)(x));}
  };
};

/* Trapezoidal distribution */
// https://en.wikipedia.org/wiki/Trapezoidal_distribution

REAL TrapezoidaldistPDF(real a=0, real b=0.25, real c=0.5, real d=1){
  return new real(real x){
    real H=2/(d+c-a-b);
    if (x >= a && x < b) { return H*(x-a)/(b-a);}
	  else if ( x >= b && x < c){ return H;}
        else { return H*(d-x)/(d-c);}
  };
};

REAL TrapezoidaldistCDF(real a=0, real b=0.25, real c=0.5, real d=1){
  return new real(real x){
    real H=1/(d+c-a-b);
    if (x >= a && x < b) { return H*(x-a)^2/(b-a);}
	  else if ( x >= b && x < c){ return H*(2*x-a-b);}
        else { return 1-H*(d-x)^2/(d-c);}
  };
};

PAIR TrapezoidaldistributionPDF(real a=0, real b=0.25, real c=0.5, real d=1){
  return new pair(real x){
    return (x,TrapezoidaldistPDF(a,b,c,d)(x));
  };
};

PAIR TrapezoidaldistributionCDF(real a=0, real b=0.25, real c=0.5, real d=1){
  return new pair(real x){
	return (x,TrapezoidaldistCDF(a,b,c,d)(x));
  };
};

/* Truncated normal distribution */
// https://en.wikipedia.org/wiki/Truncated_normal_distribution

REAL TruncatednormaldistributionPDF(real mu=0, real sigma=1, real a=0, real b=1){
  return new real(real x){
    real A=(x-mu)/sigma;
	real B=(b-mu)/sigma;
	real C=(a-mu)/sigma;
    return (1/sigma)*pdf_gaussian(A,mu=0,1)/(cdf_gaussian_P(B,mu=0,1)-cdf_gaussian_P(C,mu=0,1));
  };
};

REAL TruncatednormaldistributionCDF(real mu=0, real sigma=1, real a=0, real b=1){
  return new real(real x){
    real A=(x-mu)/sigma;
	real B=(b-mu)/sigma;
	real C=(a-mu)/sigma;
    real Z=cdf_gaussian_P(B,mu=0,1)-cdf_gaussian_P(C,mu=0,1);
    return (cdf_gaussian_P(A,mu=0,1)-cdf_gaussian_P(C,mu=0,1))/Z;
  };
};

/* U-quadratic distribution */

REAL UquadraticdistributionPDF(real alpha=0.5, real beta=1){
  return new real(real x){
    return alpha*(x-beta)^2;
  };
};

REAL UquadraticdistributionCDF(real alpha=0.5, real beta=1, real a=0.5){
  return new real(real x){
    return (alpha/3)*((x-beta)^3+(beta-a)^3);
  };
};

/* Wigner semicircle distribution */

REAL WignersemicircledistributionPDF(real R=0.5){
  return new real(real x){
    return 2*sqrt(R^2-x^2)/(pi*R^2);
  };
};

REAL WignersemicircledistributionCDF(real R=0.5){
  return new real(real x){
    return 0.5+x*sqrt(R^2-x^2)/(pi*R^2)+asin(x/R)/pi;
  };
};

/* Continuous Bernoulli distribution */

REAL ContinuousBernoullidistributionPDF(real lambda=0.5){
  return new real(real x){
    real C(real a){
	  return (a != 1/2) ? 2*atanh(1-2*a)/(1-2*a) : 2;
	}
    return C(lambda)*(lambda^x)*(1-lambda)^(1-x);
  };
};

REAL ContinuousBernoullidistributionCDF(real lambda=0.5){
  return new real(real x){
    real A= ((lambda^x)*(1-lambda)^(1-x)+lambda-1)/(2*lambda-1);
    return (lambda != 0.5) ? A : x;
  };
};

// Supported on intervals of length 2π – directional distributions

/* von Mises distribution */

REAL vonMisesdistributionPDF(real mu=0, real kappa=1){
  return new real(real x){
    return exp(kappa*cos(x-mu))/(2*pi*I0(kappa));
  };
};
/* I don't know how to write.
REAL vonMisesdistributionCDF(real alpha=1, real beta=1){
  return new real(real x){
    return beta(a=alpha,b=beta,x/(1+x));
  };
};
*/

/* Wrapped normal distribution */

// It is difficult for me.

/* Wrapped exponential distribution */

REAL WrappedexponentialdistributionPDF(real lambda=0.5){
  return new real(real x){
    return lambda*exp(-lambda*x)/(1-exp(-2*pi*lambda));
  };
};

REAL WrappedexponentialdistributionCDF(real lambda=0.5){
  return new real(real x){
    return (1-exp(-lambda*x))/(1-exp(-2*pi*lambda));
  };
};

/* Wrapped Lévy distribution */

// It is difficult for me.

/* Wrapped Cauchy distribution */

REAL WrappedCauchydistributionPDF(real mu=0, real gamma=0.5){
  return new real(real x){
    real A=(sinh(gamma))/(cosh(gamma)-cos(x-mu));
    return A/(2*pi);
  };
};

// I don't know its CDF

/* Wrapped asymmetric Laplace distribution */
// have only PDF

REAL WrappedasymmetricLaplacedistributionPDF(real m=0, real lambda=0.5, real kappa=.5){
  return new real(real x){
    real A=exp(-(x-m)*lambda*kappa);
	real B=exp((x-m)*lambda/kappa);
    return (x >= m) ? A/(1-exp(-2*pi*kappa*lambda))-B/(1-exp(2*pi*lambda/kappa)) 
	                : A/(exp(2*pi*kappa*lambda)-1)-B/(exp(-2*pi*lambda/kappa)-1) ;
  };
};


// Supported on semi-infinite intervals, usually [0,+infinity)

/* The Beta prime distribution */

REAL BetaprimedistributionPDF(real alpha=1, real beta=1){
  return new real(real x){
    return (x^(alpha-1)*(1+x)^(-alpha-beta))/beta(alpha,beta);
  };
};
REAL BetaprimedistributionCDF(real alpha=1, real beta=1){
  return new real(real x){
    return beta(a=alpha,b=beta,x/(1+x));
  };
};

/* Birnbaum–Saunders distribution (Fatigue Life Distribution) */

REAL BirnbaumSaundersdistributionPDF(real mu=0, real gamma=1, real beta=1){
  return new real(real x){
    // https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution
	// http://atomic.phys.uni-sofia.bg/local/nist-e-handbook/e-handbook/eda/section3/eda366a.htm
    real A=sqrt((x-mu)/beta);
	real B=1/(2*gamma*(x-mu));
	return (A+1/A)*B*pdf_gaussian((A-1/A)/gamma,mu=0,sigma=1);
  };
};

REAL BirnbaumSaundersdistributionCDF(real gamma=1){
  return new real(real x){
    return cdf_gaussian_P((sqrt(x)-sqrt(1/x))/gamma,mu=0,sigma=1);
  };
};

/* Chi distribution */

REAL ChidistributionPDF(real k=0.5){
  return new real(real x){
    real A=x^(k-1)*exp(-x^2/2);
	real B=2^(k/2-1)*gamma(k/2);
    return A/B;
  };
};

REAL ChidistributionCDF(real k=0.5){
  return new real(real x){
  // http://www.ntrand.com/chi-distribution/
    return gamma(k/2,x^2/2)/gamma(k/2);
  };
};

/* The Chi-squared distribution */

REAL ChisquaredistributionPDF(real k=1){
  return new real(real x){
    return pdf_chisq(x,k);
  };
};
REAL ChisquaredistributionCDF(real k=1){
  return new real(real x){
    return cdf_chisq_P(x,k);
  };
};

/* Dagum distribution */

REAL DagumdistributionPDF(real p=0.5, real a= 0.5, real b=0.5){
  return new real(real x){
    // https://en.wikipedia.org/wiki/Dagum_distribution
    real A=a*p/x;
	real B=(x/b)^(a*p);
	real C=((x/b)^a+1)^(p+1);
    return A*B/C;
  };
};
REAL DagumdistributionCDF(real p=0.5, real a= 0.5, real b=0.5){
  return new real(real x){
    return (1+(x/b)^(-a))^(-p);
  };
};

/* The Exponential distribution */

REAL ExponentialdistributionPDF(real lambda=1){
  return new real(real x){
    return pdf_exponential(x,1/lambda);
  };
};
REAL ExponentialdistributionCDF(real lambda=1){
  return new real(real x){
    return cdf_exponential_P(x,1/lambda);
  };
};

/* Exponential-logarithmic distribution */
// https://en.wikipedia.org/wiki/Exponential-logarithmic_distribution

REAL ELdistributionPDF(real p=0.5, real beta=0.5){
  return new real(real x){
    real A=1/(-log(p));
	real B=beta*(1-p)*exp(-beta*x);
	real C=1-(1-p)*exp(-beta*x);
    return A*B/C;
  };
};
REAL ELdistributionCDF(real p=0.5, real beta=0.5){
  return new real(real x){
    real A=log(1-(1-p)*exp(-beta*x));
    return 1-A/log(p);
  };
};

/* The F-distribution */

REAL FdistributionPDF(real d1=1, real d2=1){
  return new real(real x){
    return pdf_fdist(x,nu1=d1,nu2=d2);
  };
};
REAL FdistributionCDF(real d1=1, real d2=1){
  return new real(real x){
    return cdf_fdist_P(x,nu1=d1,nu2=d2);
  };
};

/* Folded normal distribution */
// https://en.wikipedia.org/wiki/Folded_normal_distribution

REAL FoldednormaldistributionPDF(real mu=0, real sigma=0.5){
  return new real(real x){
    real a=(x-mu)^2/(2*sigma*sigma);
	real b=(x+mu)^2/(2*sigma*sigma);
	real q=1/(sigma*sqrt(2*pi));
    return q*exp(-a)+q*exp(-b);
  };
};
REAL FoldednormaldistributionCDF(real mu=0, real sigma=0.5){
  return new real(real x){
    real a=(x+mu)/(sigma*sqrt(2));
	real b=(x-mu)/(sigma*sqrt(2));
    return 0.5*(erf(a)+erf(b));
  };
};

/* Fréchet distribution */

REAL FrechetdistributionPDF(real alpha=0.5, real s=1, real m=0){
  return new real(real x){
  // http://www.mathwave.com/articles/extreme-value-distributions.html#Fr%C3%A9chetDistribution
    real z=s/(x-m); // x > m
    return alpha/s*z^(alpha+1)*exp(-z^alpha); 
  };
};

REAL FrechetdistributionCDF(real alpha=0.5, real s=1, real m=0){
  return new real(real x){
    real z=s/(x-m); // x > m
    return exp(-z^alpha);
  };
};

/* The Gamma distribution */

REAL GammadistributionPDF(real k=1, real theta=1){
  return new real(real x){
    // scale parameter 'theta > 0'
	// exponent (or shape parameter) 'k > 0'
    return pdf_gamma(x,a=k,b=theta);
  };
};
REAL GammadistributionCDF(real k=1, real theta=1){
  return new real(real x){
    return cdf_gamma_P(x,a=k,b=theta);
  };
};

/* Generalized gamma distribution */
// https://en.wikipedia.org/wiki/Generalized_gamma_distribution

REAL GeneralizedgammadistributionPDF(real a=0.5, real d=0.5, real p=0.5){
  return new real(real x){
    real A=(p/(a^d))/gamma(d/p);
	real B=x^(d-1)*exp(-(x/a)^p);
    return A*B;
  };
};

REAL GeneralizedgammadistributionCDF(real a=0.5, real d=0.5, real p=0.5){
  return new real(real x){
  // https://mathworld.wolfram.com/IncompleteGammaFunction.html
    real A=gamma(d/p)-gamma(d/p,(x/a)^p);
	real B=gamma(d/p);
    return A/B;
  };
};

/* Generalized Pareto distribution */
// https://en.wikipedia.org/wiki/Generalized_Pareto_distribution

REAL GeneralizedParetodistributionPDF(real mu=0, real sigma=0.5, real ksi=0.5){
  return new real(real x){
    real a=(x-mu)/sigma;
    return (1/sigma)*(1+ksi*a)^(-1/ksi-1);
  };
};

REAL GeneralizedParetodistributionCDF(real mu=0, real sigma=0.5, real ksi=0.5){
  return new real(real x){
    real a=(x-mu)/sigma;
    return (ksi != 0) ? 1-(1+ksi*a)^(-1/ksi) : 1-exp(-a);
  };
};

/* Gamma/Gompertz distribution */

REAL GammaGompertzdistributionPDF(real b=0.5, real s=0.5, real beta=0.5){
  return new real(real x){
    real a=(beta-1+exp(b*x))^(s+1);
    return (b*s*exp(b*x)*beta^s)/a;
  };
};

REAL GammaGompertzdistributionCDF(real b=0.5, real s=0.5, real beta=0.5){
  return new real(real x){
    real a=(beta-1+exp(b*x))^s;
    return 1-(beta^s)/a;
  };
};

/* Gompertz distribution */
// https://en.wikipedia.org/wiki/Gompertz_distribution

REAL GompertzdistributionPDF(real eta=0.5, real b=0.5){
  return new real(real x){
    return b*eta*exp(eta+b*x-eta*exp(b*x));
  };
};

REAL GompertzdistributionCDF(real eta=0.5, real b=0.5){
  return new real(real x){
    return 1-exp(-eta*(exp(b*x)-1));
  };
};

/* Half-normal distribution */
// https://en.wikipedia.org/wiki/Half-normal_distribution

REAL HalfnormaldistributionPDF(real sigma=0.5){
  return new real(real x){
    return sqrt(2)/(sigma*sqrt(pi))*exp(-x^2/(2*sigma*sigma));
  };
};

REAL HalfnormaldistributionCDF(real sigma=0.5){
  return new real(real x){
    return erf(x/(sigma*sqrt(2)));
  };
};

/* The Lognormal Distribution */

REAL LognormaldistributionPDF(real mu=0, real sigma=1){
  return new real(real x){
    return pdf_lognormal(x,zeta=mu,sigma);
  };
};
REAL LognormaldistributionCDF(real mu=0, real sigma=1){
  return new real(real x){
    return cdf_lognormal_P(x,zeta=mu,sigma);
  };
};

/* The Pareto (Type I) Distribution */

REAL ParetoTypeIdistributionPDF(real alpha=1, real x_m=1){
  return new real(real x){
    return pdf_pareto(x,a=alpha,b=x_m);
  };
};
REAL ParetoTypeIdistributionCDF(real alpha=1, real x_m=1){
  return new real(real x){
    return cdf_pareto_P(x,a=alpha,b=x_m);
  };
};

/* The Rayleigh distribution */

REAL RayleighdistributionPDF(real sigma=1){
  return new real(real x){
    return pdf_rayleigh(x,sigma);
  };
};
REAL CauchydistributionCDF(real sigma=1){
  return new real(real x){
    return cdf_rayleigh_P(x,sigma);
  };
};

/* Type-2 Gumbel distribution */

REAL Type2GumbeldistributionPDF(real a=1, real b=1){
  return new real(real x){
    return pdf_gumbel2(x,a,b);
  };
};
REAL Type2GumbeldistributionCDF(real a=1, real b=1){
  return new real(real x){
    return cdf_gumbel2_P(x,a,b);
  };
};

/* The Weibull Distribution */

REAL WeibulldistributionPDF(real lambda=1, real k=1){
  return new real(real x){
    return pdf_weibull(x,a=lambda,b=k);
  };
};
REAL WeibulldistributionCDF(real lambda=1, real k=1){
  return new real(real x){
    return cdf_weibull_P(x,a=lambda,b=k);
  };
};

/* Inverse Gaussian distribution */

REAL InverseGaussiandistributionPDF(real mu=1, real lambda=1){
  return new real(real x){
    real a = -(lambda*(x-mu)^2)/(2*x*mu^2);
    return sqrt(lambda/(2*pi*x^3))*exp(a);
  };
};

REAL InverseGaussiandistributionCDF(real mu=1, real lambda=1){
  return new real(real x){
    real A = sqrt(lambda/x)*(x/mu-1);
	real B = -sqrt(lambda/x)*(x/mu+1);
    return cdf_gaussian_P(A,mu=0,sigma=1)+exp(2*lambda/mu)*cdf_gaussian_P(B,mu=0,sigma=1);
  };
};

/* Lévy distribution */

REAL LevydistributionPDF(real mu=0, real c=1){
  return new real(real x){
    real a = c/(2*(x-mu));
    return sqrt(c/(2*pi))*exp(-a)/((x-mu)^(3/2));
  };
};

REAL LevydistributionCDF(real mu=0, real c=1){
  return new real(real x){
    real a = c/(2*(x-mu));
    return erfc(sqrt(a));
  };
};

/* Log-Cauchy distribution */

REAL LogCauchydistributionPDF(real mu=0, real sigma=1){
  return new real(real x){
    real a = sigma/((log(x)-mu)^2+sigma^2);
    return a/(x*pi);
  };
};

REAL LogCauchydistributionCDF(real mu=0, real sigma=1){
  return new real(real x){
    real a = (log(x)-mu)/sigma;
    return 0.5+atan(a)/pi;
  };
};

/* Log-Laplace distribution */
// a log-laplace growth rate model pdf (google)
/*
REAL LogLaplacedistributionPDF(real mu=0, real b=1){
  return new real(real x){
    real a = -abs(log(x)-mu)/b;
    return exp(a)/(2*b*x);
  };
};

REAL LogLaplacedistributionCDF(real mu=0, real b=1){
  return new real(real x){
    return ;
  };
};
*/
/* Log-logistic distribution */

REAL LoglogisticdistributionPDF(real alpha=1, real beta=1){
  return new real(real x){
    real a = (beta/alpha)*(x/alpha)^(beta-1);
	real b = (1+(x/alpha)^beta)^2;
    return a/b;
  };
};

REAL LoglogisticdistributionCDF(real alpha=1, real beta=1){
  return new real(real x){
    return 1/(1+(x/alpha)^(-beta));
  };
};

/* Lomax distribution (Pareto Type II distribution) */

REAL LomaxdistributionPDF(real alpha=1, real lambda=1){
  return new real(real x){
    real a=1+x/lambda;
    return (alpha/lambda)*a^(-alpha-1);
  };
};

REAL LomaxdistributionCDF(real alpha=1, real lambda=1){
  return new real(real x){
    real a=1+x/lambda;
    return 1-a^(-alpha);
  };
};

/* Mittag-Leffler distribution */
// I don't know

/* Nakagami distribution */

REAL NakagamidistributionPDF(real m=1, real omega=1){
  return new real(real x){
    real a=-m*x^2/omega;
	real b=2*m^m/(gamma(m)*omega^m);
    return b*x^(2*m-1)*exp(a);
  };
};

REAL NakagamidistributionCDF(real m=1, real omega=1){
  return new real(real x){
    real a=gamma(m)-gamma(m,m*(x^2/omega)); // See GeneralizedgammadistributionCDF
    return a/gamma(m);
  };
};

// Rice distribution

REAL RicedistributionPDF(real ny=0, real sigma=1){
  return new real(real x){
    real A=-(x^2+ny^2)/(2*sigma*sigma);
	real B= J0(x*ny/(sigma^2));
    return (x/(sigma^2))*exp(A)*B;
  };
};
/* I can't write Marcum Q-function in Asymptote
REAL RicedistributionCDF(real ny=0, real sigma=1){
  return new real(real x){
    
    return a/gamma(m);
  };
};
*/

// Shifted Gompertz distribution

REAL ShiftedGompertzdistributionPDF(real b=1, real eta=1){
  return new real(real x){
    real a=exp(-b*x);
    return b*a*exp(-eta*a)*(1+eta(1-a));
  };
};

REAL ShiftedGompertzdistributionCDF(real b=1, real eta=1){
  return new real(real x){
    real a=exp(-b*x);
    return (1-a)*exp(-eta*a);
  };
};

// Supported on the whole real line

/* The Cauchy distribution */

REAL CauchydistributionPDF(real xzero=0, real gamma=1){
  return new real(real x){
    return pdf_cauchy(x-xzero,gamma);
  };
};
REAL CauchydistributionCDF(real xzero=0, real gamma=1){
  return new real(real x){
    return cdf_cauchy_P(x-xzero,gamma);
  };
};

/* The Generalized Normal (version 1) */

REAL ExponentialpowerdistributionPDF(real mu=0,real alpha=1, real beta=1){
  return new real(real x){
    // scale parameter 'alpha > 0'
	// exponent (or shape parameter) 'beta'
    return pdf_exppow(x-mu,a=alpha,b=beta);
  };
};
REAL ExponentialpowerdistributionCDF(real mu=0,real alpha=1, real beta=1){
  return new real(real x){
    return cdf_exppow_P(x-mu,a=alpha,b=beta);
  };
};

/* The Generalized Normal (version 2) */

REAL GeneralizedNormalV2PDF(real ksi=0, real alpha=1, real kappa=0.5){
  return new real(real x){
    real b=(x-ksi)/alpha;
	real a = (kappa != 0) ? -log(1-kappa*b)/kappa : b;
    real q = alpha-kappa*(x-ksi);
    return pdf_gaussian(a,mu=0,sigma=1)/q;
  };
};

REAL GeneralizedNormalV2CDF(real ksi=0, real alpha=1, real kappa=0.5){
  return new real(real x){
    real b = (x-ksi)/alpha;
    real a = (kappa != 0) ? -log(1-kappa*b)/kappa : b;
    return cdf_gaussian_P(a,mu=0,sigma=1);
  };
};

/* Gumbel distribution */

REAL GumbeldistributionPDF(real mu=0.5, real beta=1){
  return new real(real x){
    real z=(x-mu)/beta;
    return exp(-(z+exp(-z)))/beta;
  };
};

REAL GumbeldistributionCDF(real mu=0.5, real beta=1){
  return new real(real x){
    real z=(x-mu)/beta;
    return exp(-exp(-z));
  };
};

/* The Landau Distribution */

real LandaudistributionPDF(real x){
    return pdf_landau(x);
};

/* The Laplace distribution */

REAL LaplacedistributionPDF(real mu=0, real b=1){
  return new real(real x){
    return pdf_laplace(x-mu,b);
  };
};
REAL LaplacedistributionCDF(real mu=0, real b=1){
  return new real(real x){
    return cdf_laplace_P(x-mu,b);
  };
};

/* The Logistic distribution */

REAL LogisticdistributionPDF(real mu=0, real s=1){
  return new real(real x){
    return pdf_logistic(x-mu,s);
  };
};
REAL LogisticdistributionCDF(real mu=0, real s=1){
  return new real(real x){
    return cdf_logistic_P(x-mu,s);
  };
};


/* The Normal distribution (Gauss) */

REAL NormaldistributionPDF(real mu=0, real sigma=1){
  return new real(real x){
    return pdf_gaussian(x,mu=mu,sigma=sigma);
  };
};

REAL NormaldistributionCDF(real mu=0, real sigma=1){
  return new real(real x){
    return cdf_gaussian_P(x,mu=mu,sigma=sigma);
  };
};

/* The Student's t-distribution */

REAL tdistributionPDF(real ny=1){
  return new real(real x){
    return pdf_tdist(x,ny);
  };
};
REAL tdistributionCDF(real ny=1){
  return new real(real x){
    return cdf_tdist_P(x,ny);
  };
};

/* Type-1 Gumbel distribution */

REAL Type1GumbeldistributionPDF(real a=1, real b=1){
  return new real(real x){
    return pdf_gumbel1(x,a,b);
  };
};
REAL Type1GumbeldistributionCDF(real a=1, real b=1){
  return new real(real x){
    return cdf_gumbel1_P(x,a,b);
  };
};
