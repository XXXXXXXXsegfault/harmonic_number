#include <stdio.h>

// number of terms to approximate 1/(1-n)
#define TERMS 100

// storing expansions of (n-k)^m, which are used to generate the sum formula

// Note that all polynomial coefficients in this program are stored in ascending
// order.
double polynomial_diff[TERMS+1][TERMS+2];

// calculate expansions of (n-k)^m for m = 1 ... TERMS+1 with a given k
void calculate_polynomial_diff(double k)
{
	int i,j;
	// m = 1, (n-k)^m=n-k
	polynomial_diff[0][0]=-k;
	polynomial_diff[0][1]=1.0;
	// m = 2, 3, 4, ..., TERMS+1
	// (n-k)^m=(n-k)^(m-1)*(n-k)
	for(i=1;i<TERMS+1;++i)
	{
		polynomial_diff[i][0]=-polynomial_diff[i-1][0]*k;
		for(j=1;j<=i;++j)
		{
			polynomial_diff[i][j]=polynomial_diff[i-1][j-1]-polynomial_diff[i-1][j]*k;
		}
		polynomial_diff[i][i+1]=1.0;
	}
}


// This function accepts a polynomial P(n) ("polynomial") and a number k, and returns
// the polynomial p(n) ("sum_formula") so that p(n) - p(n-k) = P(n).

// The resulting p(n) can be used to calculate the sum P(n) + P(n+k) +
// P(n+2*k) + ... + P(n+m*k). When p(n) - p(n-k) = P(n), the sum above equals
// to p(n+m*k) - p(n-k).

// Note: This function will modify the contents of "polynomial". "sum_formula[0]" stores
// the coefficient of n^1.
void generate_sum_formula(double *polynomial,double k,double *sum_formula)
{
	int i,j;
	double scale;
	calculate_polynomial_diff(k);
	for(i=TERMS;i>=0;--i)
	{
		scale=polynomial[i]/polynomial_diff[i][i];
		for(j=0;j<i;++j)
		{
			polynomial[j]-=scale*polynomial_diff[i][j];
		}
		sum_formula[i]=-scale;
	}
}


// This function approximates the sum 1/a + 1/(a+1) + 1/(a+2) + ... + 1/(b-1) + 1/b.

// The sum above equals to 1/x * (1/(1-(x-a)/x) + 1/(1-(x-a-1)/x) + 1/(1-(x-a-2)/x) +
// ... + 1/(1-(x-b+1)/x) + 1/(1-(x-b)/x)) (S).
// We can approximate 1/(1-n) using 1 + n + n^2 + n^3 + ... +n^m when -1 < n < 1. If we do
// so, the error is E=n^(m+1)/(1-n). In the sum (S), if we perform this approximation,
// the error of sum (S) will not exceed E*(b-a)/x.
// To minimize the error, a should not be less than b/2.
double sum_of_fractions(unsigned long long a,unsigned long long b)
{
	double x=((double)a+(double)b+1.0)*0.5;
	double n=((double)b+1.0-x)/x;
	double n1=n,n2=-n,result=0.0;
	int i;
	double polynomial[TERMS+1],sum_formula[TERMS+1];
	for(i=0;i<TERMS+1;++i)
	{
		polynomial[i]=1.0;
	}
	generate_sum_formula(polynomial,1.0/x,sum_formula);
	// When m is even, n^m-(-n)^m=0.
	for(i=0;i<TERMS+1;i+=2)
	{
		result+=(n1-n2)*sum_formula[i];
		n1*=n*n;
		n2*=n*n;
	}
	return result/x;
}

int main(void)
{
	unsigned long long n;
	unsigned long long a=16384,b=32768;
	double result=0.0;
	unsigned long long i;
	printf("Input an integer: ");
	scanf("%llu",&n);
	if(n<16384)
	{
		for(i=1;i<=n;++i)
		{
			result+=1.0/(double)i;
		}
	}
	else
	{
		for(i=1;i<=16383;++i)
		{
			result+=1.0/(double)i;
		}
		while(a<=n/2)
		{
			result+=sum_of_fractions(a,b-1);
			a=a*2;
			b=b*2;
		}
		result+=sum_of_fractions(a,n);
	}
	printf("Approximate Harmonic(%llu) is %.10f\n",n,result);
	return 0;
}
