#include<iostream>
#include<random>
#include<math.h>
#include<chrono> 
#include<algorithm>
#include<fstream>

using namespace std;

static random_device rd;  // static seed
static mt19937 gen(rd()); // static generator

bool isBetween(double x, double a, double b) {
    return (x >= a && x <= b);
}

// Function to generate a random real number between min and max
double randomReal(double min, double max) 
{
    uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

double randomAround(double value, double disp) 
{
    uniform_real_distribution<> dis(value*(1.0-disp), value*(1.0+disp));
    return dis(gen);
}

//Notice how I change the sign of the mass to avoid complex numbers

double MinH (double m12, double lh, double lp, double ls, double l1)
{
	return sqrt(-m12)*sqrt(lp) / sqrt(( 8.0*l1*lh-lp*lp+4.0*lh*ls ));
	
}

double MinSR (double m12, double lh, double lp, double ls, double l1)
{
	return 2.0*sqrt(-m12)*sqrt(lh) / sqrt(( 8.0*l1*lh-lp*lp+4.0*lh*ls));
	
}

// a, b, c are the entries of the 2x2 matrix
// a c
// c b

double a (double m12, double lh, double lp, double ls, double l1)
{
	return (8.0*m12*lh*lp ) / (-8.0*l1*lh + lp*lp - 4.0*lh*ls) ;
}

double b (double m12, double lh, double lp, double ls, double l1)
{
	return -(8.0*m12*lh*(2.0*l1+ls)) / (8.0*l1*lh - lp*lp + 4.0*lh*ls) ;
}

double c (double m12, double lh, double lp, double ls, double l1)
{
	return (4.0*m12*sqrt(lh)*pow(lp,1.5)) / (8.0*l1*lh - lp*lp + 4.0*lh*ls) ;
}


//Tree mass eigenvalues

double mass1 (double a, double b, double c)
{
	return 0.5*(a+b-sqrt( pow(a-b,2) +4.0*c*c));
}

double mass2 (double a, double b, double c)
{
	return 0.5*(a+b+sqrt( pow(a-b,2) +4.0*c*c));
}

double mass3 (double m12, double lh, double lp, double ls, double l1)
{
	return 2.0*m12*(8.0*l1*lh + lp*lp - 4.0*lh*ls) / (8.0*l1*lh - lp*lp + 4.0*lh*ls) ;
}

//mixing angle

double theta (double a, double b, double c)
{
	return 0.5*atan(2.0*c/(a-b));
}

vector<double> calc_physical_parameters (double m12, double lh, double lp, double ls, double l1)
{

	
	double Na, Nb, Nc;
	double NmH, NmS1, NmS2;
	double NmHsq, NmS1sq, NmS2sq;
	
	double Ntheta;
	double PhysVev;
	
	//Compute numerical values for the entries of the 2x2 submatrix
	Na = a(m12, lh, lp, ls, l1);
	Nb = b(m12, lh, lp, ls, l1);
	Nc = c(m12, lh, lp, ls, l1);
	
	//Compute mass eigenvalues and angle
	
	NmHsq = mass1(Na, Nb, Nc);
	NmS1sq = mass2(Na, Nb, Nc);
	NmS2sq = mass3(m12, lh, lp, ls, l1);
	
	if(NmHsq > 0 && NmS1sq > 0  && NmS2sq > 0)
	{
		NmH = sqrt(NmHsq);
		NmS1 = sqrt(NmS1sq);
		NmS2 = sqrt(NmS2sq);
		
		//Compute mixing angle and physical Higgs vev 
		//The vev is given by the ExpVal * sqrt(2) 
	
		Ntheta = theta(Na, Nb, Nc);
		PhysVev = MinH(m12, lh, lp, ls, l1)*sqrt(2);
		
		return {NmH, NmS1, NmS2, PhysVev, Ntheta};
	}
	
	else return {0,0,0,0,0};
	
}

int main(){

	ofstream output ("valid_points.txt");
	long n_points = 0;
	double m12, lh, lp, ls, l1;
	double best_m12, best_lh, best_lp, best_ls, best_l1;
	
	vector<double> phys_par = {0,0,0,0,0};
	
	//Define the experimental value of Higgs mass and vev;
	
	double TrueMass = 0.12511;
	double TrueVev = 0.24622;
	
	//A bunch of constants I need in the code
	
	double MaxError = 0.2; 
	double GoalError = 0.001; //The percentual error I am willing to toleare around the true value 
	double CurrentError=1.0;
	double BestError=1.0;
	double dispersion = 0.1;
	long iterations = 0;
	long max_iterations = 50000;
	
	// Get start time
   	auto start = chrono::high_resolution_clock::now();
   	
   	// I look for parameters between 0.05 and 0.5 for naturalness reasons,
   	// ideally we should look for something as close to 1 as possible
   	double range_min = 0.05;
   	double range_max = 5.0;
   	
   		
    	while(n_points < 5000) {
    	
    	
    		// Initialize random initial values
    		
		m12 = randomReal(-1, -0.1);
		lh = randomReal(range_min, range_max);
		lp = randomReal(range_min, range_max);
		ls = randomReal(range_min, range_max);
		do { l1 = randomReal(-range_max, range_max);} while (!isBetween(l1, -range_max,-range_min) && !isBetween(l1,range_min,range_max) );
		
		phys_par = calc_physical_parameters(m12, lh, lp, ls, l1);
		
	
		// Check if the physical masses are all positive
		// and if the numerical value are within the error
		
		// If I find a good value (with error within 10%), I improve
		// it till I find one with error below the Goal error
		
		if( isBetween(phys_par[0], TrueMass*(1.0-MaxError), TrueMass*(1.0+MaxError)) && 
		isBetween(phys_par[3], TrueVev*(1.0-MaxError), TrueVev*(1.0+MaxError)) )
		{
			
			CurrentError = max(abs(1.0-phys_par[0]/TrueMass), abs(1.0-phys_par[3]/TrueVev));
			BestError = CurrentError;
			iterations = 0;
			best_m12 = m12;
			best_l1 = l1;
			best_lh = lh;
			best_lp = lp;
			best_ls = ls;
			while (CurrentError > GoalError && iterations < max_iterations)
			{
				//Check that th parameters have the correct sign
				do
				{
					m12 = randomAround(best_m12, dispersion);
					lh = randomAround(best_lh, dispersion);
					lp = randomAround(best_lp, dispersion);
					ls = randomAround(best_ls, dispersion);
					l1 = randomAround(best_l1, dispersion);
					iterations++;
				
				} while ( !isBetween(m12,-1.0,-0.1) || !isBetween(lh, range_min,range_max) || !isBetween(lp, range_min,range_max) || !isBetween(ls,range_min,range_max) || (!isBetween(l1,-range_max, -range_min) && !isBetween(l1,range_min,range_max) )  );
				
				iterations++;
				
		
				
				phys_par = calc_physical_parameters(m12, lh, lp, ls, l1);
				CurrentError = max(abs(1.0-phys_par[0]/TrueMass), abs(1.0-phys_par[3]/TrueVev));
				
				if(CurrentError < BestError )
				{
					BestError = CurrentError;
					best_m12 = m12;
					best_l1 = l1;
					best_lh = lh;
					best_lp = lp;
					best_ls = ls;
				}
				
			}
			
			if(CurrentError < GoalError && CurrentError == BestError) //Just to be sure I am taking the right par
			{
				output<<best_m12<<" "<<best_lh<<" "<<best_lp<<" "<<best_ls<<" "<<best_l1<<" ";
				for(auto c : phys_par )
				{
					//cout<<c<<" ";
					output<<c<<" ";
				}
				//cout<<endl;
				output<<endl;
				n_points++;
				output.flush();
				//just to check how many points I have already found
				if(n_points%50==0) cout<<"\n"<<n_points;
			}
			
		}		
		
    	}


	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> duration = end - start;

	cout << "Time taken: " << duration.count() << " ms\n";

}



/*
m12 = -0.656599;
   	lh = 0.0568641;
   	lp = 0.0342948;
   	ls = 3.57486;
   	l1 = 0.0384954;
   	
   	//Compute numerical values for the entries of the 2x2 submatrix
	Na = a(m12, lh, lp, ls, l1);
	Nb = b(m12, lh, lp, ls, l1);
	Nc = c(m12, lh, lp, ls, l1);
	
	//Compute mass eigenvalues and angle
	
	NmHsq = mass1(Na, Nb, Nc);
	NmS1sq = mass2(Na, Nb, Nc);
	NmS2sq = mass3(m12, lh, lp, ls, l1);
	
	if(NmHsq > 0 && NmS1sq > 0  && NmS2sq > 0)
	{
		NmH = sqrt(NmHsq);
		NmS1 = sqrt(NmS1sq);
		NmS2 = sqrt(NmS2sq);
	}
	
	
	Ntheta = theta(Na, Nb, Nc);
	PhysVev = MinH(m12, lh, lp, ls, l1)*cos(Ntheta)*sqrt(2);
	
	cout<<"\n"<< NmH << "  " << PhysVev;
		
	// Check if the physical masses are all positive
	// and if the numerical value are within the error

	if( isBetween(NmH, TrueMass*(1.0-MaxError), TrueMass*(1+MaxError)) && 
	isBetween(PhysVev, TrueVev*(1-MaxError), TrueVev*(1+MaxError)) )
	{
		cout<<"\n"<< NmH << "  " << PhysVev;
	}
*/

