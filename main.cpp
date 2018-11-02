// Arunabha Chakraborty(axc175630)

#include <iostream>
#include "rv.h"
#include "event.h"
#include <cmath>
#include <iomanip>

/* Queueing System description:
Consider a Web site that operates m = 2 identical replicated servers to distribute files. Requests to the site arrive
according to a Poisson process with a rate of λ requests per second. However, if the site is already handling K or
more requests (those currently being served or waiting to be served), then, to reduce the load, the site will start to
randomly block subsequent requests with 50% probability. If the number of requests reaches 2K, then all additional
requests will be blocked until the number of requests in the system falls below 2K. There are a total of m = 2 servers;
however, to reduce energy consumption, only one of the servers is operating whenever the number in the system is
less than or equal to K. Whenever the number in the system exceeds K, the second server is placed into operation.
The service time of each request is exponentially distributed with an average service time of 1/µ seconds. Assume that K >= 2.

For the below simulation, k = 4, µ = 3, m = 2, λ = ρ*m*µ, where  ρ = 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 */

int main()
{
  using namespace std;
  int i;
  double rho = 0.0;

  //Loop 10 time for all rho values(0.1,0.2,....,1.0)
  for (i =0; i<10; i++){
	  EventList Elist;                // Create event list
	  enum {ARR,DEP};                 // Define the event types

	  rho = rho + 0.1;
	  int m = 2;                      // Number of servers
	  double mu = 3;                  // Service rate
	  double lambda = rho*m*mu;       // Arrival rate
	  double rh = lambda/mu;

	  double clock = 0.0;             // System clock
	  int N = 0;                      // Number of customers in system
	  int Ndep = 0;                   // Number of departures from system

	  int k = 4;                      // Value of K
	  
	  double EN = 0.0;                // For calculating E[N]

	  int done = 0;                   // End condition satisfied?
	  double blk = 0.0;		  // To count the number of blocked/dropped customers
	  double util = 0.0; 		  // To calculate utilization
	  double random =0.0;

	  cout<<"FOR rho = "<<rho<<"\n=========================================================================================================================\n";
	  Event* CurrentEvent;

	  Elist.insert(exp_rv(lambda),ARR); 	// Generate first arrival event
	  while (!done){
		CurrentEvent = Elist.get();               // Get next Event from list
		if (CurrentEvent == 0){
			//IF the eventlist is empty, generate new arrival event.
			Elist.insert(clock+exp_rv(lambda),ARR);
			//refresh the Current Event
			delete CurrentEvent;
			CurrentEvent = Elist.get();
		}
	    	double prev = clock;			  // Store old clock value
	    	clock=CurrentEvent->time;                 // Update system clock 
		if (N < k){
	    		switch (CurrentEvent->type) {
	    		case ARR:                                 	// If arrival 
	      			EN += N*(clock-prev);				//  update system statistics
	      			N++;                                    	//  update system size
	      			Elist.insert(clock+exp_rv(lambda),ARR); 	//  generate next arrival
	      			if (N==1) {                              	//  If this is the only customer
					Elist.insert(clock+exp_rv(mu),DEP);   	//  generate its departure event
	     	 		}
	      			break;
	    		case DEP:                                 	// If departure
	      			EN += N*(clock-prev);                   	//  update system statistics
				util += 0.5*(clock-prev);
	      			N--;                                    	//  decrement system size
	      			Ndep++;                                 	//  increment num. of departures
	      			if (N > 0) {                            	//  If customers remain
					Elist.insert(clock+exp_rv(mu),DEP);   	//  generate next departure
	      			} 
	      			break;
	    		}
		}
		else if ((N>=k) && (N<(2*k))){
			switch (CurrentEvent->type){
			case ARR:
				EN += N*(clock-prev);
				N++;
				random = uni_rv();
				if (random > 0.5)
					Elist.insert(clock+exp_rv(lambda),ARR); 	//Since randomly 50% requests are blocked
				else if (random < 0.5) 
					blk++;
				break;
			case DEP:
				EN += N*(clock-prev);
				if (N == k){
					util += 0.5*(clock-prev);
				}
				else {
					util += 1*(clock-prev);
				}
				N--;
				Ndep++;
				if(N > 0){
					if (N <= k){
						Elist.insert(clock+exp_rv(mu),DEP);	//if N <= k , 1 server is used 
					}
					else{
						Elist.insert(clock+exp_rv(2*mu),DEP);	//since two servers are active now double the dep rates
					}
				}
				break;
			}

		}
		else if (N>=(2*k)){
			switch (CurrentEvent->type){
			case ARR:					//Arrivals will be blocked if N>2k
				EN += N*(clock-prev);
				//Block further arrivals.
				blk++;					//Increment the number of customers blocked
				break;
			case DEP:
				EN += N*(clock-prev);
				util += 1*(clock-prev);
				N--;
				Ndep++;
				if(N > 0){
					Elist.insert(clock+exp_rv(2*mu),DEP);	//since two servers are active, double the dep rate
				}
				break;
			}
		}
	    delete CurrentEvent;
	    if (Ndep > 100000) done=1;        // End condition
	  }
	  
	  double p0 = 1/(1+rh+(pow(rh,2))+(pow(rh,3))+(pow(rh,4))+((0.25)*(pow(rh,5)))+((0.0625)*(pow(rh,6)))+
	  		((0.015625)*(pow(rh,7)))+((0.00390625)*(pow(rh,8))));
	  double p1 = rh*p0;
	  double p2 = rh*p1;
	  double p3 = rh*p2;
	  double p4 = rh*p3;
	  double p5 = 0.25*(rh*p4);
	  double p6 = 0.25*(rh*p5);
	  double p7 = 0.25*(rh*p6);
	  double p8 = 0.25*(rh*p7);
	  double lambda_half = lambda/2;
	  
	  //Theoretical: EN[X] and average time
	  double EN_analysis = p1+(2*p2)+(3*p3)+(4*p4)+(5*p5)+(6*p6)+(7*p7)+(8*p8);
	  double lambda_avg = (lambda*p0)+(lambda*p1)+(lambda*p2)+(lambda*p3)+(lambda_half*p4)+(lambda_half*p5)+
	  		      (lambda_half*p6)+(lambda_half*p7);
	  double avg_time = EN_analysis/lambda_avg;	//Little's law
	  
	  //Theoretical : blocking probability
	  double numerator = (lambda*0.5*p4)+(lambda*0.5*p5)+(lambda*0.5*p6)+(lambda*0.5*p7)+(lambda*1*p8);
	  double denominator = (lambda*p0)+(lambda*p1)+(lambda*p2)+(lambda*p3)+(lambda_half*p4)+(lambda_half*p5)+
	  		       (lambda_half*p6)+(lambda_half*p7)+(lambda_half*p8*0);
	  double Pblock = numerator/denominator;
	  //double Pblock = (p5/2)+(p6/2)+(p7/2)+(p8);
	  
	  //Theoretical : Utilization
	  double utilization = (0*p0)+(0.5*p1)+(0.5*p2)+(0.5*p3)+(0.5*p4)+(1*p5)+(1*p6)+(1*p7)+(1*p8);

	  cout << setw(50)<<"Expected number of customers (simulation): " << EN/clock <<setw(50)<<"Expected number of customers (analysis): "<< EN_analysis<< endl;
	  cout << setw(50)<<"Average time (simulation): " << EN/(100000) <<setw(50)<<"Average time(analysis): "<<avg_time<<endl;
	  cout << setw(50)<<"Blocking probability(simulation): "<<blk/(100000-blk) <<setw(50)<<"Blocking probability(analysis): "<<Pblock<<endl;
	  cout << setw(50)<<"Utilization(simulation): "<< util/clock<<setw(50)<<"Utilization(analysis): "<<utilization<<endl;

  }
}

