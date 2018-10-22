// Jason Jue

#include <iostream>
#include "rv.h"
#include "event.h"

// Simulates an M/M/1 queueing system.  The simulation terminates
// once 100000 customers depart from the system.
int main()
{
  using namespace std;
  EventList Elist;                // Create event list
  enum {ARR,DEP};                 // Define the event types

  double lambda = 0.8;            // Arrival rate
  double mu = 1.0;                // Service rate

  double clock = 0.0;             // System clock
  int N = 0;                      // Number of customers in system
  int Ndep = 0;                   // Number of departures from system
  double EN = 0.0;                // For calculating E[N]

  int done = 0;                   // End condition satisfied?

  Event* CurrentEvent;

  Elist.insert(exp_rv(lambda),ARR); // Generate first arrival event

  while (!done)
  {
    CurrentEvent = Elist.get();               // Get next Event from list
    double prev = clock;                      // Store old clock value
    clock=CurrentEvent->time;                 // Update system clock 

    switch (CurrentEvent->type) {
    case ARR:                                 // If arrival 
      EN += N*(clock-prev);                   //  update system statistics
      N++;                                    //  update system size
      Elist.insert(clock+exp_rv(lambda),ARR); //  generate next arrival
      if (N==1) {                             // If this is the only customer
        Elist.insert(clock+exp_rv(mu),DEP);   //  generate its departure event
      }
      break;
    case DEP:                                 // If departure
      EN += N*(clock-prev);                   //  update system statistics
      N--;                                    //  decrement system size
      Ndep++;                                 //  increment num. of departures
      if (N > 0) {                            // If customers remain
        Elist.insert(clock+exp_rv(mu),DEP);   //  generate next departure
      } 
      break;
    }
    delete CurrentEvent;
    if (Ndep > 100000) done=1;        // End condition
  }
  // output simulation results for N, E[N] 
  cout << "Current number of customers in system: " << N << endl;
  cout << "Expected number of customers (simulation): " << EN/clock << endl;

  // output derived value for E[N]
  double rho = lambda/mu; 
  cout << "Expected number of customers (analysis): " << rho/(1-rho) << endl;
}

