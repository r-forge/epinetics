
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "networkapi.h"
#include "utils.h"
#include "epinetics_common.h"

#define MAXINDEX 100
#define MAXIMMUN 1000

#define EVNTMUT (uint16_t) 0
#define EVNTINF (uint16_t) 1
#define EVNTREC (uint16_t) 2
#define WIMMSTIM (uint16_t) 3
#define WEVNTINF (uint16_t) 4
#define WEVNTREC (uint16_t) 5

#define PRINTINFO 0
#define PRINTINFO2 0

/*Structures */

struct event
{
  double rate;
  int function_id;
  int ego_id;
  int altar_id;
  int strain_id;
  int phylo_id;
  struct event *next;
  struct event *prev;
};

/*Global variables*/

struct event *start = NULL;	// Pointer to first element in doubly-linked list
struct event *last = NULL;	// Pointer to last element in doubly-linked list
double sim_T;			// Baseline transmissability, i.e. Pr(transmission over entire infection)
double sim_mu, sim_gamma, sim_beta;	//Mutation rate, recovery rate, instantaneous transmission rate
double sim_time = 0;		// Simulation time
double a = 0.1;			//instantaneous rate of decay of sim_T with distance between new and old strains
double total_rates = 0;		//sum of rate of all events possible at simulation time TIME 
SEXP g;				//pointer to network object 
int case_load = -1;		//counts the number of current infected hosts
double adaptive_immunity_levels[MAXIMMUN]= {0};     //counts of specific immune responses to strains with different phylo_ids
int next_phylo_id = 1;		//next strain phylogenetic id number, as opposed to the normal id that's based on distance
double efficiency = 0.1;        //Instantaneous rate by which each immune predator eats infected prey
double conversion_factor = 0.1; //Each infectious prey consumed by the immune pred is converted into this many new immune preds

struct value *v, *found;
struct hashtable *h;
struct hashtable_itr *itr;
/*Function prototypes */

SEXP epiSimSSA_R (SEXP network, SEXP trans, SEXP infectious_period,
		  SEXP mutation_rate, SEXP max_steps, SEXP verbose,
		  SEXP immune_decay, SEXP make_time_series,
		  SEXP initial_infectious_id);

SEXP within_epiSimSSA_R (SEXP network, SEXP trans, SEXP infectious_period,
		  SEXP mutation_rate, SEXP max_steps, SEXP verbose,
		  SEXP immune_decay, SEXP make_time_series,
		  SEXP initial_infectious_id, SEXP eff, SEXP conv_fac, SEXP innoc_size);
//the functions that gets called from R

int infect (struct key *q, struct value *p);

int recover (struct key *q, struct value *p);

int mutate (struct key *q, struct value *p);

void do_next_event (int make_time_series, double *distance, double *fitness);

int winfect (struct key *q, struct value *p);

int wrecover (struct key *q, struct value *p);

int wmutate (struct key *q, struct value *p);

int wimmstim (struct key *q, struct value *p);

void within_do_next_event (int make_time_series, double *distance, double *fitness);

static unsigned int hashfromkey(void *ky);

static int equalkeys(void *k1, void *k2);


SEXP
epiSimSSA_R (SEXP network, SEXP trans, SEXP infectious_period,
	     SEXP mutation_rate, SEXP max_steps, SEXP verbose,
	     SEXP immune_decay, SEXP make_time_series,
	     SEXP initial_infectious_id)
     /*
      * This function is the main loop for a Gillespie's direct 
      * method stochastic simulation of a SIR model with mutation.
      */
{
  struct key *k;
  int num_nodes, step, max, pc = 0;
  double tolerance = 0.001, rand;
  double *distance, *fitness;

  /* reset global variables (mainly time) for secondary runs */
  struct event *start = NULL;	// Pointer to first element in doubly-linked list
  struct event *last = NULL;	// Pointer to last element in doubly-linked list
  sim_time = 0;			// Simulation time
  a = 0;			//instantaneous rate of decay of sim_T with distance between new and old strains
  total_rates = 0;		//sum of rate of all events possible at simulation time TIME 
  case_load = 0;
  next_phylo_id = 1;

  h = create_hashtable (16, hashfromkey, equalkeys);
  if (NULL == h)
    {
      error ("%s:%d: Unable to create hash table\n", __FILE__, __LINE__);
    }

  /*Verify that we were called properly */
  netRegisterFunctions ();
  g = network;
  if (!netIsNetwork (g))
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with a network object\n", __FILE__, __LINE__);
    }
  if (netIsDir (g))
    {
      error ("%s:%d: epiSimSSA_R must be \
 called with an undirected network\n", __FILE__, __LINE__);
    }
  PROTECT (trans = coerceVector (trans, REALSXP));
  pc++;
  PROTECT (infectious_period = coerceVector (infectious_period, REALSXP));
  pc++;
  PROTECT (mutation_rate = coerceVector (mutation_rate, REALSXP));
  pc++;
  PROTECT (max_steps = coerceVector (max_steps, INTSXP));
  pc++;
  PROTECT (immune_decay = coerceVector (immune_decay, REALSXP));
  pc++;

  /*Store simulation parameters as network attributes for later reference */
  sim_T = REAL (trans)[0];
  if (sim_T > 0.99)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with trans in [0,0.99]\n", __FILE__, __LINE__);
    }
  sim_mu = REAL (mutation_rate)[0];
  if (sim_mu < 0 || sim_mu > 1000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with mutation rate in [0, 1000]\n", __FILE__, __LINE__);
    }
  sim_gamma = 1 / REAL (infectious_period)[0];
  if (sim_gamma <= 0.001 || sim_gamma > 1000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with a mean infectious period in  [0.001, 1000]\n", __FILE__, __LINE__);
    }
  a = REAL (immune_decay)[0];
  if (a < 0 || a > 1000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with immune_decay in  [0, 1000]\n", __FILE__, __LINE__);
    }
  max = INTEGER (max_steps)[0];

  if (max < 0 || max > 1000000000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with max.steps in  [0, 1e9]\n", __FILE__, __LINE__);
    }
  sim_beta = -sim_T / ((1 / sim_gamma) * (sim_T - 1));

  g = netSetNetAttrib (g, "transmissability", trans);
  g = netSetNetAttrib (g, "mutation_rate", mutation_rate);
  g = netSetNetAttrib (g, "infectious_period", infectious_period);
  g = netSetNetAttrib (g, "immune_decay", immune_decay);

  /*Setup random number generator */
  GetRNGstate ();
  
  k = (struct key *) malloc (sizeof (struct key));
  if (NULL == k)
    {
      hashtable_destroy (h, 1);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
    }

  v = (struct value *) malloc (sizeof (struct value));
  if (NULL == v)
    {
      hashtable_destroy (h, 1);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
      free (k);
    }
  k->ego_id = 0;
  k->altar_id = INTEGER (initial_infectious_id)[0];
  v->strain_id = 0;
  v->phylo_id = next_phylo_id;
  next_phylo_id++;
  infect (k, v);

  step = 0;

  if (INTEGER (verbose)[0])
    {
      while (step <= max && total_rates >= tolerance)
	{
	  step++;
	  Rprintf ("steps = %2d, total_rates=%f, sim_time=%f\n", step,
		   total_rates, sim_time);
          Rprintf("\t %d items in table\n", hashtable_count(h));
	  do_next_event (INTEGER (make_time_series)[0], distance, fitness);
	}
    }
  else
    {
      while (step <= max && total_rates >= tolerance)
	{
	  step++;
	  do_next_event (INTEGER (make_time_series)[0], distance, fitness);
	}
    }
  if (step > max)
    {
      Rprintf ("Reached maximum number of steps in max.steps argument.\n");
    }
  if (hashtable_count(h))
    {
      hashtable_destroy (h, 1);
      error("Items left in table at end of simulation\n");
    }
  /*Clear protection stack and return */
  UNPROTECT (pc);
  return g;
}

SEXP
seq_epiSimSSA_R (SEXP network, SEXP trans, SEXP infectious_period,
	     SEXP mutation_rate, SEXP max_steps, SEXP verbose,
	     SEXP immune_decay, SEXP make_time_series,
	     SEXP initial_infectious_id, SEXP prop_mutant)
     /*
      * This function is the main loop for a Gillespie's direct 
      * method stochastic simulation of a SIR model with mutation.
      */
{
  struct key *k;
  int num_nodes, step, max, pc = 0, init_id, i;
  double tolerance = 0.001, rand, prop;
  double *distance, *fitness;

  /* reset global variables (mainly time) for secondary runs */
  struct event *start = NULL;	// Pointer to first element in doubly-linked list
  struct event *last = NULL;	// Pointer to last element in doubly-linked list
  sim_time = 0;			// Simulation time
  a = 0;			//instantaneous rate of decay of sim_T with distance between new and old strains
  total_rates = 0;		//sum of rate of all events possible at simulation time TIME 
  case_load = 0;
  next_phylo_id = 1;

  h = create_hashtable (16, hashfromkey, equalkeys);
  if (NULL == h)
    {
      error ("%s:%d: Unable to create hash table\n", __FILE__, __LINE__);
    }
  
  /*Verify that we were called properly */
  netRegisterFunctions ();
  g = network;
  if (!netIsNetwork (g))
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with a network object\n", __FILE__, __LINE__);
    }
  if (netIsDir (g))
    {
      error ("%s:%d: epiSimSSA_R must be \
 called with an undirected network\n", __FILE__, __LINE__);
    }
  PROTECT (trans = coerceVector (trans, REALSXP));
  pc++;
  PROTECT (infectious_period = coerceVector (infectious_period, REALSXP));
  pc++;
  PROTECT (mutation_rate = coerceVector (mutation_rate, REALSXP));
  pc++;
  PROTECT (max_steps = coerceVector (max_steps, INTSXP));
  pc++;
  PROTECT (immune_decay = coerceVector (immune_decay, REALSXP));
  pc++;

  /*Store simulation parameters as network attributes for later reference */
  sim_T = REAL (trans)[0];
  if (sim_T > 0.99)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with trans in [0,0.99]\n", __FILE__, __LINE__);
    }
  sim_mu = REAL (mutation_rate)[0];
  if (sim_mu < 0 || sim_mu > 1000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with mutation rate in [0, 1000]\n", __FILE__, __LINE__);
    }
  sim_gamma = 1 / REAL (infectious_period)[0];
  if (sim_gamma <= 0.001 || sim_gamma > 1000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with a mean infectious period in  [0.001, 1000]\n", __FILE__, __LINE__);
    }
  a = REAL (immune_decay)[0];
  if (a < 0 || a > 1000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with immune_decay in  [0, 1000]\n", __FILE__, __LINE__);
    }
  max = INTEGER (max_steps)[0];

  if (max < 0 || max > 1000000000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with max.steps in  [0, 1e9]\n", __FILE__, __LINE__);
    }
  sim_beta = -sim_T / ((1 / sim_gamma) * (sim_T - 1));

  g = netSetNetAttrib (g, "transmissability", trans);
  g = netSetNetAttrib (g, "mutation_rate", mutation_rate);
  g = netSetNetAttrib (g, "infectious_period", infectious_period);
  g = netSetNetAttrib (g, "immune_decay", immune_decay);

  /*Setup random number generator */
  GetRNGstate ();
  
  k = (struct key *) malloc (sizeof (struct key));
  if (NULL == k)
    {
      hashtable_destroy (h, 1);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
    }

  v = (struct value *) malloc (sizeof (struct value));
  if (NULL == v)
    {
      hashtable_destroy (h, 1);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
      free (k);
    }

  init_id = INTEGER (initial_infectious_id)[0];
  prop = REAL (prop_mutant)[0];
  for (i = init_id; i < prop*(init_id + 10); i++){
      k->ego_id = 0;
      k->altar_id = i;
      v->strain_id = 53;
      v->phylo_id = next_phylo_id;
      next_phylo_id++;
      infect (k, v);
  }
  for (; i < init_id + 10; i++){
      k->ego_id = 0;
      k->altar_id = i;
      v->strain_id = 0;
      v->phylo_id = next_phylo_id;
      next_phylo_id++;
      infect (k, v);
  }



  step = 0;

  if (INTEGER (verbose)[0])
    {
      while (step <= max && total_rates >= tolerance)
	{
	  step++;
	  Rprintf ("steps = %2d, total_rates=%f, sim_time=%f\n", step,
		   total_rates, sim_time);
          Rprintf("\t %d items in table\n", hashtable_count(h));
	  do_next_event (INTEGER (make_time_series)[0], distance, fitness);
	}
    }
  else
    {
      while (step <= max && total_rates >= tolerance)
	{
	  step++;
	  do_next_event (INTEGER (make_time_series)[0], distance, fitness);
	}
    }
  if (step > max)
    {
      Rprintf ("Reached maximum number of steps in max.steps argument.\n");
    }
  if (hashtable_count(h))
    {
      hashtable_destroy (h, 1);
      error("Items left in table at end of simulation\n");
    }
  /*Clear protection stack and return */
  UNPROTECT (pc);
  return g;
}


int
infect (struct key *q, struct value *p)
{
  /* host EGO_ID infects host ALTAR_ID with strain STRAIN_ID at time SIM_TIME */

  struct key *k, *kk;
  struct value *found;
  SEXP sim_time_p, strain_p, status_p, vall, val, infected_status_p;
  SEXP hood, infector_id_p, sus_status_p, altar_nb_strain_id_p;
  SEXP phylo_id_p;
  int pc = 0, i, altar_nb_id, distance, altar_id, strain_id, phylo_id;
  double tau, rate;

  /*Assignments to help avoid excessive pointer confusion
   * and simply some lines of code */
  strain_id = p->strain_id;
  phylo_id = p->phylo_id;
  altar_id = q->altar_id;

  /*Add vertex attributes */
  PROTECT (sim_time_p = allocVector (REALSXP, 1));
  pc++;
  PROTECT (strain_p = allocVector (INTSXP, 1));
  pc++;
  PROTECT (infector_id_p = allocVector (INTSXP, 1));
  pc++;
  PROTECT (phylo_id_p = allocVector (INTSXP, 1));
  pc++;
  PROTECT (infected_status_p = allocVector (STRSXP, 1));
  pc++;
  PROTECT (sus_status_p = allocVector (STRSXP, 1));
  pc++;
  REAL (sim_time_p)[0] = sim_time;
  INTEGER (infector_id_p)[0] = q->ego_id;
  INTEGER (strain_p)[0] = strain_id;
  INTEGER (phylo_id_p)[0] = phylo_id;
  SET_STRING_ELT (infected_status_p, 0, mkChar ("infectious"));
  SET_STRING_ELT (sus_status_p, 0, mkChar ("susceptible"));
  g = netSetVertexAttrib (g, "time_infected", sim_time_p, altar_id);
  g = netSetVertexAttrib (g, "infection_history", strain_p, altar_id);
  g = netSetVertexAttrib (g, "status", infected_status_p, altar_id);
  g = netSetVertexAttrib (g, "infector_id", infector_id_p, altar_id);
  g = netSetVertexAttrib (g, "phylo_id", phylo_id_p, altar_id);

  /*update case load counter */
  case_load++;
  /*Add event of strain mutating to the table */
  k = (struct key *) malloc (sizeof (struct key));
  if (NULL == k)
    {
      hashtable_destroy (h, 1);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
    }
  k->ego_id = altar_id;
  k->altar_id = 0;
  k->function_id = EVNTMUT;
  k->two_port = 0;

  v = (struct value *) malloc (sizeof (struct value));
  if (NULL == v)
    {
      hashtable_destroy (h, 1);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
      free(k);
    }
  v->rate = sim_mu;
  total_rates += v->rate;
  if (!insert_some (h, k, v))
    {
      hashtable_destroy (h, 1);
      free (k);
      error ("%s:%d: Insertion into hash table failed\n", __FILE__, __LINE__);
    }

  /*Add event of host recovering to the table */
  k = (struct key *) malloc (sizeof (struct key));
  if (NULL == k)
    {
      hashtable_destroy (h, 1);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
    }
  k->ego_id = altar_id;
  k->altar_id = 0;
  k->function_id = EVNTREC;
  k->two_port = 0;

  v = (struct value *) malloc (sizeof (struct value));
  if (NULL == v)
    {
      hashtable_destroy (h, 1);
      free (k);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
    }
  v->rate = sim_gamma;
  total_rates += v->rate;

  if (!insert_some (h, k, v))
    {
      hashtable_destroy (h, 1);
      free (k);
      error ("%s:%d: Insertion into hash table failed\n", __FILE__, __LINE__);
    }

  /* Add event of infected host infecting suceptible neighbors */
  PROTECT (hood = netGetNeighborhood (g, altar_id, "out", 1));
  pc++;
  PROTECT (val = getListElement (g, "val"));
  pc++;
  for (i = 0; i < length (hood); i++)
    {
      altar_nb_id = INTEGER (hood)[i];
      PROTECT (vall =
	       getListElement (VECTOR_ELT (val, altar_nb_id - 1), "status"));
      pc++;
     if (vecEq (vall, sus_status_p))
	{
	  k = (struct key *) malloc (sizeof (struct key));
	  if (NULL == k)
	    {
	      hashtable_destroy (h, 1);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
	  k->ego_id = altar_id;
	  k->altar_id = altar_nb_id;
	  k->function_id = EVNTINF;
	  k->two_port = 0;

	  v = (struct value *) malloc (sizeof (struct value));
	  if (NULL == v)
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
		v->strain_id = strain_id;
		v->phylo_id = phylo_id;
/*This does strange things for some reason, so assign RHS to ints above
 * v->strain_id = p->strain_id;
 *	  	  v->phylo_id = p->phylo_id;
 */
	  PROTECT (vall =
		   getListElement (VECTOR_ELT (val, altar_nb_id - 1),
				   "immune_memory"));
	  pc++;
	  distance = p->strain_id + fabs (INTEGER (vall)[0]);
	  tau = exp (-a * distance);
	  rate = sim_beta * (1 - tau);
	  v->rate = rate;
	  total_rates += v->rate;

	  if (!insert_some (h, k, v))
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s:%d: Insertion into hash table failed\n", __FILE__,
		     __LINE__);
	    }
	}
    }
  /* Remove event of infected host geting infected by infectious neighbors */
  /* For some reason, doing this before adding the infection events switches
   * the ego_id in q, so it is done last*/
  for (i = 0; i < length (hood); i++)
    {
      altar_nb_id = INTEGER (hood)[i];
      PROTECT (vall =
	       getListElement (VECTOR_ELT (val, altar_nb_id - 1), "status"));
      pc++;
      k = (struct key *) malloc (sizeof (struct key));
      if (NULL == k)
        {
          hashtable_destroy (h, 1);
          error ("%s: %d: Ran out of memory allocating a key\n",
                 __FILE__, __LINE__);
        }
      if (vecEq (vall, infected_status_p))
	{
	  k->function_id = EVNTINF;
	  k->ego_id = altar_nb_id;
	  k->altar_id = altar_id;
          k->two_port = 0;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
	      free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free(found);
	}
    }
  UNPROTECT (pc);
  return (0);
}

int
mutate (struct key *q, struct value *p)
{
  /* host EGO_ID has his infection history incremented  at time SIM_TIME */


  struct key *k;
  SEXP sim_time_p, vall, val;
  SEXP hood, sus_status_p;
  SEXP strain_id_p;
  SEXP new_phylo_id_p;
  SEXP old_mutation_times_vector_p;
  SEXP new_mutation_times_vector_p;
  SEXP old_phylo_id_vector_p;
  SEXP new_phylo_id_vector_p;
  int pc = 0, i, ego_id, ego_nb_id, Rc_int = 0, new_strain_id;
  int new_phylo_id;
  double distance, tau, rate;

  /*Assignments to help avoid pointer confusion*/
  ego_id = q->ego_id;

  /*Add and update vertex attributes */
  PROTECT (sim_time_p = allocVector (REALSXP, 1));
  pc++;
  PROTECT (sus_status_p = allocVector (STRSXP, 1));
  pc++;
  REAL (sim_time_p)[0] = sim_time;
  SET_STRING_ELT (sus_status_p, 0, mkChar ("susceptible"));


  PROTECT (val = getListElement (g, "val"));
  pc++;				//get list of vertex attributes

  PROTECT (old_mutation_times_vector_p =
	   getListElement (VECTOR_ELT (val, q->ego_id - 1), "time_mutated"));
  pc++;
  PROTECT (new_mutation_times_vector_p =
	   vecAppend (old_mutation_times_vector_p, sim_time_p));
  pc++;

  g =
    netSetVertexAttrib (g, "time_mutated", new_mutation_times_vector_p,
			q->ego_id);

  PROTECT (old_phylo_id_vector_p =
	   getListElement (VECTOR_ELT (val, q->ego_id - 1), "phylo_id"));
  pc++;
  new_phylo_id = next_phylo_id;
  next_phylo_id++;
  PROTECT (new_phylo_id_p = allocVector (INTSXP, 1));
  pc++;
  INTEGER (new_phylo_id_p)[0] = new_phylo_id;
  PROTECT (new_phylo_id_vector_p =
	   vecAppend (old_phylo_id_vector_p, new_phylo_id_p));
  pc++;

  g = netSetVertexAttrib (g, "phylo_id", new_phylo_id_vector_p, q->ego_id);

  PROTECT (strain_id_p =
	   getListElement (VECTOR_ELT (val, q->ego_id - 1),
			   "infection_history"));
  pc++;
  new_strain_id = INTEGER (strain_id_p)[0] + 1;
  INTEGER (strain_id_p)[0] = new_strain_id;
  g = netSetVertexAttrib (g, "infection_history", strain_id_p, q->ego_id);


  /* Update the rate and strain_id elements for the events of this newly
   * mutated strain spreading from EGO_ID */

  PROTECT (hood = netGetNeighborhood (g, q->ego_id, "out", 1));
  pc++;
  for (i = 0; i < length (hood); i++)
    {
      ego_nb_id = INTEGER (hood)[i];
      PROTECT (vall =
	       getListElement (VECTOR_ELT (val, ego_nb_id - 1), "status"));
      pc++;

      if (vecEq (vall, sus_status_p))
	{
	  k = (struct key *) malloc (sizeof (struct key));
	  if (NULL == k)
	    {
	      hashtable_destroy (h, 1);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
	  k->ego_id = ego_id;
	  k->altar_id = ego_nb_id;
	  k->function_id = EVNTINF;
	  k->two_port = 0;

        if (NULL == (v = search_some(h,k))) {
	      error ("%s:%d: Search for key in hash table failed\n", __FILE__,
		     __LINE__);
        }
	  v->strain_id = new_strain_id;
	  v->phylo_id = new_phylo_id;

          /* update the rate */
	  PROTECT (vall =
		   getListElement (VECTOR_ELT (val, ego_nb_id - 1),
				   "immune_memory"));
	  pc++;
	  distance = new_strain_id + fabs (INTEGER (vall)[0]);
	  tau = exp (-a * distance);
	  rate = sim_beta * (1 - tau);
          total_rates -= v->rate;
	  v->rate = rate;
	  total_rates += rate;

/*	  if (!hashtable_change (h, k, v))
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s:%d: Modification of key in hash table failed\n", __FILE__,
		     __LINE__);
	    }*/
       }
    }
  UNPROTECT (pc);
  return (0);
}


int
recover (struct key *q, struct value *p)
{
  /* host EGO_ID recovers at time SIM_TIME */

  SEXP sim_time_p, status_p, vall, val, infected_status_p;
  SEXP hood, infector_id_p, sus_status_p, rec_status_p;
  SEXP strain_id_p;
  int pc = 0, i, ego_nb_id, ego_id;
  struct key *k;

  ego_id = q->ego_id;

  /*Add and update vertex attributes */
  PROTECT (sim_time_p = allocVector (REALSXP, 1));
  pc++;
  PROTECT (infected_status_p = allocVector (STRSXP, 1));
  pc++;
  PROTECT (sus_status_p = allocVector (STRSXP, 1));
  pc++;
  PROTECT (rec_status_p = allocVector (STRSXP, 1));
  pc++;
  REAL (sim_time_p)[0] = sim_time;
  SET_STRING_ELT (infected_status_p, 0, mkChar ("infectious"));
  SET_STRING_ELT (sus_status_p, 0, mkChar ("susceptible"));
  SET_STRING_ELT (rec_status_p, 0, mkChar ("recovered"));
  g = netSetVertexAttrib (g, "time_recovered", sim_time_p, q->ego_id);
  g = netSetVertexAttrib (g, "status", rec_status_p, q->ego_id);

  /*update case load counter */
  case_load--;

  /*Remove event of the newly recovered host again recovering from the dl list */
	  k = (struct key *) malloc (sizeof (struct key));
	  if (NULL == k)
	    {
	      hashtable_destroy (h, 1);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
	  k->ego_id = ego_id;
	  k->altar_id = 0;
	  k->function_id = EVNTREC;
          k->two_port = 0;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
	      free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free (found);

          /*Remove event of strain mutating to the events array and doubly linked list */
	  k = (struct key *) malloc (sizeof (struct key));
          if (NULL == k)
            {
              hashtable_destroy (h, 1);
              error ("%s: %d: Ran out of memory allocating a key\n",
                     __FILE__, __LINE__);
            }
	  k->ego_id = ego_id;
	  k->altar_id = 0;
	  k->function_id = EVNTMUT;
          k->two_port = 0;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }

          total_rates -= found->rate;
          free (found);

          /* Remove event of newly recovered host infecting suceptible neighbors from pointer array and dl list */
  PROTECT (hood = netGetNeighborhood (g, ego_id, "out", 1));
  pc++;
  PROTECT (val = getListElement (g, "val"));
  pc++;				//get list of vertex attributes
  PROTECT (strain_id_p =
	   getListElement (VECTOR_ELT (val, ego_id - 1),
			   "infection_history"));
  pc++;
  for (i = 0; i < length (hood); i++)
    {
      ego_nb_id = INTEGER (hood)[i];
      PROTECT (vall =
	       getListElement (VECTOR_ELT (val, ego_nb_id - 1), "status"));
      pc++;
      if (vecEq (vall, sus_status_p))
	{
	  k = (struct key *) malloc (sizeof (struct key));
          if (NULL == k)
            {
              hashtable_destroy (h, 1);
              error ("%s: %d: Ran out of memory allocating a key\n",
                     __FILE__, __LINE__);
            }
	  k->ego_id = ego_id;
	  k->altar_id = ego_nb_id;
	  k->function_id = EVNTINF;
          k->two_port = 0;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free (found);
	}

    }
  UNPROTECT (pc);
  return (0);
}

void
simulator_free (void)
{
/* the if (h) does not avoid double freeing hash tables, so I'll just turn this 
 * off until I figure out the right way to do it
  	if (h)
    {
      Rprintf("At termination, %d items in table\n", hashtable_count(h));
      hashtable_destroy (h, 1);
    }
    */
}

void
do_next_event (int make_time_series, double *distance, double *fitness)
{
  int infector_id, number_infectious, i, pc = 0;
  double rand_one, rand_two, rate_sum, time_interval;
  double fitness_sum, distance_sum, dist_by_fit_sum;
  double dist_fit_cov;
  struct key *k;

  SEXP old_time_series_vector_p;
  SEXP new_time_series_vector_p;
  SEXP old_covariance_series_vector_p;
  SEXP new_covariance_series_vector_p;
  SEXP old_case_series_vector_p;
  SEXP new_case_series_vector_p;
  SEXP old_fitness_sum_series_vector_p;
  SEXP new_fitness_sum_series_vector_p;
  SEXP old_mean_distance_series_vector_p;
  SEXP new_mean_distance_series_vector_p;
  SEXP sim_time_p;
  SEXP cov_p;
  SEXP case_p;
  SEXP fitness_p;
  SEXP distance_p;
  int num_infectious_crit = make_time_series;
  rand_one = runif (0.0, 1.0);
  time_interval = -log (rand_one) / total_rates;
  sim_time += time_interval;
  rand_two = runif (0.0, total_rates);

/*
  rate_sum = 0;
  itr = hashtable_iterator(h);
  do {
        k = hashtable_iterator_key(itr);
        v = hashtable_iterator_value(itr);
        rate_sum += v->rate;
      } while(hashtable_iterator_advance(itr));
  free(itr);
*/


  rate_sum = 0;
    itr = hashtable_iterator(h);
    if (hashtable_count(h) > 0)
    {
        do {
            k = hashtable_iterator_key(itr);
            v = hashtable_iterator_value(itr);
            rate_sum += v->rate;
        } while (rate_sum < rand_two && hashtable_iterator_advance(itr));
    }
  if (!hashtable_iterator_advance(itr) && rate_sum < total_rates - 0.01)
    {
      error ("%s:%d: Reached end of dl list before (rate_sum = %g) >= (total_rates - 0.01 = %g)\n",
	     __FILE__, __LINE__, rate_sum, total_rates - 0.01);
    }
  
  free(itr);
  switch (k->function_id)
    {
    case EVNTINF:
#if PRINTINFO2
      Rprintf ("%d infects %d\n", k->ego_id, k->altar_id);
#endif
      infect (k, v);
      break;
    case EVNTMUT:
#if PRINTINFO2
      Rprintf ("Mutation of strain in %d\n", k->ego_id);
#endif
      mutate (k, v);
      break;
    case EVNTREC:
#if PRINTINFO2
      Rprintf ("Recovery of  %d\n", k->ego_id);
#endif
      recover (k, v);
      break;
    default:
      error ("%s:%d: invalid function_id in key\n", __FILE__,
	     __LINE__);
    }
  /* Put away random number generator. */
  PutRNGstate ();
  UNPROTECT (pc);
}

static unsigned int
hashfromkey(void *ky)
{
    struct key *k = (struct key *)ky;
    return (((k->ego_id << 17) | (k->ego_id >> 15)) ^ k->altar_id) +
            (k->function_id * 17) + (k->two_port * 13 * 29);
}

static int
equalkeys(void *k1, void *k2)
{
    return (0 == memcmp(k1,k2,sizeof(struct key)));
}

/******************************************************************************/

SEXP
within_epiSimSSA_R (SEXP network, SEXP trans, SEXP infectious_period,
	     SEXP mutation_rate, SEXP max_steps, SEXP verbose,
	     SEXP immune_decay, SEXP make_time_series,
	     SEXP initial_infectious_id, SEXP eff, SEXP conv_fac, SEXP innoc_size)
     /*
      * This function is the main loop for a Gillespie's direct 
      * method stochastic simulation of a within-host infection model with mutation.
      */
{
  struct key *k;
  int num_nodes, step, max, pc = 0, i, init_id;
  double tolerance = 0.001, rand;
  double *distance, *fitness;
  SEXP imm_levels_p;

  /* reset global variables for secondary runs */
  struct event *start = NULL;	// Pointer to first element in doubly-linked list
  struct event *last = NULL;	// Pointer to last element in doubly-linked list
  sim_time = 0;			// Simulation time
  a = 0;			//instantaneous rate of decay of sim_T with distance between new and old strains
  total_rates = 0;		//sum of rate of all events possible at simulation time TIME 
  case_load = 0;

  for (i = 0; i < MAXIMMUN; i++)
    {
      adaptive_immunity_levels[i]=1;
    }
  next_phylo_id = 1;

  h = create_hashtable (16, hashfromkey, equalkeys);
  if (NULL == h)
    {
      error ("%s:%d: Unable to create hash table\n", __FILE__, __LINE__);
    }

  /*Verify that we were called properly */
  netRegisterFunctions ();
  g = network;
  if (!netIsNetwork (g))
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with a network object\n", __FILE__, __LINE__);
    }
  if (netIsDir (g))
    {
      error ("%s:%d: epiSimSSA_R must be \
 called with an undirected network\n", __FILE__, __LINE__);
    }
  PROTECT (trans = coerceVector (trans, REALSXP));
  pc++;
  PROTECT (infectious_period = coerceVector (infectious_period, REALSXP));
  pc++;
  PROTECT (mutation_rate = coerceVector (mutation_rate, REALSXP));
  pc++;
  PROTECT (max_steps = coerceVector (max_steps, INTSXP));
  pc++;
  PROTECT (immune_decay = coerceVector (immune_decay, REALSXP));
  pc++;
  PROTECT (eff = coerceVector (eff, REALSXP));
  pc++;
  PROTECT (conv_fac = coerceVector (conv_fac, REALSXP));
  pc++;
  PROTECT (innoc_size = coerceVector (innoc_size, INTSXP));
  pc++;

  /*Store simulation parameters as network attributes for later reference */
  sim_T = REAL (trans)[0];
  if (sim_T > 0.99)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with trans in [0,0.99]\n", __FILE__, __LINE__);
    }
  sim_mu = REAL (mutation_rate)[0];
  if (sim_mu < 0 || sim_mu > 1000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with mutation rate in [0, 1000]\n", __FILE__, __LINE__);
    }
  sim_gamma = 1 / REAL (infectious_period)[0];
  if (sim_gamma <= 0.001 || sim_gamma > 1000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with a mean infectious period in  [0.001, 1000]\n", __FILE__, __LINE__);
    }
  a = REAL (immune_decay)[0];
  if (a < 0 || a > 1000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with immune_decay in  [0, 1000]\n", __FILE__, __LINE__);
    }
  max = INTEGER (max_steps)[0];

  if (max < 0 || max > 1000000000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with max.steps in  [0, 1e9]\n", __FILE__, __LINE__);
    }
  sim_beta = -sim_T / ((1 / sim_gamma) * (sim_T - 1));

  efficiency = REAL(eff)[0];
  if (efficiency < 0 || efficiency > 1000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with efficiency in [0, 1000]\n", __FILE__, __LINE__);
    }
  
  conversion_factor = REAL(conv_fac)[0];
  if (conversion_factor < 0 || conversion_factor > 1000)
    {
      error ("%s:%d: epiSimSSA_R must be\
 called with conversion_factor in [0, 1000]\n", __FILE__, __LINE__);
    }


  g = netSetNetAttrib (g, "transmissability", trans);
  g = netSetNetAttrib (g, "mutation_rate", mutation_rate);
  g = netSetNetAttrib (g, "infectious_period", infectious_period);
  g = netSetNetAttrib (g, "immune_decay", immune_decay);
  g = netSetNetAttrib (g, "efficiency", eff);
  g = netSetNetAttrib (g, "conversion_factor", conv_fac);
 
  /*Setup random number generator */
  GetRNGstate ();
  
  k = (struct key *) malloc (sizeof (struct key));
  if (NULL == k)
    {
      hashtable_destroy (h, 1);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
    }

  v = (struct value *) malloc (sizeof (struct value));
  if (NULL == v)
    {
      hashtable_destroy (h, 1);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
      free (k);
    }

  init_id = INTEGER (initial_infectious_id)[0];
  for (i = init_id; i < init_id + INTEGER(innoc_size)[0]; i++){
      k->ego_id = 0;
      k->altar_id = i;
      v->strain_id = 0;
      v->phylo_id = next_phylo_id;
      winfect (k, v);
  }
  next_phylo_id++;

  step = 0;

  if (INTEGER (verbose)[0])
    {
      while (step <= max && total_rates >= tolerance)
	{
	  step++;
	  Rprintf ("steps = %2d, total_rates=%f, sim_time=%f\n", step,
		   total_rates, sim_time);
          Rprintf("\t %d items in table\n", hashtable_count(h));
	  within_do_next_event (INTEGER (make_time_series)[0], distance, fitness);
	}
    }
  else
    {
      while (step <= max && total_rates >= tolerance)
	{
	  step++;
	  within_do_next_event (INTEGER (make_time_series)[0], distance, fitness);
	}
    }
  if (step > max)
    {
      Rprintf ("Reached maximum number of steps in max.steps argument.\n");
    }
  if (hashtable_count(h))
    {
      itr = hashtable_iterator(h);
      do {
          k = hashtable_iterator_key(itr);
          Rprintf("fun: %d, ego: %d, altar %d, twoport: %d\n", 
                  k->function_id, k->ego_id, k->altar_id, k->two_port);
      } while (hashtable_iterator_advance(itr));

      hashtable_destroy (h, 1);
      error("Items with above keys left in table at end of simulation\n");
    }
  
  PROTECT (imm_levels_p = allocVector (REALSXP, next_phylo_id - 1));
  pc++;
  for (i = 1; i < next_phylo_id; i++)
    {
      REAL (imm_levels_p)[i-1] = adaptive_immunity_levels[i];
    }
  PROTECT (netSetNetAttrib (g, "imm_levels", imm_levels_p));
  pc++;
  /*Clear protection stack and return */
  UNPROTECT (pc);
  return g;
}

void
within_do_next_event (int make_time_series, double *distance, double *fitness)
{
  int infector_id, number_infectious, i, pc = 0;
  double rand_one, rand_two, rate_sum, time_interval;
  double fitness_sum, distance_sum, dist_by_fit_sum;
  double dist_fit_cov;
  struct key *k;
  
  
  int num_infectious_crit = make_time_series;
  rand_one = runif (0.0, 1.0);
  time_interval = -log (rand_one) / total_rates;
  sim_time += time_interval;
  rand_two = runif (0.0, total_rates);

/*
  rate_sum = 0;
  itr = hashtable_iterator(h);
  do {
        k = hashtable_iterator_key(itr);
        v = hashtable_iterator_value(itr);
        rate_sum += v->rate;
      } while(hashtable_iterator_advance(itr));
  free(itr);
*/


  rate_sum = 0;
    itr = hashtable_iterator(h);
    if (hashtable_count(h) > 0)
    {
        do {
            k = hashtable_iterator_key(itr);
            v = hashtable_iterator_value(itr);
            rate_sum += v->rate;
        } while (rate_sum < rand_two && hashtable_iterator_advance(itr));
    }
  if (!hashtable_iterator_advance(itr) && rate_sum < total_rates - 0.01)
    {
      error ("%s:%d: Reached end of dl list before (rate_sum = %g) >= (total_rates - 0.01 = %g)\n",
	     __FILE__, __LINE__, rate_sum, total_rates - 0.01);
    }
  
  free(itr);
  switch (k->function_id)
    {
    case EVNTINF:
#if PRINTINFO2
      Rprintf ("%d infects %d\n", k->ego_id, k->altar_id);
#endif
      winfect (k, v);
      break;
    case WIMMSTIM:
#if PRINTINFO2
      Rprintf ("Stimulation of immunity by strain in %d\n", k->ego_id);
#endif
      wimmstim(k, v);
      break;
    case EVNTREC:
#if PRINTINFO2
      Rprintf ("Recovery of  %d\n", k->ego_id);
#endif
      wrecover (k, v);
      break;
    default:
      error ("%s:%d: invalid function_id in key\n", __FILE__,
	     __LINE__);
    }
  /* Put away random number generator. */
  PutRNGstate ();
  UNPROTECT (pc);
}

int
winfect (struct key *q, struct value *p)
{
  /* host EGO_ID infects host ALTAR_ID with strain STRAIN_ID at time SIM_TIME */

  struct key *k, *kk;
  struct value *found;
  SEXP sim_time_p, strain_p, status_p, vall, val, infected_status_p;
  SEXP hood, infector_id_p, sus_status_p, altar_nb_strain_id_p;
  SEXP phylo_id_p;
  int pc = 0, i, altar_nb_id, distance, altar_id, strain_id, phylo_id, two_port;
  double tau, rate;

  /*Assignments to help avoid excessive pointer confusion
   * and simply some lines of code */
  strain_id = p->strain_id;
  phylo_id = p->phylo_id;
  altar_id = q->altar_id;
  two_port = q->two_port;

  /*Add vertex attributes */
  PROTECT (sim_time_p = allocVector (REALSXP, 1));
  pc++;
  PROTECT (strain_p = allocVector (INTSXP, 1));
  pc++;
  PROTECT (infector_id_p = allocVector (INTSXP, 1));
  pc++;
  PROTECT (phylo_id_p = allocVector (INTSXP, 1));
  pc++;
  PROTECT (infected_status_p = allocVector (STRSXP, 1));
  pc++;
  PROTECT (sus_status_p = allocVector (STRSXP, 1));
  pc++;
  REAL (sim_time_p)[0] = sim_time;
  INTEGER (infector_id_p)[0] = q->ego_id;
  INTEGER (strain_p)[0] = strain_id;
  INTEGER (phylo_id_p)[0] = phylo_id;
  SET_STRING_ELT (infected_status_p, 0, mkChar ("infectious"));
  SET_STRING_ELT (sus_status_p, 0, mkChar ("susceptible"));
  g = netSetVertexAttrib (g, "time_infected", sim_time_p, altar_id);
  g = netSetVertexAttrib (g, "infection_history", strain_p, altar_id);
  g = netSetVertexAttrib (g, "status", infected_status_p, altar_id);
  g = netSetVertexAttrib (g, "infector_id", infector_id_p, altar_id);
  g = netSetVertexAttrib (g, "phylo_id", phylo_id_p, altar_id);

  if (phylo_id > MAXIMMUN)
    {
      error("%s: %d: MAXIMMUN is not larger than phylo_id\n");
    }

  /*Mutate if specified*/
  if (two_port == 1)
    {
      phylo_id = next_phylo_id;
      INTEGER (phylo_id_p)[0] = phylo_id;
      next_phylo_id++;
      g = netSetVertexAttrib (g, "phylo_id", phylo_id_p, altar_id);
      INTEGER (strain_p)[0] = ++strain_id;
      g = netSetVertexAttrib (g, "infection_history", strain_p, altar_id);
    }
   else
     {
       g = netSetVertexAttrib (g, "infection_history", strain_p, altar_id);
       g = netSetVertexAttrib (g, "phylo_id", phylo_id_p, altar_id);
     }
  
  /*update case load counter */
  case_load++;
  
  /*Add event of host stimulating immunity to the table */
  k = (struct key *) malloc (sizeof (struct key));
  if (NULL == k)
    {
      hashtable_destroy (h, 1);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
    }
  k->ego_id = altar_id;
  k->altar_id = 0;
  k->function_id = WIMMSTIM;
  k->two_port = 0;

  v = (struct value *) malloc (sizeof (struct value));
  if (NULL == v)
    {
      hashtable_destroy (h, 1);
      free (k);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
    }
  v->phylo_id = phylo_id;
  v->rate = efficiency * adaptive_immunity_levels[phylo_id];
  total_rates += v->rate;

  if (!insert_some (h, k, v))
    {
      hashtable_destroy (h, 1);
      free (k);
      error ("%s:%d: Insertion into hash table failed\n", __FILE__, __LINE__);
    }

  /*Add event of host recovering to the table */
  k = (struct key *) malloc (sizeof (struct key));
  if (NULL == k)
    {
      hashtable_destroy (h, 1);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
    }
  k->ego_id = altar_id;
  k->altar_id = 0;
  k->function_id = EVNTREC;
  k->two_port = 0;

  v = (struct value *) malloc (sizeof (struct value));
  if (NULL == v)
    {
      hashtable_destroy (h, 1);
      free (k);
      error ("%s: %d: Ran out of memory allocating a key\n",
	     __FILE__, __LINE__);
    }
  v->rate = sim_gamma;
  total_rates += v->rate;

  if (!insert_some (h, k, v))
    {
      hashtable_destroy (h, 1);
      free (k);
      error ("%s:%d: Insertion into hash table failed\n", __FILE__, __LINE__);
    }
  

  /* Add event of infected host infecting susceptible neighbor cells */
  PROTECT (hood = netGetNeighborhood (g, altar_id, "out", 1));
  pc++;
  PROTECT (val = getListElement (g, "val"));
  pc++;
  for (i = 0; i < length (hood); i++)
    {
      altar_nb_id = INTEGER (hood)[i];
      PROTECT (vall =
	       getListElement (VECTOR_ELT (val, altar_nb_id - 1), "status"));
      pc++;
     if (vecEq (vall, sus_status_p))
	{
	  k = (struct key *) malloc (sizeof (struct key));
	  if (NULL == k)
	    {
	      hashtable_destroy (h, 1);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
	  k->ego_id = altar_id;
	  k->altar_id = altar_nb_id;
	  k->function_id = EVNTINF;
	  k->two_port = 0;

	  v = (struct value *) malloc (sizeof (struct value));
	  if (NULL == v)
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
		v->strain_id = strain_id;
		v->phylo_id = phylo_id;
/*This does strange things for some reason, so assign RHS to ints above
 * v->strain_id = p->strain_id;
 *	  	  v->phylo_id = p->phylo_id;
 */
	  PROTECT (vall =
		   getListElement (VECTOR_ELT (val, altar_nb_id - 1),
				   "immune_memory"));
	  pc++;
	  distance = p->strain_id + fabs (INTEGER (vall)[0]);
	  tau = exp (-a * distance);
	  rate = sim_beta * (1 - tau);
	  v->rate = rate;
	  total_rates += v->rate;

	  if (!insert_some (h, k, v))
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s:%d: Insertion into hash table failed\n", __FILE__,
		     __LINE__);
	    }
	}
    }
  
  /* Add event of infected host infecting susceptible neighbor cells with a mutant*/
  for (i = 0; i < length (hood); i++)
    {

      altar_nb_id = INTEGER (hood)[i];
      PROTECT (vall =
	       getListElement (VECTOR_ELT (val, altar_nb_id - 1), "status"));
      pc++;
     if (vecEq (vall, sus_status_p))
	{
	  k = (struct key *) malloc (sizeof (struct key));
	  if (NULL == k)
	    {
	      hashtable_destroy (h, 1);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
	  k->ego_id = altar_id;
	  k->altar_id = altar_nb_id;
	  k->function_id = EVNTINF;
	  k->two_port = 1;

	  v = (struct value *) malloc (sizeof (struct value));
	  if (NULL == v)
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
		v->strain_id = strain_id;
		v->phylo_id = phylo_id;
/*This does strange things for some reason, so assign RHS to ints above
 * v->strain_id = p->strain_id;
 *	  	  v->phylo_id = p->phylo_id;
 */
	  PROTECT (vall =
		   getListElement (VECTOR_ELT (val, altar_nb_id - 1),
				   "immune_memory"));
	  pc++;
	  distance = p->strain_id + fabs (INTEGER (vall)[0]);
	  tau = exp (-a * distance);
	  rate = sim_beta * (1 - tau) * sim_mu;
	  v->rate = rate;
	  total_rates += v->rate;

	  if (!insert_some (h, k, v))
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s:%d: Insertion into hash table failed\n", __FILE__,
		     __LINE__);
	    }
	}
    }

  /* Remove event of infected host geting infected by infectious neighbors, either with or without a mutant */
  /* For some reason, doing this before adding the infection events switches
   * the ego_id in q, so it is done last*/
  for (i = 0; i < length (hood); i++)
    {
      altar_nb_id = INTEGER (hood)[i];
      PROTECT (vall =
	       getListElement (VECTOR_ELT (val, altar_nb_id - 1), "status"));
      pc++;
      k = (struct key *) malloc (sizeof (struct key));
      if (NULL == k)
        {
          hashtable_destroy (h, 1);
          error ("%s: %d: Ran out of memory allocating a key\n",
                 __FILE__, __LINE__);
        }
      if (vecEq (vall, infected_status_p))
	{
	  k->function_id = EVNTINF;
	  k->ego_id = altar_nb_id;
	  k->altar_id = altar_id;
          k->two_port = 0;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
	      free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free(found);
	  
          k->function_id = EVNTINF;
	  k->ego_id = altar_nb_id;
	  k->altar_id = altar_id;
          k->two_port = 1;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
	      free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free(found);
	}
    }
  UNPROTECT (pc);
  return (0);
}

int
wimmstim (struct key *q, struct value *p)
{
  /* host EGO_ID is burst at time SIM_TIME, and stimulates adaptive immunity */

  SEXP sim_time_p, status_p, vall, val, infected_status_p;
  SEXP hood, infector_id_p, sus_status_p, rec_status_p;
  SEXP strain_id_p;
  int pc = 0, i, ego_nb_id, ego_id, phylo_id, n;
  struct key *k;

  ego_id = q->ego_id;
  phylo_id = p->phylo_id;
  adaptive_immunity_levels[phylo_id] += conversion_factor;
/*  Rprintf("adaptive_immunity_levels[%d] = %g\n", phylo_id, adaptive_immunity_levels[phylo_id]);*/

  /*Increase the immune stimulation rates with new level of immunity */
  n = netNetSize(g);
  k = (struct key *)malloc(sizeof(struct key));
      if (NULL == k) {
          printf("ran out of memory allocating a key\n");
          return 1;
      }  
   for (i = 1; i <= n; i++)
     {
       k->ego_id = i;
       k->altar_id = 0;
       k->function_id = WIMMSTIM;
       k->two_port = 0;

       if (NULL != (found = search_some(h,k))) {
           if (found->phylo_id == phylo_id){
               total_rates -= found->rate;
               found->rate = adaptive_immunity_levels[phylo_id] * efficiency;
               total_rates += found->rate;
           }
       }
     }



  /*Add and update vertex attributes */
  PROTECT (sim_time_p = allocVector (REALSXP, 1));
  pc++;
  PROTECT (infected_status_p = allocVector (STRSXP, 1));
  pc++;
  PROTECT (sus_status_p = allocVector (STRSXP, 1));
  pc++;
  PROTECT (rec_status_p = allocVector (STRSXP, 1));
  pc++;
  REAL (sim_time_p)[0] = sim_time;
  SET_STRING_ELT (infected_status_p, 0, mkChar ("infectious"));
  SET_STRING_ELT (sus_status_p, 0, mkChar ("susceptible"));
  SET_STRING_ELT (rec_status_p, 0, mkChar ("recovered"));
  g = netSetVertexAttrib (g, "time_recovered", sim_time_p, q->ego_id);
  g = netSetVertexAttrib (g, "status", rec_status_p, q->ego_id);

  /*update case load counter */
  case_load--;

  /*Remove event of the newly recovered host again recovering from the dl list */
	  k = (struct key *) malloc (sizeof (struct key));
	  if (NULL == k)
	    {
	      hashtable_destroy (h, 1);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
	  k->ego_id = ego_id;
	  k->altar_id = 0;
	  k->function_id = EVNTREC;
          k->two_port = 0;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
	      free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free (found);
  
  /*Remove event of the newly recovered host again stimulating immunity */
	  k = (struct key *) malloc (sizeof (struct key));
	  if (NULL == k)
	    {
	      hashtable_destroy (h, 1);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
	  k->ego_id = ego_id;
	  k->altar_id = 0;
	  k->function_id = WIMMSTIM;
          k->two_port = 0;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
	      free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free (found);

  /* Remove event of newly recovered host infecting suceptible neighbors from pointer array and dl list */
  PROTECT (hood = netGetNeighborhood (g, ego_id, "out", 1));
  pc++;
  PROTECT (val = getListElement (g, "val"));
  pc++;				//get list of vertex attributes
  PROTECT (strain_id_p =
	   getListElement (VECTOR_ELT (val, ego_id - 1),
			   "infection_history"));
  pc++;
  for (i = 0; i < length (hood); i++)
    {
      ego_nb_id = INTEGER (hood)[i];
      PROTECT (vall =
	       getListElement (VECTOR_ELT (val, ego_nb_id - 1), "status"));
      pc++;
      if (vecEq (vall, sus_status_p))
	{
	  k = (struct key *) malloc (sizeof (struct key));
          if (NULL == k)
            {
              hashtable_destroy (h, 1);
              error ("%s: %d: Ran out of memory allocating a key\n",
                     __FILE__, __LINE__);
            }
	  k->ego_id = ego_id;
	  k->altar_id = ego_nb_id;
	  k->function_id = EVNTINF;
          k->two_port = 0;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free (found);
	  
          k = (struct key *) malloc (sizeof (struct key));
          if (NULL == k)
            {
              hashtable_destroy (h, 1);
              error ("%s: %d: Ran out of memory allocating a key\n",
                     __FILE__, __LINE__);
            }
	  k->ego_id = ego_id;
	  k->altar_id = ego_nb_id;
	  k->function_id = EVNTINF;
          k->two_port = 1;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free (found);
	}

    }
  UNPROTECT (pc);
  return (0);
}


int
wrecover (struct key *q, struct value *p)
{
  /* host EGO_ID recovers at time SIM_TIME */

  SEXP sim_time_p, status_p, vall, val, infected_status_p;
  SEXP hood, infector_id_p, sus_status_p, rec_status_p;
  SEXP strain_id_p;
  int pc = 0, i, ego_nb_id, ego_id;
  struct key *k;

  ego_id = q->ego_id;

  /*Add and update vertex attributes */
  PROTECT (sim_time_p = allocVector (REALSXP, 1));
  pc++;
  PROTECT (infected_status_p = allocVector (STRSXP, 1));
  pc++;
  PROTECT (sus_status_p = allocVector (STRSXP, 1));
  pc++;
  PROTECT (rec_status_p = allocVector (STRSXP, 1));
  pc++;
  REAL (sim_time_p)[0] = sim_time;
  SET_STRING_ELT (infected_status_p, 0, mkChar ("infectious"));
  SET_STRING_ELT (sus_status_p, 0, mkChar ("susceptible"));
  SET_STRING_ELT (rec_status_p, 0, mkChar ("recovered"));
  g = netSetVertexAttrib (g, "time_recovered", sim_time_p, q->ego_id);
  g = netSetVertexAttrib (g, "status", rec_status_p, q->ego_id);

  /*update case load counter */
  case_load--;
  
  /*Remove event of the newly recovered host again stimulating immunity */
	  k = (struct key *) malloc (sizeof (struct key));
	  if (NULL == k)
	    {
	      hashtable_destroy (h, 1);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
	  k->ego_id = ego_id;
	  k->altar_id = 0;
	  k->function_id = WIMMSTIM;
          k->two_port = 0;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
	      free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free (found);

  /*Remove event of the newly recovered host again recovering from the dl list */
	  k = (struct key *) malloc (sizeof (struct key));
	  if (NULL == k)
	    {
	      hashtable_destroy (h, 1);
	      error ("%s: %d: Ran out of memory allocating a key\n",
		     __FILE__, __LINE__);
	    }
	  k->ego_id = ego_id;
	  k->altar_id = 0;
	  k->function_id = EVNTREC;
          k->two_port = 0;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
	      free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free (found);

  /* Remove event of newly recovered host infecting suceptible neighbors from pointer array and dl list */
  PROTECT (hood = netGetNeighborhood (g, ego_id, "out", 1));
  pc++;
  PROTECT (val = getListElement (g, "val"));
  pc++;				//get list of vertex attributes
  PROTECT (strain_id_p =
	   getListElement (VECTOR_ELT (val, ego_id - 1),
			   "infection_history"));
  pc++;
  for (i = 0; i < length (hood); i++)
    {
      ego_nb_id = INTEGER (hood)[i];
      PROTECT (vall =
	       getListElement (VECTOR_ELT (val, ego_nb_id - 1), "status"));
      pc++;
      if (vecEq (vall, sus_status_p))
	{
	  k = (struct key *) malloc (sizeof (struct key));
          if (NULL == k)
            {
              hashtable_destroy (h, 1);
              error ("%s: %d: Ran out of memory allocating a key\n",
                     __FILE__, __LINE__);
            }
	  k->ego_id = ego_id;
	  k->altar_id = ego_nb_id;
	  k->function_id = EVNTINF;
          k->two_port = 0;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free (found);
	  
          k = (struct key *) malloc (sizeof (struct key));
          if (NULL == k)
            {
              hashtable_destroy (h, 1);
              error ("%s: %d: Ran out of memory allocating a key\n",
                     __FILE__, __LINE__);
            }
	  k->ego_id = ego_id;
	  k->altar_id = ego_nb_id;
	  k->function_id = EVNTINF;
          k->two_port = 1;
	  if (NULL == (found = remove_some (h, k)))
	    {
	      hashtable_destroy (h, 1);
              free (k);
	      error ("%s:%d: Key not found for removal\n",
		     __FILE__, __LINE__);
	    }
          total_rates -= found->rate;
          free (found);
	}

    }
  UNPROTECT (pc);
  return (0);
}
