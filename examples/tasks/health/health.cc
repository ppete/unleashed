/*
 * Code modified to run with the UCL task framework
 *
 * The Unleashed Concurrency Library's Task testing framework
 *
 * Copyright (c) 2019, LLNL
 *
 * - ported to support multiple task-parallel models (incl. Cilk, UCL)
 * - ported code partially to C++11
 * - replaced use of locks with lock-free insert
 */

/**********************************************************************************************/
/*  This program is part of the Barcelona OpenMP Tasks Suite                                  */
/*  Copyright (C) 2009 Barcelona Supercomputing Center - Centro Nacional de Supercomputacion  */
/*  Copyright (C) 2009 Universitat Politecnica de Catalunya                                   */
/**********************************************************************************************/

/* OLDEN parallel C for dynamic structures: compiler, runtime system
 * and benchmarks
 *
 * Copyright (C) 1994-1996 by Anne Rogers (amr@cs.princeton.edu) and
 * Martin Carlisle (mcc@cs.princeton.edu)
 * ALL RIGHTS RESERVED.
 *
 * OLDEN is distributed under the following conditions:
 *
 * You may make copies of OLDEN for your own use and modify those copies.
 *
 * All copies of OLDEN must retain our names and copyright notice.
 *
 * You may not sell OLDEN or distribute OLDEN in conjunction with a
 * commercial product or service without the expressed written consent of
 * Anne Rogers and Martin Carlisle.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE.
 *
 */


/*******************************************************************
 *  Health.c : Model of the Colombian Health Care System           *
 *******************************************************************/
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <cassert>
#include <atomic>
#include <chrono>

#include "../common/common-includes.hpp"
#include "../common/bots.hpp"


#include "health.h"


//~ #define CUTOFF (2)

#define CUTOFF (0)

/* global variables */
int sim_level;
int sim_cities;
int sim_population_ratio;
int sim_time;
int sim_assess_time;
int sim_convalescence_time;
int32_t sim_seed;
float sim_get_sick_p;
float sim_convalescence_p;
float sim_realloc_p;
int sim_pid = 0;

int res_population;
int res_hospitals;
int res_personnel;
int res_checkin;
int res_village;
int res_waiting;
int res_assess;
int res_inside;
float res_avg_stay;

/**********************************************************
 * Handles math routines for health.c                     *
 **********************************************************/
float my_rand(int32_t *seed)
{
   int32_t k;
   int32_t idum = *seed;

   idum ^= MASK;
   k = idum / IQ;
   idum = IA * (idum - k * IQ) - IR * k;
   idum ^= MASK;
   if (idum < 0) idum  += IM;
   *seed = idum * IM;
   return (float) AM * idum;
}

/********************************************************************
 * Handles lists.                                                   *
 ********************************************************************/
void addList(std::atomic<Patient*>* list, Patient* patient)
{
  Patient* pred = list->load(std::memory_order_relaxed);

  patient->forward.store(nullptr, std::memory_order_relaxed);

  if (pred == nullptr)
  {
    list->store(patient, std::memory_order_relaxed);

    patient->back.store(nullptr, std::memory_order_relaxed);
  }
  else
  {
    Patient* curr = pred->forward.load(std::memory_order_relaxed);
    while (curr != nullptr)
    {
      pred = curr;
      curr = pred->forward.load(std::memory_order_relaxed);
    }

    pred->forward.store(patient, std::memory_order_relaxed);
    patient->back.store(pred, std::memory_order_relaxed);
  }
}


// simple lock-free add to replace use of locks in original code
//   - ONLY SAFE with other concurrent addList_c calls,
//   - BUT NOT with concurrent removes.
void c_addList(std::atomic<Patient*>* list, Patient* patient)
{
  Patient* pred = list->load(std::memory_order_relaxed);

  patient->forward.store(nullptr, std::memory_order_relaxed);

  if (pred == nullptr)
  {
    patient->back = NULL;

    if (list->compare_exchange_strong( pred,
                                       patient,
                                       std::memory_order_relaxed,
                                       std::memory_order_relaxed
                                     )
       )
      return;
  }

  assert(pred != nullptr);

  Patient* curr = pred->forward.load(std::memory_order_relaxed);

  while (curr != nullptr) { pred = curr; curr = pred->forward.load(std::memory_order_relaxed); }

  patient->back.store(pred, std::memory_order_relaxed);

  while (!pred->forward.compare_exchange_strong( curr,
                                                 patient,
                                                 std::memory_order_relaxed,
                                                 std::memory_order_relaxed
                                               )
        )
  {
    pred = curr;
    patient->back.store(pred, std::memory_order_relaxed);
    curr = nullptr;
  }
}


void removeList(std::atomic<Patient*> *list, Patient *patient)
{
  Patient* pred = patient->back.load(std::memory_order_relaxed);
  Patient* succ = patient->forward.load(std::memory_order_relaxed);

  if (pred != nullptr)
    pred->forward.store(succ, std::memory_order_relaxed);
  else
    list->store(succ, std::memory_order_relaxed);

  if (succ != nullptr)
    succ->back.store(pred, std::memory_order_relaxed);
}

/**********************************************************************/
void allocate_village(Village** capital, Village* back,Village* next, int level, int32_t vid)
{
  if (level == 0) { *capital = NULL; return; }

  int personnel = (int) pow(2, level);
  int population = personnel * sim_population_ratio;

  /* Allocate Village + Hospital */
  *capital = new Village(vid, back, next, level, vid * (IQ + sim_seed), personnel);
  /* Initialize Village */
  (*capital)->back  = back;
  (*capital)->next  = next;
  (*capital)->level = level;
  (*capital)->id    = vid;
  (*capital)->seed  = vid * (IQ + sim_seed);

  for (int i=0; i<population; i++)
  {
    Patient* patient = new Patient(sim_pid++, (*capital)->seed, *capital);

    // changes seed for capital:
    my_rand(&((*capital)->seed));
    addList(&((*capital)->population), patient);
  }

  // Create Cities (lower level)
  Village* current = nullptr;

  for (int i = sim_cities; i>0; i--)
  {
    Village* inext = current;

    allocate_village(&current, *capital, inext, level-1, (vid * (int32_t) sim_cities)+ (int32_t) i);
  }

  (*capital)->forward = current;
}

/**********************************************************************/
struct Results get_results(struct Village *village)
{
   struct Village *vlist;
   struct Patient *p;
   struct Results t_res, p_res;

   t_res.hosps_number     = 0.0;
   t_res.hosps_personnel  = 0.0;
   t_res.total_patients   = 0.0;
   t_res.total_in_village = 0.0;
   t_res.total_waiting    = 0.0;
   t_res.total_assess     = 0.0;
   t_res.total_inside     = 0.0;
   t_res.total_hosps_v    = 0.0;
   t_res.total_time       = 0.0;

   if (village == NULL) return t_res;

   /* Traverse village hierarchy (lower level first)*/
   vlist = village->forward;
   while(vlist)
   {
      p_res = get_results(vlist);
      t_res.hosps_number     += p_res.hosps_number;
      t_res.hosps_personnel  += p_res.hosps_personnel;
      t_res.total_patients   += p_res.total_patients;
      t_res.total_in_village += p_res.total_in_village;
      t_res.total_waiting    += p_res.total_waiting;
      t_res.total_assess     += p_res.total_assess;
      t_res.total_inside     += p_res.total_inside;
      t_res.total_hosps_v    += p_res.total_hosps_v;
      t_res.total_time       += p_res.total_time;
      vlist = vlist->next;
   }
   t_res.hosps_number     += 1.0;
   t_res.hosps_personnel  += village->hosp.personnel;

   // Patients in the village
   p = village->population;
   while (p != NULL)
   {
      t_res.total_patients   += 1.0;
      t_res.total_in_village += 1.0;
      t_res.total_hosps_v    += ((float)p->hosps_visited);
      t_res.total_time       += (float)(p->time);
      p = p->forward;
   }
   // Patients in hospital: waiting
   p = village->hosp.waiting;
   while (p != NULL)
   {
      t_res.total_patients += 1.0;
      t_res.total_waiting  += 1.0;
      t_res.total_hosps_v  += (float)(p->hosps_visited);
      t_res.total_time     += (float)(p->time);
      p = p->forward;
   }
   // Patients in hospital: assess
   p = village->hosp.assess;
   while (p != NULL)
   {
      t_res.total_patients += 1.0;
      t_res.total_assess   += 1.0;
      t_res.total_hosps_v  += (float)(p->hosps_visited);
      t_res.total_time     += (float)(p->time);
      p = p->forward;
   }
   // Patients in hospital: inside
   p = village->hosp.inside;
   while (p != NULL)
   {
      t_res.total_patients += 1.0;
      t_res.total_inside   += 1.0;
      t_res.total_hosps_v  += (float)(p->hosps_visited);
      t_res.total_time     += (float)(p->time);
      p = p->forward;
   }

   return t_res;
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void check_patients_inside(struct Village *village)
{
   Patient *list = village->hosp.inside.load(std::memory_order_relaxed);

   while (list != NULL)
   {
      Patient *p = list;
      list = list->forward.load(std::memory_order_relaxed);
      p->time_left--;
      if (p->time_left == 0)
      {
         village->hosp.free_personnel++;
         removeList(&(village->hosp.inside), p);
         addList(&(village->population), p);
      }
   }
}

/**********************************************************************/
void check_patients_assess_par(struct Village *village)
{
   Patient *list = village->hosp.assess.load(std::memory_order_relaxed);

   while (list != NULL)
   {
      Patient* p = list;
      list = list->forward.load(std::memory_order_relaxed);
      p->time_left--;

      if (p->time_left == 0)
      {
         float rand = my_rand(&(p->seed));
         /* sim_covalescense_p % */
         if (rand < sim_convalescence_p)
         {
            rand = my_rand(&(p->seed));
            /* !sim_realloc_p % or root hospital */
            if (rand > sim_realloc_p || village->level == sim_level)
            {
               removeList(&(village->hosp.assess), p);
               addList(&(village->hosp.inside), p);
               p->time_left = sim_convalescence_time;
               p->time += p->time_left;
            }
            else /* move to upper level hospital !!! */
            {
               village->hosp.free_personnel++;
               removeList(&(village->hosp.assess), p);
               c_addList(&(village->back->hosp.realloc), p);
            }
         }
         else /* move to village */
         {
            village->hosp.free_personnel++;
            removeList(&(village->hosp.assess), p);

            // \pp should this not be
            // c_addList(&(p->home_village->population), p);
            addList(&(village->population), p);
         }
      }
   }
}

/**********************************************************************/
void check_patients_waiting(struct Village *village)
{
   struct Patient *list = village->hosp.waiting;
   struct Patient *p;

   while (list != NULL)
   {
      p = list;
      list = list->forward.load(std::memory_order_relaxed);
      if (village->hosp.free_personnel > 0)
      {
         village->hosp.free_personnel--;
         p->time_left = sim_assess_time;
         p->time += p->time_left;
         removeList(&(village->hosp.waiting), p);
         addList(&(village->hosp.assess), p);
      }
      else
      {
         p->time++;
      }
   }
}

/**********************************************************************/
void check_patients_realloc(struct Village *village)
{
   while (village->hosp.realloc != NULL)
   {
      Patient *s = village->hosp.realloc.load(std::memory_order_relaxed);
      Patient *p = s;

      while (p != nullptr)
      {
         if (p->id < s->id) s = p;
         p = p->forward.load(std::memory_order_relaxed);
      }

      removeList(&(village->hosp.realloc), s);
      put_in_hosp(&(village->hosp), s);
   }
}

/**********************************************************************/
void check_patients_population(struct Village *village)
{
   Patient *list = village->population.load(std::memory_order_relaxed);

   while (list != NULL)
   {
      Patient *p = list;
      list = list->forward.load(std::memory_order_relaxed);
      /* randomize in patient */
      float rand = my_rand(&(p->seed));
      if (rand < sim_get_sick_p)
      {
         removeList(&(village->population), p);
         put_in_hosp(&(village->hosp), p);
      }
   }
}

/**********************************************************************/
void put_in_hosp(struct Hosp *hosp, struct Patient *patient)
{
   (patient->hosps_visited)++;

   if (hosp->free_personnel > 0)
   {
      hosp->free_personnel--;
      addList(&(hosp->assess), patient);
      patient->time_left = sim_assess_time;
      patient->time += patient->time_left;
   }
   else
   {
      addList(&(hosp->waiting), patient);
   }
}

#if WOMP_VERSION

/**********************************************************************/
void sim_village_par(struct Village *village)
{
   struct Village *vlist;

   // lowest level returns nothing
   // only for sim_village first call with village = NULL
   // recursive call cannot occurs
   if (village == NULL) return;

   /* Traverse village hierarchy (lower level first)*/
   vlist = village->forward;
   while(vlist)
   {
      #pragma omp task untied if (!CUTOFF || (sim_level - village->level) < CUTOFF)
      sim_village_par(vlist);
      vlist = vlist->next;
   }

   /* Uses lists v->hosp->inside, and v->return */
   check_patients_inside(village);

   /* Uses lists v->hosp->assess, v->hosp->inside, v->population and (v->back->hosp->realloc) !!! */
   check_patients_assess_par(village);

   /* Uses lists v->hosp->waiting, and v->hosp->assess */
   check_patients_waiting(village);

#pragma omp taskwait

   /* Uses lists v->hosp->realloc, v->hosp->asses and v->hosp->waiting */
   check_patients_realloc(village);

   /* Uses list v->population, v->hosp->asses and v->h->waiting */
   check_patients_population(village);
}

/**********************************************************************/
void sim_village_main_par(Village *top, size_t numthreads)
{
#pragma omp parallel num_threads(numthreads)
#pragma omp single
#pragma omp task untied
   for (long i = 0; i < sim_time; i++) sim_village_par(top);
}

#endif /* WOMP_VERSION */


#if OMP_VERSION

/**********************************************************************/
void sim_village_par(Village* village, ucl::continuation outer)
{
  typedef ucl::continuation::counter_t counter_t;

  // lowest level returns nothing
  // only for sim_village first call with village = NULL
  // recursive call cannot occurs
  if (village == NULL) return;

  Village*          vlist = village->forward;
  counter_t         cont_external_cnt = ucl::continuation::counter_max();
  ucl::continuation cont( [ village, outercont = std::move(outer) ]() -> void
                          {
                            /* Uses lists v->hosp->realloc, v->hosp->asses and v->hosp->waiting */
                            check_patients_realloc(village);

                            /* Uses list v->population, v->hosp->asses and v->h->waiting */
                            check_patients_population(village);
                          },
                          cont_external_cnt
                        );

  /* Traverse village hierarchy (lower level first)*/
  vlist = village->forward;
  while(vlist)
  {
    if (!CUTOFF || (sim_level - village->level) < CUTOFF)
    {
      --cont_external_cnt;

      ucl::count_fun_base* ptr = cont.omp_ptr();

      #pragma omp task untied firstprivate(vlist, ptr)
      sim_village_par(vlist, ucl::continuation(ptr));

      assert(cont_external_cnt > 1);
    }
    else
    {
      sim_village_par(vlist, ucl::continuation());
    }

    vlist = vlist->next;
  }

  /* Uses lists v->hosp->inside, and v->return */
  check_patients_inside(village);

  /* Uses lists v->hosp->assess, v->hosp->inside, v->population and (v->back->hosp->realloc) !!! */
  check_patients_assess_par(village);

  /* Uses lists v->hosp->waiting, and v->hosp->assess */
  check_patients_waiting(village);

  if (cont_external_cnt == ucl::continuation::counter_max())
  {
    cont.release_local(cont_external_cnt);
  }
  else
  {
    cont.release_external(cont_external_cnt);
  }
}

/**********************************************************************/
void sim_village_main_par(Village *top, size_t numthreads)
{
  #pragma omp parallel num_threads(numthreads)
  #pragma omp single
  for (long i = 0; i < sim_time; i++)
  {
    #pragma omp taskgroup
    {
      #pragma omp task untied firstprivate(top)
      sim_village_par(top, ucl::continuation());
    }
  }
}

#endif /* OMP_VERSION */


#if UCL_VERSION

/**********************************************************************/

struct sim_task
{
  sim_task()
  : village(nullptr), after()
  {}

  sim_task(Village* v, ucl::continuation&& cont)
  : village(v), after(std::move(cont))
  {}

  sim_task(sim_task&& other)
  : village(other.village), after(std::move(other.after))
  {}

  sim_task& operator=(sim_task&& other)
  {
    sim_task tmp(std::move(*this));

    village = other.village;
    after   = std::move(other.after);
    return *this;
  }

  sim_task(const sim_task& other)
  : village(other.village), after(other.after)
  {}

  sim_task& operator=(const sim_task& other) = delete;

  Village*          village;
  ucl::continuation after;
};

struct simulator
{
  template <class Pool>
  ucl::Void operator()(Pool& pool, sim_task t)
  {
    typedef ucl::continuation::counter_t counter_t;

    // lowest level returns nothing
    // only for sim_village first call with village = NULL
    // recursive call cannot occurs
    if (t.village == NULL) return ucl::Void();

    /* Traverse village hierarchy (lower level first) */
    Village*          village           = t.village;
    Village*          vlist             = village->forward;
    counter_t         cont_external_cnt = ucl::continuation::counter_max();

    ucl::continuation cont( [ tt = std::move(t) ]() -> void
                            {
                              /* Uses lists v->hosp->realloc, v->hosp->asses and v->hosp->waiting */
                              check_patients_realloc(tt.village);

                              /* Uses list v->population, v->hosp->asses and v->h->waiting */
                              check_patients_population(tt.village);
                            },
                            cont_external_cnt
                          );

     while (vlist)
     {
        if (!CUTOFF || (sim_level - village->level) < CUTOFF)
        {
          --cont_external_cnt;
          pool.enq(sim_task(vlist, cont.copy_external_count()));
          assert(cont_external_cnt > 1);
        }
        else
        {
          (*this)(pool, sim_task(vlist, ucl::continuation()));
        }

        vlist = vlist->next;
     }

     /* Uses lists v->hosp->inside, and v->return */
     check_patients_inside(village);

     /* Uses lists v->hosp->assess, v->hosp->inside, v->population and (v->back->hosp->realloc) !!! */
     check_patients_assess_par(village);

     /* Uses lists v->hosp->waiting, and v->hosp->assess */
     check_patients_waiting(village);

     if (cont_external_cnt == ucl::continuation::counter_max())
     {
       cont.release_local(cont_external_cnt);
     }
     else
     {
       cont.release_external(cont_external_cnt);
     }

     return ucl::Void();
  }
};

/**********************************************************************/
void sim_village_main_par(Village* top, size_t numthreads)
{
#if 1
  long                        i = 0;
  ucl::default_pool<sim_task> pool{numthreads};
  std::function<void()>       action{ [&i, &pool, &action, top]() -> void
                                      {
                                        if (++i < sim_time)
                                          pool.enq(sim_task{top, ucl::continuation{action, 1}});
                                      }
                                    };
  pool.enq(sim_task{top, ucl::continuation{action, 1}});

  ucl::execute_pool(simulator(), pool);

#else
  // \todo
  // there seems to be a performance bug when tasks pools are frequently created
  //   and their threads pinned to a specific processor.
  //   The performance degrades by >10% for small tasks. Interestingly, the
  //    degradation does not occur without thread pinning (using generic_arch).

  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  for (long i = 0; i < sim_time; i++)
  {
    time_point starttime   = std::chrono::system_clock::now();
    ucl::execute_tasks_x(numthreads, simulator(), sim_task{ top, ucl::continuation() });
    time_point endtime     = std::chrono::system_clock::now();

    std::cerr << i << ": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count()
              << std::endl;
  }
#endif
}

#endif /* UCL_VERSION */


#if CILK_VERSION

/**********************************************************************/
void sim_village_par(struct Village *village)
{
   struct Village *vlist;

   // lowest level returns nothing
   // only for sim_village first call with village = NULL
   // recursive call cannot occurs
   if (village == NULL) return;

   /* Traverse village hierarchy (lower level first)*/
   vlist = village->forward;
   while(vlist)
   {
      if (!CUTOFF || (sim_level - village->level) < CUTOFF)
      {
        cilk_spawn sim_village_par(vlist);
      }
      else
      {
        sim_village_par(vlist);
      }

      vlist = vlist->next;
   }

   /* Uses lists v->hosp->inside, and v->return */
   check_patients_inside(village);

   /* Uses lists v->hosp->assess, v->hosp->inside, v->population and (v->back->hosp->realloc) !!! */
   check_patients_assess_par(village);

   /* Uses lists v->hosp->waiting, and v->hosp->assess */
   check_patients_waiting(village);

   cilk_sync;

   /* Uses lists v->hosp->realloc, v->hosp->asses and v->hosp->waiting */
   check_patients_realloc(village);

   /* Uses list v->population, v->hosp->asses and v->h->waiting */
   check_patients_population(village);
}

/**********************************************************************/
void sim_village_main_par(Village *top, size_t numthreads)
{
  set_cilk_workers(numthreads);

  for (long i = 0; i < sim_time; i++) sim_village_par(top);
}

#endif /* CILK_VERSION */


/**********************************************************************/
void my_print(struct Village *village)
{
   struct Village *vlist;
   struct Patient *plist;

   if (village == NULL) return;

   /* Traverse village hierarchy (lower level first)*/
   vlist = village->forward;
   while(vlist) {
      my_print(vlist);
      vlist = vlist->next;
   }

   plist = village->population;

   while (plist != NULL) {
      bots_debug("[pid:%d]",plist->id);
      plist = plist->forward;
   }
   bots_debug("[vid:%d]\n",village->id);

}
/**********************************************************************/
void read_input_data(const char* filename)
{
   FILE *fin;
   int res;
   long int xseed;

   if ((fin = fopen(filename, "r")) == NULL) {
      bots_message("Could not open sequence file (%s)\n", filename);
      exit (-1);
   }

   res = fscanf(fin,"%d %d %d %d %d %d %ld %f %f %f %d %d %d %d %d %d %d %d %f",
             &sim_level,
             &sim_cities,
             &sim_population_ratio,
             &sim_time,
             &sim_assess_time,
             &sim_convalescence_time,
             &xseed,
             &sim_get_sick_p,
             &sim_convalescence_p,
             &sim_realloc_p,
             &res_population,
             &res_hospitals,
             &res_personnel,
             &res_checkin,
             &res_village,
             &res_waiting,
             &res_assess,
             &res_inside,
             &res_avg_stay
   );

   sim_seed = xseed;
   if ( res == EOF ) {
      bots_message("Bogus input file (%s)\n", filename);
      exit(-1);
   }
   fclose(fin);

      // Printing input data
   bots_message("\n");
   bots_message("Number of levels    = %d\n", (int) sim_level);
   bots_message("Cities per level    = %d\n", (int) sim_cities);
   bots_message("Population ratio    = %d\n", (int) sim_population_ratio);
   bots_message("Simulation time     = %d\n", (int) sim_time);
   bots_message("Assess time         = %d\n", (int) sim_assess_time);
   bots_message("Convalescence time  = %d\n", (int) sim_convalescence_time);
   bots_message("Initial seed        = %d\n", (int) sim_seed);
   bots_message("Get sick prob.      = %f\n", (float) sim_get_sick_p);
   bots_message("Convalescence prob. = %f\n", (float) sim_convalescence_p);
   bots_message("Realloc prob.       = %f\n", (float) sim_realloc_p);
}



int check_village(struct Village *top)
{
   struct Results result = get_results(top);
   int answer = BOTS_RESULT_SUCCESSFUL;

   if (res_population != result.total_patients) answer = BOTS_RESULT_UNSUCCESSFUL;
   if (res_hospitals != result.hosps_number) answer = BOTS_RESULT_UNSUCCESSFUL;
   if (res_personnel != result.hosps_personnel) answer = BOTS_RESULT_UNSUCCESSFUL;
   if (res_checkin != result.total_hosps_v) answer = BOTS_RESULT_UNSUCCESSFUL;
   if (res_village != result.total_in_village) answer = BOTS_RESULT_UNSUCCESSFUL;
   if (res_waiting != result.total_waiting) answer = BOTS_RESULT_UNSUCCESSFUL;
   if (res_assess != result.total_assess) answer = BOTS_RESULT_UNSUCCESSFUL;
   if (res_inside != result.total_inside) answer = BOTS_RESULT_UNSUCCESSFUL;

   bots_message("\n");
   bots_message("Sim. Variables      = expect / result\n");
   bots_message("Total population    = %6d / %6d people\n", (int)   res_population, (int) result.total_patients);
   bots_message("Hospitals           = %6d / %6d people\n", (int)   res_hospitals, (int) result.hosps_number);
   bots_message("Personnel           = %6d / %6d people\n", (int)   res_personnel, (int) result.hosps_personnel);
   bots_message("Check-in's          = %6d / %6d people\n", (int)   res_checkin, (int) result.total_hosps_v);
   bots_message("In Villages         = %6d / %6d people\n", (int)   res_village, (int) result.total_in_village);
   bots_message("In Waiting List     = %6d / %6d people\n", (int)   res_waiting, (int) result.total_waiting);
   bots_message("In Assess           = %6d / %6d people\n", (int)   res_assess, (int) result.total_assess);
   bots_message("Inside Hospital     = %6d / %6d people\n", (int)   res_inside, (int) result.total_inside);
   bots_message("Average Stay        = %6f / %6f u/time\n", (float) res_avg_stay,(float) result.total_time/result.total_patients);

   my_print(top);

   return answer;
}


int main(int argc, char** argv)
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  size_t      num_threads = NUMTHREADS;
  std::string inputfile = "data/medium.input";
  Village*    top = nullptr;

  if (argc > 1) num_threads = aux::as<size_t>(*(argv+1));
  if (argc > 2) inputfile = argv[2];

  std::cout << "loading " << inputfile << std::endl;
  read_input_data(inputfile.c_str());

  allocate_village(&top, NULL, NULL, sim_level, 0);

#if QTHREADS_VERSION
  init_qthreads(num_threads);
#endif /* QTHREADS_VERSION */

  time_point starttime   = std::chrono::system_clock::now();
  sim_village_main_par(top, num_threads);
  time_point endtime     = std::chrono::system_clock::now();
  int        elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  if (check_village(top) != BOTS_RESULT_SUCCESSFUL)
  {
    std::cout << "result mismatch" << std::endl;
    std::cerr << "result mismatch" << std::endl;
    return 1;
  }

  // print timing
  std::cout << "time = " << elapsedtime << "ms" << std::endl;
  std::cerr << elapsedtime << std::endl;

  return 0;
}
