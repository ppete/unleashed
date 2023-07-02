/**********************************************************************************************/
/*  This program is part of the Barcelona OpenMP Tasks Suite                       */
/*  Copyright (C) 2009 Barcelona Supercomputing Center - Centro Nacional de Supercomputacion  */
/*  Copyright (C) 2009 Universitat Politecnica de Catalunya                        */
/*                                                              */
/*  This program is free software; you can redistribute it and/or modify               */
/*  it under the terms of the GNU General Public License as published by               */
/*  the Free Software Foundation; either version 2 of the License, or                 */
/*  (at your option) any later version.                                     */
/*                                                              */
/*  This program is distributed in the hope that it will be useful,                  */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of                   */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    */
/*  GNU General Public License for more details.                               */
/*                                                              */
/*  You should have received a copy of the GNU General Public License                 */
/*  along with this program; if not, write to the Free Software                     */
/*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA        */
/**********************************************************************************************/

#ifndef _HEALTH_H
#define _HEALTH_H

#include <atomic>

/* random defines */
#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

struct Results {
  long hosps_number;
  long hosps_personnel;
  long total_patients;
  long total_in_village;
  long total_waiting;
  long total_assess;
  long total_inside;
  long total_time;
  long total_hosps_v;
};

extern int sim_level;

struct Village;

struct Patient
{
  Patient(int pid, int32_t pseed, Village* home)
  : id(pid), seed(pseed), time(0), time_left(0), hosps_visited(0),
    home_village(home), back(nullptr), forward(nullptr)
  {}

  int                   id;
  int32_t               seed;
  int                   time;
  int                   time_left;
  int                   hosps_visited;
  Village*              home_village;
  std::atomic<Patient*> back;
  std::atomic<Patient*> forward;
};

struct Hosp
{
  explicit
  Hosp(int pers)
  : personnel(pers), free_personnel(pers), waiting(nullptr), assess(nullptr),
    inside(nullptr), realloc(nullptr)
  {}

  int                   personnel;
  int                   free_personnel;
  std::atomic<Patient*> waiting;
  std::atomic<Patient*> assess;
  std::atomic<Patient*> inside;
  std::atomic<Patient*> realloc;
};

struct Village
{
  Village(int vid, Village* pback, Village* pnext, int lv, int32_t pseed, int hosp_pers)
  : id(vid), back(pback), next(pnext), forward(nullptr), population(nullptr),
    level(lv), seed(pseed), hosp(hosp_pers)
  {}

  int                   id;
  Village*              back;
  Village*              next;
  Village*              forward;
  std::atomic<Patient*> population;
  int                   level;
  int32_t               seed;
  Hosp                  hosp;
};

float my_rand(int32_t *seed);

struct Patient *generate_patient(struct Village *village);
void put_in_hosp(struct Hosp *hosp, struct Patient *patient);

void addList(std::atomic<Patient*>* list, Patient* patient);
void c_addList(std::atomic<Patient*>* list, Patient* patient);
void removeList(std::atomic<Patient*>* list, Patient* patient);

void check_patients_inside(struct Village *village);
void check_patients_waiting(struct Village *village);
void check_patients_realloc(struct Village *village);

void check_patients_assess_par(struct Village *village);

float get_num_people(struct Village *village);
float get_total_time(struct Village *village);
float get_total_hosps(struct Village *village);

struct Results get_results(struct Village *village);

void read_input_data(char *filename);
void allocate_village( struct Village **capital, struct Village *back, struct Village *next, int level, int32_t vid);

void sim_village_main_par(Village *top, size_t numthreads);

void sim_village_par(struct Village *village);
int check_village(struct Village *top);

void check_patients_assess(struct Village *village);
void check_patients_population(struct Village *village);
void sim_village(struct Village *village);
void my_print(struct Village *village);

#endif
