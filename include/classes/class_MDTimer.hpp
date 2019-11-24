/*! \file md_timer
 * \brief Timer to well, time the simulation.
 */

#pragma once

#include <time.h>
#include <sys/time.h>

struct MDTimer {
    /*! \struct MDTimer
     * \brief Timer for simulation
     */

    struct timeval clockHolder; //!<
    struct timeval duration;
};

void initialize_timer (struct MD_Timer* time);
void start_timer (struct MD_Timer* time);
void stop_timer (struct MD_Timer* time);
double timer_duration(const struct MD_Timer& time);
