#ifndef LU_LU_TIMER_H
#define LU_LU_TIMER_H 1

#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>
double ATL_walltime(void)
{
   struct tms ts;
   static double ClockTick=0.0;

   if (ClockTick == 0.0) ClockTick = 1.0 / ((double) sysconf(_SC_CLK_TCK));
   return( ((double) times(&ts)) * ClockTick);
}

#endif // LU_LU_TIMER_H
