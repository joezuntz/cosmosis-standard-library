#ifndef TICKER_H
#define TICKER_H

#include "globaldef.h"

#include <sys/time.h>
#include <iostream>

/** Estimates remaining time for calculation
 * Original code by Patrick Simon
 */
class ticker
{
 public:

  ticker();
  ~ticker();

  /// sets counter to zero
  void reset();

  /// prints the remaining time, to be called with fraction of work already done
  void tick(number);

  /// returns the overall runtime in minutes
  number runtime();

 private:

  int             started;
};

#endif
