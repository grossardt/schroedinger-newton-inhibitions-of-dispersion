/* This code is part of Schroedinger Newton Inhibitions of Dispersion
   (C) Copyright Andre Grossardt 2010-2023
   https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
   
   This code is licensed under the MIT License (see LICENSE.txt for details)
 */

/* Parameter file
   All parameters are defined here
 */

#define W           500.0L          // width in nm (long double)
#define M           50.0e9L         // mass in u (long double)
#define N           5100            // grid size (int)
#define DR          0.6L            // dr in nm (long double)
#define DT          1000000.0L      // dt in ns (long double)
#define MAXT        20000000UL      // number of time steps (unsigned long)
#define COUP        1.0L            // coupling constant for potential (long double)
#define SAVEEVERY   1000UL          // save every X time steps (long double)
#define OUTDIR      "/tmp/test"     // directory for output (string)
#define WAVEFUNCT   'g'             // initial wave function type (char):
                                    //   'g' = gaussian (default)
                                    //   'r' = rectangular
                                    //   'e' = exponential with hole in the middle

/* Note: Total time in ns should be about > 1000 * M (in u) */