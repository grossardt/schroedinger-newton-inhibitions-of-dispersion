/* This code is part of Schroedinger Newton Inhibitions of Dispersion
   (C) Copyright Andre Grossardt 2010-2023
   https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
   
   This code is licensed under the MIT License (see LICENSE.txt for details)
 */

/* Header file for helpers.c */

#ifndef MODULE_HELPERS_H
#define MODULE_HELPERS_H
#include "param.h"
void make_outpath ( char opath[256], const char *path );
void save_settings( const char *path, unsigned long t );
void save_wf ( unsigned long t, long double complex psi[N], const char *path );
void load_wf ( unsigned long t, long double complex psi[N], const char *path );
void progress ( unsigned long t );
void cont_notify ( unsigned long t, const char *path );
#endif /* MODULE_HELPERS_H */
