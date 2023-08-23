/* This code is part of Schroedinger Newton Inhibitions of Dispersion
   (C) Copyright Andre Grossardt 2010-2023
   https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
   
   This code is licensed under the MIT License (see LICENSE.txt for details)
 */

/* 
 * This file contains the definition of wave function shapes.
 */

/* Includes */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#define M_2_SQRTPIl 1.1283791670955125738961589031215452L /* 2/sqrt(pi) */

/* Parameters for the run are in param.h */
#include "param.h"

/* All wave functions are in nm^(-3/2) */

/* Gaussian wave package */
static void gaussian_wf ( long double complex psi[N] )
{
    /* (pi w^2)^(-3/4) * exp(-r^2 / (2 w^2)) 
     * Normalisation: int |psi(x)|^2 d3x = 4 pi int |psi(r)|^2 r^2 dr = 1 / nm^3
     */
    static long double prefact;
    static const long double exp_pre = -.5L * DR / W * DR / W;
    int i;
    
    prefact = sqrtl ( M_2_SQRTPIl / W * M_2_SQRTPIl / W * M_2_SQRTPIl / W / 8.0L );
    for ( i=0; i<N; ++i )
    {
        psi[i] = ( long double complex ) ( prefact * expl( ( ( long double ) ( i*i ) ) * exp_pre ) );
    }
}

/* Rect wave package */
static void rect_wf ( long double complex psi[N] )
{
    /* Theta function, psi 0 to rect_index are rect_value */
    static const int rect_index = ( int ) ( W / DR );
    static long double complex rect_value;
    rect_value = 0.48860251190291992 / W / sqrtl ( W );
    int i = 0;
    while ( i<=rect_index )
    {
        psi[i++] = rect_value;
    }
    while ( i<N )
    {
        psi[i++] = 0.0L;
    }
}

/* Exponential wave package with hole in the middle
   psi(r) = 1/sqrt(9375 pi) w^(-5/2) r exp(-5r/w)
   width is just roughly estimated!
 */
static void exp_ball_wf ( long double complex psi[N] )
{
    static long double prefact;
    static const long double exp_pre = -5.0L * DR / W;
    int i;
    
    prefact = 5.82692496315775504289805709761433814768e-3L / sqrtl ( W * W * W * W * W );
    for ( i=0; i<N; ++i )
    {
        psi[i] = ( long double complex ) ( prefact * i * expl( ( ( long double ) i ) * exp_pre ) );
    }
}



/* Current wave function */
void wave_function ( long double complex psi[N] )
{
    switch ( WAVEFUNCT )
    {
        case 'r':
            rect_wf( psi );
            break;
        case 'b':
            exp_ball_wf( psi );
            break;
        default:
            gaussian_wf( psi );
    }
}
