/* This code is part of Schroedinger Newton Inhibitions of Dispersion
   (C) Copyright Andre Grossardt 2010-2023
   https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
   
   This code is licensed under the MIT License (see LICENSE.txt for details)
 */

/* 
 * This is the main file of the program.
 * Runs the simulation of the Schrödinger-Newton equation for the
 * parameters defined in param.h and outputs wave function files
 * to the specified path.
 */

/* Includes */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <fenv.h>

/* Parameters for the run are in param.h */
#include "param.h"
/* Wave function shapes */
#include "wf.h"
/* Helper functions */
#include "helpers.h"

/* Define some constants for relevant pre-factors
 * These numbers are in SI units whereas parameters are given in
 * atomic mass units (u), nanometers (nm), and nanoseconds (ns).
 */
/* pi * G / hbar * (1u in kg)^2 */
#define PIGOHBAR 5.4824699260461014e-30L
/* -hbar / 8 * 10^9 / (1u in kg) */
#define MHBAROEI -7.9384748449675167L

/* Activate floating point exceptions */
#pragma STDC FENV_ACCESS ON
/* test only for: FE_INVALID, FE_DIVBYZERO, FE_OVERFLOW, FE_UNDERFLOW*/
#define MY_FE_EXCEPT 29

/* prefactor: -i hbar/(8m) dt/(dr)^2 */
static const long double complex pre_beta = MHBAROEI * I / M / DR / DR * DT;
/* prefactors for diagonal elements b */
static const long double complex b_pre = 0.5L - 2.0L * MHBAROEI * I / M / DR / DR * DT;
static const long double complex bn_pre = 0.5L - 6.0L * MHBAROEI * I / M / DR / DR * DT;
/* prefactor for gravitation potential: V = (-) i * (pi G / hbar) dt (m dr)^2 * coupling constant */
static const long double complex v_pre = I * COUP * PIGOHBAR * DT * M * DR * M * DR;




/* Calculate the potential */
void grav_potential ( long double complex psi[N], long double complex b[N] )
{
    long double v[N];
    long double psisq;
    long double qi_sum;
    long double ldi;
    int i;
    
    /* Calculate v - v0 in first loop making use of
       v_j = v_0 + 1/j sum psi^2 i^2 - sum psi^2 i
     */
    psisq = cabsl ( psi[1] );
    psisq *= psisq;
    qi_sum = 0.5L * psisq;
    v[0] = psisq;
    v[1] = 0.0L;
    for ( i=2; i<=(N-2); ++i )
    {
        ldi = (long double) i;
        psisq = cabsl ( psi[i] );
        psisq *= psisq;
        psisq *= ldi;
        v[i] = qi_sum - v[0];
        v[0] += psisq;
        ldi /= ( ldi + 1.0L );
        qi_sum += psisq;
        qi_sum *= ldi;
    }
    psisq = cabsl ( psi[N-1] );
    psisq *= psisq;
    v[N-1] = qi_sum - v[0];
    v[0] += psisq * (long double) (N-1);
    
    /* Calculate diagonal elements b, the real v_i is v_i + v_0 */
    b[0] = bn_pre - v_pre * v[0];
    for ( i=1; i<N; ++i )
    {
        b[i] = b_pre - v_pre * ( v[0] + v[i] );
    }
    /* Check if v is big enough to make a difference
     * i.e. make sure that b_0 is not equal to bn_pre.
     * As v_0 > v_j for all j>0 considering b_0 is enough.
     * Due to efficiency, do this only if CHECK_OFF is not set.
     */
    #ifndef CHECK_OFF
    static int check_potential = 1;
    if ( check_potential && b[0] == bn_pre )
    {
        printf ( "Potential too weak to be represented numerically!\n" );
        check_potential = 0;
    }
    #endif /* CHECK_OFF */
}

/* Calculate inital Q matrix off-diagonals */
void q_init ( long double complex a[N], long double complex c[N] )
{
    long double iinv;
    int i;
    
    /* a (subdiagonal), b (diagonal), c (superdiagonal) */
    c[0] = 6.0L * pre_beta;
    for ( i=1; i<N; ++i )
    {
        iinv = 1.0L / (long double) i;
        a[i] = pre_beta * (1.0L - iinv);
        c[i] = pre_beta * (1.0L + iinv);
    }
}

/* Solve the system of linear equations with a tridiagonal matrix
 * arguments: a (subdiagonal), b (diagonal), c (superdiagonal)
 * cf. Background section in README
 */
void solve_linear_system ( long double complex a[N], long double complex b[N], long double complex c[N], long double complex psi[N] )
{
    long double complex bb[N];
    long double complex d[N];
    long double complex x[N];
    int i;
    
    /* Transform the matrix */
    bb[0] = b[0];
    d[0] = psi[0];
    for ( i=1; i<N; ++i )
    {
        bb[i] = b[i] - a[i] * c[i-1] / bb[i-1];
        d[i] = psi[i] - a[i] * d[i-1] / bb[i-1];
    }
    /* Form the solution x and psi_new = x - psi */
    x[N-1] = d[N-1] / bb[N-1];
    psi[N-1] = x[N-1] - psi[N-1];
    for ( i=N-2; i>=0; --i )
    {
        x[i] = ( d[i] - c[i]*x[i+1] ) / bb[i];
        psi[i] = x[i] - psi[i];
    }
}


/* Main routine */
int main ( int argc, char *argv[] )
{
    long double complex psi[N];
    long double complex a[N];
    long double complex b[N];
    long double complex c[N];
    unsigned long t = 0;
    char path[238];
    int cont = 0;
    
    /* Should we continue a former calculation? */
    if ( argc > 1 )
    {
        /* only c and s are allowed as first argument and continue mode needs three arguments */
        if ( ( argv[1][0] != 'c' && argv[1][0] != 's' ) || ( argv[1][0] == 'c' && argc < 4 ) )
        {
            printf ( "Error: Wrong number or type of arguments.\n" );
            printf ( "First argument must be either 's' (start new calculation) or 'c' (continue), " );
            printf ( "second argument must be the subpath (e.g. '20100101-0000') and " );
            printf ( "third argument must be the time step for continue mode.\nQuitting...\n" );
            exit( 0 );
        }
        else if ( argv[1][0] == 'c' ) /* continue calculation */
        {
            cont = 1;
            /* rebuild name of output path */
            outpath_name ( path, OUTDIR, argv[2] );
            /* read time where to continue */
            t = ( unsigned long ) atol ( argv[3] );
            /* Load wave function from saved file */
            load_wf ( t, psi, path );
            /* Notify and ask if we should continue */
            cont_notify ( t, path );
        }
    }
    
    if ( ! cont )
    {
        /* Create output directory */
        make_outpath ( path, OUTDIR );
        /* Initialise wave function */
        wave_function ( psi );
    }
    
    /* save the parameters of this run (append if continue) */
    save_settings ( path, t );
    
    /* Initialise Q matrix */
    q_init ( a, c );
    
    /* iterate wave function */
    while ( t++ < MAXT )
    {
        grav_potential ( psi, b );
        solve_linear_system ( a, b, c, psi );
        if ( ! ( int ) ( t % SAVEEVERY ) )
        {
            /* Save this step */
            save_wf ( t, psi, path );
            /* Print progress */
            progress ( t );
        }
        /* Check for exceptions */
        #ifndef CHECK_OFF
        if ( fetestexcept( MY_FE_EXCEPT ) )
        {
            printf ( "WARNING: A floating point exception occured for t=%lu: ", t );
            if ( fetestexcept( FE_INVALID ) ) printf( "INVALID " );
            if ( fetestexcept( FE_DIVBYZERO ) ) printf( "DIVBYZERO " );
            if ( fetestexcept( FE_OVERFLOW ) ) printf( "OVERFLOW " );
            if ( fetestexcept( FE_UNDERFLOW ) ) printf( "UNDERFLOW " );
            printf ( "\n" );
            feclearexcept( FE_ALL_EXCEPT );
        }
        #endif /* CHECK_OFF */
    }
    
    /* Save final wavefunction if not already done */
    if ( ( int ) ( --t % SAVEEVERY ) )
    {
        save_wf ( t, psi, path );
    }
    
    /* finished */
    progress ( MAXT );
    printf ( "\nDone.\n" );
    return ( 0 );
}

