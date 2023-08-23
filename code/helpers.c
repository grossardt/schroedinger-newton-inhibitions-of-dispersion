/* This code is part of Schroedinger Newton Inhibitions of Dispersion
   (C) Copyright Andre Grossardt 2010-2023
   https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
   
   This code is licensed under the MIT License (see LICENSE.txt for details)
 */

/* 
 * This file contains helper functions for file handling and output.
 */

/* Includes */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>


/* Parameters for the run are in param.h */
#include "param.h"


/* Build path name */
void outpath_name ( char opath[238], const char *path, const char timestr[14] )
{
    char first[222];
    size_t len;
    
    strncpy ( first, path, sizeof( first ) );
    len = strlen ( first );
    /* cut final / */
    if ( opath[len - 1] == '/' )
    {
        opath[len - 1] = '\0';
    }
    /* append timestring */
    sprintf ( opath, "%s/%s", first, timestr );
}

/* Create path for output */
void make_outpath ( char opath[238], const char *path )
{
    char *p;
    char timestr[14];
    char dpath[243];
    size_t len;
    time_t rawtime;
    int i = 0;
    
    /* get date and time */
    time ( &rawtime );
    strftime ( timestr, 14, "%Y%m%d-%H%M", localtime ( &rawtime ) );
    /* append date and time to path */
    outpath_name( opath, path, timestr );
    len = strlen ( opath );
    /* make sure that two instances of the program will not use the same path */
    while ( ! access ( opath, F_OK ) ) /* as long as directory exists */
    {
        /* append number to path */
        opath[len] = '\0';
        sprintf ( opath, "%s-%1d", opath, ++i );
        if ( i>9 )
        {
            printf ( "Error: Too many identical directory names. Quitting." );
            exit( 0 );
        }
    }
    /* subpath for data */
    sprintf ( dpath, "%s/data", opath );
    /* create subdirectories */
    for ( p = opath; *p; p++ )
    {
        if ( *p == '/' )
        {
            *p = '\0';
            if ( access ( opath, F_OK ) ) /* returns 0 if successfully accessed */
            {
                mkdir ( opath, S_IRWXU ); /* read, write, execute by user */
            }
            *p = '/';
        }
    }
    /* create final directory */
    if ( access ( opath, F_OK ) )
    {
        mkdir ( opath, S_IRWXU );
    }
    /* create data directory */
    if ( access ( dpath, F_OK ) )
    {
        mkdir ( dpath, S_IRWXU );
    }
}

/* Write the parameters to param.txt */
void save_settings( const char *path, unsigned long t )
{
    FILE *fp;
    char filename[256];
    char timestr[14];
    time_t rawtime;
    static const long double totaltime = 1e-9L * DT * MAXT;
    
    sprintf ( filename, "%s/param.txt" , path );
    time ( &rawtime );
    strftime ( timestr, 14, "%Y%m%d-%H%M", localtime ( &rawtime ) );
    /* Open file to append, for each continuation the settings will be saved */
    if ( ( fp = fopen ( filename, "a" ) ) == NULL )
    {
        printf ( "Cannot open file. Quitting.\n" );
        exit( 0 );   
    }
    else /* write data */
    {
        fprintf ( fp, "*********************************\n" );
        fprintf ( fp, " SETTINGS FOR RUN @ %s\n", timestr);
        fprintf ( fp, " STARTING WITH t = %14lu\n", t );
        fprintf ( fp, "*********************************\n" );
        fprintf ( fp, "width in nm: %20Lg\n", W );
        fprintf ( fp, "mass in u  : %20Lg\n", M );
        fprintf ( fp, "grid size  : %20d\n",  N );
        fprintf ( fp, "dr in nm   : %20Lg\n", DR );
        fprintf ( fp, "dt in ns   : %20Lg\n", DT );
        fprintf ( fp, "max. time  : %20lu\n", MAXT );
        fprintf ( fp, "save every : %20lu\n", SAVEEVERY );
        fprintf ( fp, "coupling   : %20Lg\n", COUP );
        fprintf ( fp, "wave funct.: %20c\n",  WAVEFUNCT );
        if ( WAVEFUNCT == 'g' )
            fprintf ( fp, "movie cmd  : python movie.py %s %d %lu %lu \"%5Lg u, %Lg s\" %Lg %Lg %Lg %Lg\n", path, N, MAXT, SAVEEVERY, M, totaltime, W, M, DR, DT );
        else
            fprintf ( fp, "movie cmd  : python movie.py %s %d %lu %lu \"%5Lg u, %Lg s\"\n", path, N, MAXT, SAVEEVERY, M, totaltime );
        #ifdef SALZMAN_ERROR
        fprintf ( fp, "!Calculated with Salzman's error!\n" );
        #endif
        #ifdef CHECK_OFF
        fprintf ( fp, "!results obtained in UNSAFE mode!\n" );
        #endif
        fprintf ( fp, "\n");
    }
    fclose ( fp );
}

/* Open wave function file */
static FILE *open_wf_file ( int write, unsigned long t, const char *path )
{
    FILE *fp;
    char filename[256];
    char rwmode[3];
    
    strcpy ( rwmode, write ? "wb" : "rb");
    sprintf ( filename, "%s/data/w%014lu.dat" , path, t );
    if ( ( fp = fopen ( filename, rwmode ) ) == NULL )
    {
        printf ( "Cannot open file. Quitting.\n" );
        exit( 0 );
    }
    return ( fp );
}

/* Save wave function */
void save_wf ( unsigned long t, long double complex psi[N], const char *path )
{
    FILE *fp;
    fp = open_wf_file ( 1, t, path );
    if ( fwrite ( psi, sizeof ( long double complex ), N, fp ) != N )
    {
        printf ( "File read error. Quitting.\n" );
        exit( 0 );
    }
    fclose ( fp );
}

/* Load wave function */
void load_wf ( unsigned long t, long double complex psi[N], const char *path )
{
    FILE *fp;
    fp = open_wf_file ( 0, t, path );

    if ( fread ( psi, sizeof ( long double complex ), N, fp ) != N )
    {
        if ( feof ( fp ) )
        {
            printf("Premature end of file. Quitting.\n");
        }
        else
        {
            printf("File read error. Quitting.\n");
        }
        exit( 0 );
    }
    fclose ( fp );
}

/* Show progress */
void progress ( unsigned long t )
{
    printf ( "Progress %lu%%", 100UL * t / MAXT );
    printf ( "\r" );
    fflush ( stdout );
}

/* Notify when continuing */
void cont_notify ( unsigned long t, const char *path )
{
    printf ( "Will continue writing to directory %s from t=%lu.\n", path, t );
    printf ( "All data beyond this time will be overwritten!\n" );
    printf ( "Press [Enter] to continue . . ." );
    fflush ( stdout );
    getchar();
}

