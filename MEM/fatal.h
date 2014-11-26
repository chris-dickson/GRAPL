/*--------------------------------------------------------------------------------------*/
/* fatal.h																				*/
/*	Used for error handling																*/
/*--------------------------------------------------------------------------------------*/
#ifndef FATAL_H_
#define FATAL_H_

#include <stdio.h>
#include <stdlib.h>

#define Error( Str )        FatalError( Str )
#define FatalError( Str )   fprintf( stderr, "%s\n", Str ), exit( 1 )

#endif
