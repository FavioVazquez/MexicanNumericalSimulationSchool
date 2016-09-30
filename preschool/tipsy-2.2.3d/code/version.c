#include "defs.h"

static char version_string[] = "2.2.3d";

void
version(job)
    char job[MAXCOMM] ;
{
    printf("version %s:  15 DEC 2012\n", version_string) ;
    printf("MAXBOX = %d\n", MAXBOX) ;
}
