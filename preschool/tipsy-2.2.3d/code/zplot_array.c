/*
 * $Header: /share/src/hpcc/tipsy/code/zplot_array.c,v 1.1.1.1 1995/01/10 22:57:36 trq Exp $
 * $Log: zplot_array.c,v $
 * Revision 1.1.1.1  1995/01/10 22:57:36  trq
 * Import to CVS
 *
 * Revision 1.3  94/04/20  08:46:38  trq
 * Added title variable.
 * 
 * Revision 1.1  94/02/16  13:53:25  trq
 * Initial revision
 * 
 */
#include "defs.h"
#include "fdefs.h"

void
zplot_array(job)
    char *job;
{
    if (binary_loaded) {
	if(active_box == 0 && !boxes_loaded[0]) {
	    loadall() ;
	}
	clear_rot() ;
	view_array(job) ;
    }
    else {
	printf("<sorry, no binary file is loaded, %s>\n",title) ;
    }
}
