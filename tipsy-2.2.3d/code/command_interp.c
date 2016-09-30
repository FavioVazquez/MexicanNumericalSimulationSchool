/* $Header: /share/src/hpcc/tipsy/code/command_interp.c,v 1.1.1.1 1995/01/10 22:57:24 trq Exp $
 * $Log: command_interp.c,v $
 * Revision 1.1.1.1  1995/01/10 22:57:24  trq
 * Import to CVS
 *
 * Revision 2.3  94/04/20  08:45:27  trq
 * Added title variable.
 * 
 */
#include "defs.h"
#include "comm.h"

void
command_interp(job)
    char job[MAXCOMM] ;
{
    char command[MAXCOMM] ;
    int i;

    if(sscanf(job,"%s",command) != 1) {
	return ;
    }
    for(i=0;i<ncomm;i++) {
        if(strcmp(command,c_list[i].name) == 0) {
	    if(c_list[i].display && !display) {
		printf("<sorry not in display mode, %s>\n",title) ;
	    }
	    else {
		(*(c_list[i].funk))(job) ;
	    }
	    return;
	}
    }
    printf("<sorry %s is not a command, %s>\n",command,title) ;
}
