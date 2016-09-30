/* $Header: /share/src/hpcc/tipsy/code/shell_sub.c,v 1.1.1.1 1995/01/10 22:57:37 trq Exp $
 * $Log: shell_sub.c,v $
 * Revision 1.1.1.1  1995/01/10 22:57:37  trq
 * Import to CVS
 *
 * Revision 1.1  94/04/19  17:56:03  trq
 * Initial revision
 * 
 */

#include "defs.h"
#include "fdefs.h"
#include <ctype.h>

void
shell_sub(job)
    char *job;
{
  char *command;
  
  command = job;
  while(!isspace((int) *command))
    command++;
  if(*command)
    {
      system(command);
    }
  else
    {
      input_error("shell");
    }
}
