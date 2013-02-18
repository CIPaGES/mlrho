/***** queue.c ************************************
 * Description: FIFO queue array implementation
 *   taken from p. 157 of Sedgewick, R. (1998).
 *   Algorithms in C, Parts 1-4. Third Edition. 
 *   Addision-Wesley, Boston.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Nov 18 17:15:30 2009
 **************************************************/
#include <stdlib.h>
#include "queue.h"
#include "eprintf.h"

static QueueItem **q;
static int n, head, tail;

void queueInit(int maxN)
{
  q = (QueueItem **)emalloc((maxN+1)*sizeof(QueueItem *));
  n = maxN+1;
  head = n;
  tail = 0;
}

int queueEmpty()
{
  return head % n == tail;
}

void queuePut(QueueItem *item)
{
  q[tail++] = item;
  tail = tail % n;
}

QueueItem *queueGet()
{
  head = head % n;
  return q[head++];
}

QueueItem *queuePeek()
{
  head = head % n;
  return q[head];
}

QueueItem *newQueueItem()
{
  QueueItem *qi;
  int i;

  qi = (QueueItem *)emalloc(sizeof(QueueItem));
  qi->profile = (char **)emalloc(4*sizeof(char *));
  for(i=0;i<4;i++)
    qi->profile[i] = (char *)emalloc(10*sizeof(char));

  return qi;
}

void fillQueueItem(char *line, QueueItem *qi)
{
  sscanf(line,"%d\t%s\t%s\t%s\t%s",&qi->pos,qi->profile[0],
	  qi->profile[1], qi->profile[2], qi->profile[3]);
}

void freeQueueItem(QueueItem *qi)
{
  int i;

  for(i=0;i<4;i++)
    free(qi->profile[i]);
  free(qi->profile);
  free(qi);
}

void freeQueue()
{
    while(!queueEmpty())
	freeQueueItem(queueGet());
    free(q);
}
