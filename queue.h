/***** queue.h ************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Nov 18 17:15:33 2009
 **************************************************/

typedef struct queueItem{ /* the queue item: */
  int pos;                /* position */
  char **profile;         /* allele profiles at pos */
}QueueItem;

void queueInit(int maxN);
int queueEmpty();
void queuePut(QueueItem *item);
QueueItem *queueGet();
QueueItem *queuePeek();
QueueItem *newQueueItem();
void fillQueueItem(char *line, QueueItem *q);
void freeQueueItem(QueueItem *q);
void freeQueue();
