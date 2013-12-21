#ifndef IDLE_QUEUE_H
#define IDLE_QUEUE_H

typedef struct idle_element IDLE_ELEMENT;

struct idle_element {
    int data;
    IDLE_ELEMENT *link;
};

IDLE_ELEMENT* addslv (IDLE_ELEMENT * lp, int data);
/* Adds an item to the queue. */

IDLE_ELEMENT* removeslv (IDLE_ELEMENT * lp);
/* Returns the pointer to the next available item in the queue. */

int num_idle(IDLE_ELEMENT * lp);
/* Returns the number of idle processors in the queue. */

#endif
