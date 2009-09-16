
#include "hashtable.h"
#include "hashtable_itr.h"
#include "hashtable_utility.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>		/* for memcmp */


typedef unsigned int uint32_t;
typedef unsigned short uint16_t;

struct key
{
  uint32_t ego_id;
  uint32_t altar_id;
  uint16_t function_id;
  uint16_t two_port;
};

struct value
{
  double rate;
  int function_id;
  int ego_id;
  int altar_id;
  int strain_id;
  int phylo_id;
};

DEFINE_HASHTABLE_INSERT (insert_some, struct key, struct value);
DEFINE_HASHTABLE_SEARCH (search_some, struct key, struct value);
DEFINE_HASHTABLE_REMOVE (remove_some, struct key, struct value);
DEFINE_HASHTABLE_ITERATOR_SEARCH (search_itr_some, struct key);

