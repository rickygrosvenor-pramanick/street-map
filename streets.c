#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "streets.h"


struct node {
    double lat;
    double lon;
    int nodeID;
    int numWays;
    int *wayIDs;
};

struct way {
    char *name;
    int wayID;
    bool oneWay;
    float maxSpeed;
    int numNodes;
    int *nodeIDs;
};

struct ssmap {
    int numNodes;
    int numWays;
    struct node *allNodes;
    struct way *allWays;
};

struct ssmap *
ssmap_create(int nr_nodes, int nr_ways)
{
    // Assumption: Both parameters will be a non-negative number, but can be zero.
    if (nr_nodes == 0 || nr_ways == 0) {
        return NULL;
    }

    struct ssmap *map = malloc(sizeof(struct ssmap));
    if (map == NULL) {
      // CASE: Out of Memory
       return NULL;
    }

    map->numNodes = nr_nodes;
    map->numWays = nr_ways;
    map->allNodes = malloc(sizeof(struct node) * nr_nodes);
    map->allWays = malloc(sizeof(struct way) * nr_ways);

    // CASE: Out of Memory
    if (map->allNodes == NULL || map->allWays == NULL) {
        if (map->allNodes != NULL) {
            free(map->allNodes);
        }
        if (map->allWays != NULL) {
            free(map->allWays);
        }
        free(map);
        return NULL;
    }
    return map;
}

bool
ssmap_initialize(struct ssmap * m)
{
    return true;
}

void
ssmap_destroy(struct ssmap * m)
{
    for (int i = 0; i < m->numWays; i++) {
        free(m->allWays[i].name);
        free(m->allWays[i].nodeIDs);
    }
    free(m->allWays);

    for (int i = 0; i < m->numNodes; i++) {
        free(m->allNodes[i].wayIDs);
    }
    free(m->allNodes);
    free(m);
}

struct way *
ssmap_add_way(struct ssmap * m, int id, const char * name, float maxspeed, bool oneway,
              int num_nodes, const int node_ids[num_nodes])
{
    struct way newWay = m->allWays[id];
    newWay.name = strdup(name);

    if (newWay.name == NULL) {
        return NULL;
    }

    newWay.wayID = id;
    newWay.maxSpeed = maxspeed;
    newWay.oneWay = oneway;
    newWay.numNodes = num_nodes;
    newWay.nodeIDs = malloc(sizeof(int) * num_nodes);

    if (newWay.nodeIDs == NULL) {
        free(newWay.name);
        return NULL;
    }

    for (int i = 0; i < num_nodes; i++) {
        newWay.nodeIDs[i] = node_ids[i];
    }

    m->allWays[id] = newWay;
    return &(m->allWays[id]);
}

struct node *
ssmap_add_node(struct ssmap * m, int id, double lat, double lon,
               int num_ways, const int way_ids[num_ways])
{
    struct node newNode = m->allNodes[id];

    newNode.nodeID = id;
    newNode.lat = lat;
    newNode.lon = lon;
    newNode.numWays = num_ways;
    newNode.wayIDs = malloc(sizeof(int) * num_ways);

    if (newNode.wayIDs == NULL) {
        return NULL;
    }

    for (int i = 0; i < num_ways; i++) {
        newNode.wayIDs[i] = way_ids[i];
    }
    m->allNodes[id] = newNode;
    return &(m->allNodes[id]);
}

void
ssmap_print_way(const struct ssmap * m, int id)
{
    if (id < 0 || id >= m->numWays) {
        // as 0 <= id <= m->numWays - 1
        printf("error: way %d does not exist.\n", id);
        return;
    }

    const struct way *currWay = &m->allWays[id];

    printf("Way %d: %s\n", id, currWay->name);
}

void
ssmap_print_node(const struct ssmap * m, int id)
{
    if (id < 0 || id >= m->numNodes) {
        // as 0 <= id <= m->numNodes - 1
        printf("error: node %d does not exist.\n", id);
        return;
    }
    printf("Node %d: (%0.7f, %0.7f)\n", id, m->allNodes[id].lat, m->allNodes[id].lon);

}

void
ssmap_find_way_by_name(const struct ssmap * m, const char * name)
{
    for (int i = 0; i < m->numWays; i++) {
        const struct way *currWay = &m->allWays[i];
        if (strstr(currWay->name, name) != NULL) {
            printf("%d ", i);
        }
    }
    printf("\n");
}

void
ssmap_find_node_by_names(const struct ssmap * m, const char * name1, const char * name2)
{
    for (int i = 0; i < m->numNodes; i++) {
        int name1PriorityFlag = 0;
        const struct node currNode = m->allNodes[i];

        int name1Index = INVALID_ID;
        int name2Index = INVALID_ID;

      for (int j = 0; j < currNode.numWays; j++) {
          const struct way currWay = m->allWays[currNode.wayIDs[j]];
          if ((strstr(currWay.name, name1) != NULL) && (name1Index == -1) && (name2Index != currWay.wayID)) {
              name1Index = currWay.wayID;
          }
          if ((name2 != NULL) && (strstr(currWay.name, name2) != NULL) && (name2Index == INVALID_ID) && (name1Index != currWay.wayID)) {
              name2Index = currWay.wayID;
          }
      }
      if (name1Index != INVALID_ID && name2 == NULL) {
          name1PriorityFlag = 1;
      } else if (name1Index != INVALID_ID && name2Index != INVALID_ID) {
          name1PriorityFlag = 1;
      }

      if (name1PriorityFlag == 0) {
          int name1Index = INVALID_ID;
          int name2Index = INVALID_ID;

        for (int j = 0; j < currNode.numWays; j++) {
            const struct way currWay = m->allWays[currNode.wayIDs[j]];
            if ((name2 != NULL) && (strstr(currWay.name, name2) != NULL) && (name2Index == INVALID_ID) && (name1Index != currWay.wayID)) {
                name2Index = currWay.wayID;
            }
            if ((strstr(currWay.name, name1) != NULL) && (name1Index == INVALID_ID) && (name2Index != currWay.wayID)) {
                name1Index = currWay.wayID;
            }
        }

        if (name1Index != INVALID_ID && name2 == NULL) {
            name1PriorityFlag = 1;
        } else if (name1Index != INVALID_ID && name2Index != INVALID_ID) {
            name1PriorityFlag = 1;
        }

      }
      if (name1PriorityFlag == 1) {
          printf("%d ", currNode.nodeID);
      }
    }
    printf("\n");
}

/**
 * Converts from degree to radian
 *
 * @param deg The angle in degrees.
 * @return the equivalent value in radian
 */
#define d2r(deg) ((deg) * M_PI/180.)

/**
 * Calculates the distance between two nodes using the Haversine formula.
 *
 * @param x The first node.
 * @param y the second node.
 * @return the distance between two nodes, in kilometre.
 */
static double
distance_between_nodes(const struct node * x, const struct node * y) {
    double R = 6371.;
    double lat1 = x->lat;
    double lon1 = x->lon;
    double lat2 = y->lat;
    double lon2 = y->lon;
    double dlat = d2r(lat2-lat1);
    double dlon = d2r(lon2-lon1);
    double a = pow(sin(dlat/2), 2) + cos(d2r(lat1)) * cos(d2r(lat2)) * pow(sin(dlon/2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return R * c;
}

double
ssmap_path_travel_time(const struct ssmap * m, int size, int node_ids[size])
{
    double totalTravelTime = 0.0;

    // Each node id specified by the user must be valid.
    for (int i = 0; i < size; i++) {
        int currNodeID = node_ids[i];
        if (currNodeID < 0 || currNodeID >= m->numNodes) {
            printf("error: node %d does not exist.\n", currNodeID);
            return -1.0;
      }
    }

    // It should not be possible to see the same node twice in a reasonable path.
    bool visitedNodes[m->numNodes];
    for (int i = 0; i < m->numNodes; i++){
        visitedNodes[i] = false;
    }

    for (int i = 0; i < size; i++) {
        int currNodeID = node_ids[i];
        if (visitedNodes[currNodeID] == true) {
            printf("error: node %d appeared more than once.\n", currNodeID);
            return -1.0;
        }
        visitedNodes[currNodeID] = true;
    }

    // Calculating Travel time
    for (int i = 0; i < size - 1; i++) {
        int nodeID1 = node_ids[i];
        int nodeID2 = node_ids[i + 1];

        const struct node *node1 = &m->allNodes[nodeID1];
        const struct node *node2 = &m->allNodes[nodeID2];

      // There must be a way object between two adjacent nodes in the array.
        int wayID = INVALID_ID;

        for (int j = 0; j < node1->numWays; j++) {
            const struct way *currWay = &m->allWays[node1->wayIDs[j]];
            for (int k = 0; k < node2-> numWays; k++) {
                const struct way *currWayNode2 = &m->allWays[node2->wayIDs[k]];
                if (currWay->wayID == currWayNode2->wayID) {
                    wayID = currWay->wayID;
                    break;
                }
            }

        if (wayID != INVALID_ID) {
            break;
        }
        }

      if (wayID == INVALID_ID) {
          printf("error: there are no roads between node %d and node %d.\n", nodeID1, nodeID2);
          return -1.0;
      }
      // Two adjacent nodes in a path must also be adjacent in the way object they belong to.

      const struct way *wayWithBothNodes = &m->allWays[wayID];
      int indexNode1 = INVALID_ID;
      int indexNode2 = INVALID_ID;
      int oppositeFlag = INVALID_ID;

      for (int k = 0; k < wayWithBothNodes->numNodes - 1; k++) {
          if (wayWithBothNodes->nodeIDs[k] == nodeID1 && wayWithBothNodes->nodeIDs[k + 1] == nodeID2) {
              indexNode1 = k;
              indexNode2 = k + 1;
              break;
          } else if (wayWithBothNodes->nodeIDs[k] == nodeID2 && wayWithBothNodes->nodeIDs[k + 1] == nodeID1) {
              indexNode1 = k + 1;
              indexNode2 = k;
              oppositeFlag = indexNode2;
              break;
          }
      }

      if (indexNode1 == INVALID_ID || indexNode2 == INVALID_ID) {
          printf("error: cannot go directly from node %d to node %d.\n", nodeID1, nodeID2);
          return -1.0;
      }

      // If a way object is one-way, the adjacent nodes in the path sharing that way object must also be
      // in the same order.

      if (wayWithBothNodes->oneWay == true && oppositeFlag != INVALID_ID) {
          printf("error: cannot go in reverse from node %d to node %d.\n", nodeID1, nodeID2);
          return -1.0;
      }

      double distance = distance_between_nodes(node1, node2);
      double travel_time = distance / (m->allWays[wayID].maxSpeed);
      totalTravelTime += travel_time;
    }

    totalTravelTime *= 60;
    return totalTravelTime;
  }


struct adjacencyNode {
  int numNodes;
  int *allNodesAdjacent;
};
typedef struct adjacencyNode AdjacencyNode;

struct adjacencyList {
    // Index-based array where the ID of the node in the adjacencyList, corresponds its
    // adjacencyNode item.
   AdjacencyNode *adjacencyList;
};
typedef struct adjacencyList AdjacencyList;


// This function makes the linked list implementation redundant (no need to do struct adjacencyNode *next).
// This is because we now know how many adjacent nodes there are.
int adjacentNodeCount(const struct ssmap *m, const struct node *currNode) {
  int countSoFar = 0;

  for (int i = 0; i < currNode->numWays; i++) {
      struct way currWay = m->allWays[currNode->wayIDs[i]];
      for (int j = 0; j < currWay.numNodes; j++) {
          if (currNode->nodeID == currWay.nodeIDs[j] && currWay.oneWay == true) {
              if (j < currWay.numNodes - 1) {
                  countSoFar += 1;
              }
          } else if (currNode->nodeID == currWay.nodeIDs[j] && currWay.oneWay == false) {
              if (j < currWay.numNodes - 1 && j > 0) {
                  countSoFar += 2;
              } else if (j > 0) {
                  countSoFar += 1;
              } else if (j < currWay.numNodes - 1) {
                  countSoFar += 1;
              }

              }
          }
      }
  return countSoFar;
}

// general idea is that each node has a list of nodes that are adjacent to it
// then when we choose a path, we don't need to worry about if a node is in a way with
// other nodes as we already have it in the adjacency list.
AdjacencyList *initialiseAdjacencyList (const struct ssmap *m) {
    AdjacencyList *aList = malloc(sizeof(AdjacencyList));
    aList->adjacencyList = malloc((m->numNodes) * sizeof(AdjacencyNode));

    for (int i = 0; i < m->numNodes; i++) {
        struct node currNode = m->allNodes[i];
        AdjacencyNode adjacencyNode;
        adjacencyNode.numNodes = adjacentNodeCount(m, &currNode);
        adjacencyNode.allNodesAdjacent = malloc(sizeof(int) * adjacentNodeCount(m, &currNode));

        int insert_index = 0;
        for (int j = 0; j < currNode.numWays; j++) {
            struct way currWay = m->allWays[currNode.wayIDs[j]];
            for (int k = 0; k < currWay.numNodes; k++) {
                if (currWay.nodeIDs[k] == currNode.nodeID && currWay.oneWay == true) {
                    if (k < currWay.numNodes - 1) {
                        adjacencyNode.allNodesAdjacent[insert_index] = currWay.nodeIDs[k + 1];
                        insert_index += 1;
                    }
                } else if (currWay.nodeIDs[k] == currNode.nodeID && currWay.oneWay == false) {
                    if (k < currWay.numNodes - 1 && k > 0) {
                      adjacencyNode.allNodesAdjacent[insert_index] = currWay.nodeIDs[k + 1];
                      adjacencyNode.allNodesAdjacent[insert_index + 1] = currWay.nodeIDs[k - 1];
                      insert_index += 2;
                    } else if (k > 0) {
                        adjacencyNode.allNodesAdjacent[insert_index] = currWay.nodeIDs[k - 1];
                        insert_index += 1;
                    } else if (k < currWay.numNodes - 1) {
                        adjacencyNode.allNodesAdjacent[insert_index] = currWay.nodeIDs[k + 1];
                        insert_index += 1;
                    }
                }

            }
        }
        aList->adjacencyList[currNode.nodeID] = adjacencyNode;
    }
    return aList;
}

struct minHeapItem {
  int nodeID;
  double dist;
};
typedef struct minHeapItem MinHeapItem;

struct minHeap {
  MinHeapItem *nodes;
  int maxNumElements;
  int currNumElements;
};
typedef struct minHeap MinHeap;

MinHeap *addWithPriority(MinHeap *minHeap, int nodeID, double newDist) {
    MinHeapItem newHeapItem;
    newHeapItem.nodeID = nodeID;
    newHeapItem.dist = newDist;

    MinHeapItem *nodes = minHeap->nodes;
    int size = minHeap->currNumElements;

    nodes[size] = newHeapItem;
    minHeap->currNumElements += 1;

  // Bubble Up
    int i = size;
    while ((nodes[i].dist < nodes[(i - 1) / 2].dist) && i > 0) {
        MinHeapItem temp = nodes[(i - 1) / 2];
        // Swap Parent and Child
        nodes[(i - 1) / 2] = nodes[i];
        // Swap Child and Parent
        nodes[i] = temp;

        i = (i - 1) / 2;
    }

    minHeap->nodes = nodes;
    return minHeap;
}

MinHeap *minHeapInitialise(const struct ssmap *m) {
    MinHeap *minHeap = malloc(sizeof(MinHeap));
    minHeap->nodes = malloc(sizeof(MinHeapItem) * m->numNodes);
    minHeap->maxNumElements = m->numNodes;
    minHeap->currNumElements = 0;

    for (int i = 0; i < m->numNodes; i++) {
        MinHeapItem initHeapItem;
        initHeapItem.nodeID = INVALID_ID;
        initHeapItem.dist = 0.0;
        minHeap->nodes[i] = initHeapItem;
    }
    return minHeap;
}

MinHeap *bubbleDown(MinHeap *minHeap, int i) {
    int size = minHeap->currNumElements;

    // If the Left Child has distance less than the current minHeapItem,
    // which is stored in minHeap->nodes[i], then Right Child has higher priority,
    // as in a minHeap, priority is determined by having a lower distance.
    // If both, left and right child have a distance less than current minHeapItem's
    // distance, then the switch would take place with the left child.
    if ((2 * i + 1) < size && (minHeap->nodes[i].dist > minHeap->nodes[2 * i + 1].dist)) {
        MinHeapItem temp = minHeap->nodes[2 * i + 1];
        minHeap->nodes[2 * i + 1] = minHeap->nodes[i];
        minHeap->nodes[i] = temp;
        return bubbleDown(minHeap, 2 * i + 1);
    } else if ((2 * i + 2) < size && (minHeap->nodes[i].dist > minHeap->nodes[2 * i + 2].dist)) {
        // If the Right Child has distance less than the current minHeapItem,
        // which is stored in minHeap->nodes[i], then Right Child has higher priority,
        // as in a minHeap, priority is determined by having a lower distance.
        MinHeapItem temp = minHeap->nodes[2 * i + 2];
        minHeap->nodes[2 * i + 2] = minHeap->nodes[i];
        minHeap->nodes[i] = temp;
        return bubbleDown(minHeap, 2 * i + 2);
    } else {
        return minHeap;
    }
}

MinHeapItem extractMin(MinHeap *minHeap) {
  MinHeapItem head = minHeap->nodes[0];

  int size = minHeap->currNumElements;

  // Switch the last item in the heap with the head of the minHeap
  // and then bubble down.
  minHeap->nodes[0] = minHeap->nodes[size - 1];

  size = size - 1;

  MinHeapItem initHeapItem;
  initHeapItem.nodeID = INVALID_ID;
  initHeapItem.dist = 0.0;

  minHeap->nodes[size] = initHeapItem;

  minHeap->currNumElements = size;
  bubbleDown(minHeap, 0);
  return head;
}

/*
This function returns the distance between two nodes, given they are adjacent.
This is required for the Djikstra's Algorithm.
*/
double getEdge (const struct ssmap *m, int nodeID1, int nodeID2) {
    struct node node1 = m->allNodes[nodeID1];
    struct node node2 = m->allNodes[nodeID2];

    for (int i = 0; i < m->numWays; i ++) {
        struct way currWay = m->allWays[i];
    for (int j = 0; j < currWay.numNodes - 1; j++) {
        if (currWay.nodeIDs[j] == nodeID1 && currWay.nodeIDs[j + 1] == nodeID2) {
            double dist = distance_between_nodes(&node1, &node2);
            return dist / (currWay.maxSpeed);
      } else if ((currWay.nodeIDs[j] == nodeID2 && currWay.nodeIDs[j + 1] == nodeID1) && currWay.oneWay == false) {
            double dist = distance_between_nodes(&node1, &node2);
            return dist / (currWay.maxSpeed);
        }
    }
    }
  return 0.0;
}

void djikstraShortestPath(const struct ssmap * m, int start_id, int end_id) {
    MinHeap *minHeap = minHeapInitialise(m);
    int maxSize = minHeap->maxNumElements;
    int prev[maxSize];
    double dist[maxSize];

    dist[start_id] = 0.0;
    // As start_id is the start of our path, start_id does not have a previous.
    // Thus, we make the previous of the start_id an INVALID_ID.
    prev[start_id] = INVALID_ID;
    AdjacencyList *adjacencyList = initialiseAdjacencyList(m);
    for (int i = 0; i < m->numNodes; i++) {
        if (i != start_id) {
            dist[i] = INFINITY;
            prev[i] = INVALID_ID;
        }
    }

    addWithPriority(minHeap, start_id, dist[start_id]);
    while (minHeap->currNumElements > 0) {
        MinHeapItem extracted = extractMin(minHeap);
        int u = extracted.nodeID;

        if (u == end_id) {
          break;
        }
        for (int i = 0; i < adjacentNodeCount(m, &m->allNodes[u]); i++) {
            int v = adjacencyList->adjacencyList[u].allNodesAdjacent[i];
            double alt = dist[u] + getEdge(m, u, v);
            if (alt < dist[v]) {
                dist[v] = alt;
                prev[v] = u;
                /*
                Instead of filling the priority queue with all nodes in the initialization phase,
                it is also possible to initialize it to contain only source; then, inside the if
                alt < dist[v] block, the decrease_priority() becomes an add_with_priority()
                operation if the node is not already in the queue.
                */
                addWithPriority(minHeap, v, alt);
            }
        }
    }
    int S[minHeap->maxNumElements];
    int u = end_id;
    int counter = 0;
    if (prev[u] != INVALID_ID || u == start_id) {
        while (u != INVALID_ID) {
            S[counter] = u;
            u = prev[u];
            counter += 1;
        }
    }

    for (int i = counter - 1; i > -1; i -= 1) {
        printf("%d ", S[i]);
    }
    printf("\n");

    for (int i = 0; i < m->numNodes; i++) {
        free(adjacencyList->adjacencyList[i].allNodesAdjacent);
    }
    free(adjacencyList->adjacencyList);
    free(adjacencyList);
    free(minHeap->nodes);
    free(minHeap);
}

void
ssmap_path_create(const struct ssmap * m, int start_id, int end_id)
{
    if (start_id < 0 || start_id >= m->numNodes) {
        printf("error: node %d does not exist.\n", start_id);
    }

    if (end_id < 0 || end_id >= m->numNodes) {
        printf("error: node %d does not exist.\n", end_id);
    }

    djikstraShortestPath(m, start_id, end_id);
}
