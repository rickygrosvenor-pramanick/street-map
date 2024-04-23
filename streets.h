#ifndef _STREETS_H_
#define _STREETS_H_

/**
 * Valid node or way IDs start from 0, so we use -1 to denote invalid ID.
 */
#define INVALID_ID (-1)

struct ssmap;
struct node;
struct way;
struct path;

/**
 * Create a new ssmap data structure.
 *
 * @param nr_nodes Maximum number of nodes that will be added.
 * @param nr_ways Maximum number of ways that will be added.
 * @return A heap-allocated ssmap structure or NULL if malloc fails.
 */
struct ssmap * ssmap_create(int nr_nodes, int nr_ways);

/**
 * Perform any other initialization after ways and nodes have been added.
 *
 * @param m The ssmap data structure to be initialized.
 * @return true upon successful initialization, false otherwise.
 */
bool ssmap_initialize(struct ssmap * m);

/**
 * Destroy an existing ssmap data structure.
 *
 * @param m The ssmap structure to destroy. It should free ALL dynamically
 * allocated memory associated with the structure, including the nodes and
 * ways that have been added.
 */
void ssmap_destroy(struct ssmap * m);

/**
 * Add a new way object to the ssmap data structure.
 *
 * @param m The ssmap structure for which the way object should be added.
 * @param id The id of the way object.
 * @param name The name of the way object. Make sure you make a full copy 
 *             of it instead of just having a pointer to it.
 * @param maxspeed The speed limit of the street, in km/hr.
 * @param oneway Whether the street is one way.
 * @param num_nodes The number of nodes associated with this way object.
 * @param node_ids An array of node ids associated with this way object.
 * @return The way object that was just added.
 */
struct way * ssmap_add_way(struct ssmap * m, int id, const char * name, 
                           float maxspeed, bool oneway, int num_nodes, 
                           const int node_ids[num_nodes]);

/**
 * Add a new node object to the ssmap data structure.
 *
 * @param m The ssmap structure for which the node object should be added.
 * @param id The id of the node object.
 * @param lat The latitude of this node.
 * @param lon The longitude of this node.
 * @param num_ways The number of ways associated with this node object.
 * @param way_ids An array of way ids assocaited with this node object.
 * @return The node object that was just added. 
 */
struct node * ssmap_add_node(struct ssmap * m, int id, double lat, double lon, 
                             int num_ways, const int way_ids[num_ways]);

/**
 * Find a way object by id, then print its information
 * 
 * The first line of output should have the format:
 * Way <id>: <name>
 *
 * You have the freedom to add any additional output as you see fit.
 *
 * You must check that the way id is valid. If not, print "error: way
 * <id> does not exist."
 * @param m The ssmap structure where the way object should be located.
 * @param id The id of the way object to be printed.
 */
void ssmap_print_way(const struct ssmap * m, int id);

/**
 * Find a node object by id, then print its information
 * 
 * The first line of output should have the format:
 * Node <id>: (<lat>, <lon>)
 *
 * You have the freedom to add any additional output as you see fit.
 *
 * You must check that the node id is valid. If not, print "error: node
 * <id> does not exist."
 * @param m The ssmap structure where the node object should be located.
 * @param id The id of the node object to be printed.
 */
void ssmap_print_node(const struct ssmap * m, int id);

/**
 * Find all way objects with a particular keyword in its name and print them.
 *
 * The format of the output should be space separated way ids. e.g.
 * 9 32 105 5 782 96
 *
 * You may assume name will never be a null pointer.
 * 
 * @param m The ssmap structure where the way objects should be located.
 * @param name The keyword that should be used to search for way objects.
 */
void ssmap_find_way_by_name(const struct ssmap * m, const char * name);

/**
 * Find all node objects that are associated with way objects that have the 
 * specified keywords in their names. The main purpose of this function
 * is to find a particular intersection, e.g., Yonge and Bloor. A node should
 * only be printed if each keyword matches a distinct way object, e.g., Yonge
 * matches way 42 and Bloor matches way 69. 
 * 
 * The format of the output should be space separated node ids. e.g.
 * 9 32 105 5 782 96
 *
 * You may assume name1 will never be a null pointer, but not name2.
 * 
 * @param m The ssmap structure where the node objects should be located.
 * @param name1 The first keyword that should be used to search for a node.
 * @param name2 The second keyword that should be used to search for a node.
 *              This parameter is allowed to NULL. In which case, the function
 *              displays all nodes with ways that have name1 in their names.
 */
void ssmap_find_node_by_names(const struct ssmap * m, const char * name1, const char * name2);

/**
 * Calculate the travel time of a path (an ordered array of node ids)
 *
 * You may make the following simplifications
 * 1. There are no traffic. You are the only one driving on the road.
 * 2. You are always able to drive at the speed limit of the road you are on.
 * 3. There are no acceleration or deceleration. You instantly drive at the speed
 *    limit of the road, and you can turn instantaneously.
 * 4. There are no turn restrictions, i.e., you can turn in any direction at any 
 *    intersection. 
 *
 * You must check for the following errors:
 * 1. Each node id must be valid. If not, print "error: node <id> does not exist."
 * 2. There must be a way object between two adjacent nodes in the array. If not,
 *    print "error: there are no roads between node <a> and node <b>."
 * 3. If a way object contains multiple nodes [a, b, c, d, ...], you must not
 *    skip any nodes in between. If violated, print "error: cannot go directly from
 *     node <a> to node <b>."
 * 4. If a way object is one-way, the nodes in the path cannot go in reverse order.
 *    e.g., suppose a way object has nodes [a, b, c, d, ...] and it is one-way,
 *    then it would be a violation to go from node c to node b. If this happens,
 *    print "error: cannot go in reverse from node <a> to node <b>."
 *
 * If an error occurs, print the error message and return -1.0.
 * 
 * @param m The ssmap structure where the path lies.
 * @param size The size of the path array
 * @param node_ids The path array. Each element of the array is a node id, and
 * the array represents the order in which a traveler goes from one node to
 * another. E.g., suppose the path is [a, b, c, d], it means the travelers starts
 * from node a, goes to node b, c, and finaly node d in that order. For each 
 * adjacent set of nodes, there must be a valid way object that connects them.
 * @return the total travel time of a path in MINUTES, which can be calculated as 
 * the sum of the distance between each node divided by the speed limit of each 
 * segment. E.g., suppose the path is [a, b, c] and way x and y connects a and b,
 * b and c, respectively. Then the travel time would be distance(a, b) / maxspeed(x)
 * + distance(b, c) / maxspeed(y). Hint: use the distance function provided.
 */
double ssmap_path_travel_time(const struct ssmap * m, int size, int node_ids[size]);

/**
 * Compute a path from one node to another.
 *
 * You may make all of the assumptions in the ssmap_path_travel_time, and you
 * must ensure the absence of all errors described under ssmap_path_travel_time.
 *
 * The format of the output should be space-separated node ids. e.g.
 * 9 32 105 5 782 96
 *
 * Note: you should try to optimize this function so that the travel time is 
 * minimized.
 * 
 * @param m The ssmap structure where the path will be created.
 * @param start_id the starting node id 
 * @param end_id the destination node id
 */
void ssmap_path_create(const struct ssmap * m, int start_id, int end_id);

#endif /* _STREETS_H_ */