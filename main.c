#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "streets.h"

// use for reading from file and stdin
#define BUFSIZE 32768
char buffer[BUFSIZE];

#define RET_OK(expr, expected, label) do { \
    if ((expr) != (expected)) goto label; \
} while(0)

static bool
load_int_array(int size, int arr[size], FILE * f)
{
    for (int i = 0; i < size; i++) {
        if (fscanf(f, "%d", arr+i) != 1) {
            return false;
        }
    }
    fscanf(f, "\n");
    return true;
}

static void
remove_newline(char * string)
{
    char * newline = strchr(string, '\n');
    if (newline) {
        *newline = '\0';
    }
}

static struct ssmap * 
load_map(const char * filename)
{
    FILE * f = fopen(filename, "rt");
    struct ssmap * map = NULL;
    int nr_nodes, nr_ways;

    if (f == NULL) {
        fprintf(stderr, "error: could not open %s\n", filename);
        return NULL;
    }

    if (fgets(buffer, BUFSIZE, f) == NULL) {
        goto done;
    }

    RET_OK(strcmp(buffer, "Simple Street Map\n"), 0, invalid);
    RET_OK(fscanf(f, "%d ways\n", &nr_ways), 1, invalid);
    RET_OK(fscanf(f, "%d nodes\n", &nr_nodes), 1, invalid);

    map = ssmap_create(nr_nodes, nr_ways);
    if (map == NULL) {
        fprintf(stderr, "error: could not create ssmap\n");
        goto done;
    }

    for (int i = 0; i < nr_ways; i++) {
        int id, num_nodes;
        float maxspeed;
        char which_way[8];

        /* note: we are intentionally not loading the OSM id */
        RET_OK(fscanf(f, "way %d %*d ", &id), 1, cleanup);
        RET_OK(fgets(buffer, BUFSIZE, f), buffer, cleanup);
        RET_OK(fscanf(f, " %f %7s %d\n", &maxspeed, which_way, &num_nodes), 3, cleanup);

        remove_newline(buffer);
        bool oneway = strcmp(which_way, "oneway") == 0;
        struct way * way = NULL;

        if (num_nodes > 0) {
            int node_ids[num_nodes];
            RET_OK(load_int_array(num_nodes, node_ids, f), true, cleanup);
            way = ssmap_add_way(map, id, buffer, maxspeed, oneway, num_nodes, node_ids);
        }

        if (way == NULL) {
            goto cleanup;
        }
    }

    for (int i = 0; i < nr_nodes; i++) {
        int id, num_ways;
        double lat, lon;

        /* note: we are intentionally not loading the OSM id */
        RET_OK(fscanf(f, "node %d %*d %lf %lf %d\n", &id, &lat, &lon, &num_ways), 4, cleanup);
        struct node * node = NULL;
        
        if (num_ways > 0) {
            int way_ids[num_ways];
            RET_OK(load_int_array(num_ways, way_ids, f), true, cleanup);
            node = ssmap_add_node(map, id, lat, lon, num_ways, way_ids);
        }
         
        if (node == NULL) {
            goto cleanup;
        }
    }

    // custom initialization after all nodes and ways have been added
    if (!ssmap_initialize(map)) {
        goto cleanup;
    }

    printf("%s successfully loaded. %d nodes, %d ways.\n", filename, nr_nodes, nr_ways);
    goto done;
cleanup:
    ssmap_destroy(map);
invalid:
    fprintf(stderr, "error: %s has invalid file format\n", filename);
done:
    fclose(f);
    return map;
}

static bool
get_integer_argument(char * line, int * iptr)
{
    remove_newline(line);

    if (sscanf(line, "%d", iptr) == 1) {
        return true;
    }

    printf("error: '%s' is not an integer.\n", line);
    return false;
}

static void
handle_find(char * line, struct ssmap * map)
{
    char * command = strtok_r(line, " \t\r\n\v\f", &line);
    char * first = strtok_r(line, " \t\r\n\v\f", &line);
    char * second = strtok_r(line, " \t\r\n\v\f", &line);
    char * third = strtok_r(line, " \t\r\n\v\f", &line);

    if (command == NULL) {
        /* fall through */
    }
    else if (strcmp(command, "node") == 0) {
        if (first == NULL || third != NULL) {
            printf("error: invalid number of arguments.\n");
        }
        else {
            ssmap_find_node_by_names(map, first, second);
            return;
        }
    }
    else if (strcmp(command, "way") == 0) {
        if (first == NULL || second != NULL) {
            printf("error: invalid number of arguments.\n");
        }
        else {
            ssmap_find_way_by_name(map, first);
            return;
        }    
    }
    else {
        printf("error: first argument must be either node or way.\n");
    }

    printf("usage: find way keyword | find node keyword [keyword]\n");
}

static bool
handle_path_travel_time(char * line, struct ssmap * map)
{
    int capacity = 1;
    int n = 0;

    // use number of space characters to determine approximate array size
    for (int i = 0; line[i] != '\0'; i++) {
        if (isspace((int)line[i])) {
            capacity++;
        }
    }

    int node_ids[capacity];
    while(true) {
        char * token = strtok_r(line, " \t\r\n\v\f", &line);
        char * endptr;

        if (token == NULL)
            break;

        node_ids[n++] = strtol(token, &endptr, 10);
        if (endptr && *endptr != '\0') {
            printf("error: %s is not an integer.\n", token);
            return false;
        }
    }

    if (n < 2) {
        printf("error: must specify at least two nodes.\n");
        return false;
    }

    double result = ssmap_path_travel_time(map, n, node_ids);
    if (result >= 0.) {
        printf("%.4f minutes\n", result);
    }
    
    return true;
}

static bool
handle_path_create(char * line, struct ssmap * map)
{
    char * start = strtok_r(line, " \t\r\n\v\f", &line);
    char * finish = strtok_r(line, " \t\r\n\v\f", &line);
    char * endptr;

    if (start == NULL || finish == NULL) {
        printf("error: must specify start node and finish node.\n");
        return false;
    }

    int start_id = strtol(start, &endptr, 10);
    if (endptr && *endptr != '\0') {
        printf("error: %s is not an integer.\n", start);
        return false;
    }

    int end_id = strtol(finish, &endptr, 10);
    if (endptr && *endptr != '\0') {
        printf("error: %s is not an integer.\n", finish);
        return false;
    }

    ssmap_path_create(map, start_id, end_id);
    return true;
}

static void
handle_path(char * line, struct ssmap * map)
{
    char * command = strtok_r(line, " \t\r\n\v\f", &line);

    if (command == NULL) {
        /* fall through */
    }
    else if (strcmp(command, "time") == 0) {
        if (handle_path_travel_time(line, map))
            return;
    }
    else if (strcmp(command, "create") == 0) {
        if (handle_path_create(line, map))
            return;
    }
    else {
        printf("error: first argument must be either time or create.\n");
    }

    printf("usage: path create start finish | path time node1 node2 [nodes...]\n");
}

int 
main(int argc, const char * argv[])
{
    if (argc != 2) {
        fprintf(stderr, "usage: %s FILE\n", argv[0]);
        return 0;
    }

    struct ssmap * map = load_map(argv[1]);
    if (map == NULL) {     
        return 1;
    }

    while(true) {
        printf(">> ");
        fflush(stdout);
        char * ptr = fgets(buffer, BUFSIZE, stdin);
        if (ptr == NULL) {
            break;
        }

        char * command = strtok_r(buffer, " \t\r\n\v\f", &ptr);

        if (command == NULL) {
            /* fall through */
        }
        else if (strcmp(command, "quit") == 0) {
            break;
        }
        else if (strcmp(command, "node") == 0) {
            int id;
            if (get_integer_argument(ptr, &id)) {
                ssmap_print_node(map, id);
            }
        }
        else if (strcmp(command, "way") == 0) {
            int id;
            if (get_integer_argument(ptr, &id)) {
                ssmap_print_way(map, id);
            }
        }
        else if (strcmp(command, "find") == 0) {
            handle_find(ptr, map);
        }
        else if (strcmp(command, "path") == 0) {
            handle_path(ptr, map);
        }
        else {
            printf("error: unknown command %s. Available commands are:\n"
                   "\tnode, way, find, path, quit\n", command);
        }
    }
    
    ssmap_destroy(map);
    return 0;
}