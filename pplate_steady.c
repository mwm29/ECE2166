#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
//#include <mpi.h>

//Point is an xy pair
struct Point {
   double  x;
   double  y;
};
//Edge is 2 indicies to Points and a type
struct Edge {
   int  point1;
   int  point2;
   char type;// (s)ource, si(n)k, (w)all, or (i)nner
};
//Cell is 4 indicies to Edges
struct Cell {
   int  edge1;
   int  edge2;
   int  edge3;
   int  edge4;
};

void print_usage(char *filename){
    printf("\nIncorrect number of args:\n");
    printf("usage:\n");
    printf("\t%s <input_file> <velocity_in> <pressure_in> <pressure_out>\n\n", filename);
}

void validate_inputs(int argc, char **argv){
    if(argc < 5){
        print_usage(argv[0]);
        exit(1);
    }
}

void add_point(int id, double x, double y, struct Point **point_list, int *num_points, int *max_points){
    if(*num_points == *max_points){
        *max_points *= 2;
        (*point_list) = realloc((*point_list), (*max_points)*sizeof(struct Point));
    }
    (*point_list)[*num_points].x = x;
    (*point_list)[*num_points].y = y;
    *num_points += 1;
}

void add_edge(int id, char type, int index1, int index2, struct Edge **edge_list, int *num_edges, int *max_edges){
    if(*num_edges == *max_edges){
        *max_edges *= 2;
        (*edge_list) = realloc((*edge_list), (*max_edges)*sizeof(struct Edge));
    }
    (*edge_list)[*num_edges].point1 = index1;
    (*edge_list)[*num_edges].point2 = index2;
    (*edge_list)[*num_edges].type = type;
    *num_edges += 1;
}

void add_cell(int id, int index1, int index2, int index3, int index4, struct Cell **cell_list, int *num_cells, int *max_cells){
    if(*num_cells == *max_cells){
        *max_cells *= 2;
        (*cell_list) = realloc((*cell_list), (*max_cells)*sizeof(struct Cell));
    }
    (*cell_list)[*num_cells].edge1 = index1;
    (*cell_list)[*num_cells].edge2 = index2;
    (*cell_list)[*num_cells].edge3 = index3;
    (*cell_list)[*num_cells].edge4 = index4;
    *num_cells += 1;
}

void read_inputs(const char *filename, struct Cell **cell_list, int *num_cells, struct Edge **edge_list, int *num_edges, struct Point **point_list, int *num_points){
    FILE *infile;
    infile = fopen(filename, "r");

    *num_points = 0;
    int max_points = 10;
    *point_list = malloc(max_points*sizeof(struct Point));
    *num_edges = 0;
    int max_edges = 10;
    *edge_list = malloc(max_edges*sizeof(struct Edge));
    *num_cells = 0;
    int max_cells = 10;
    *cell_list = malloc(max_cells*sizeof(struct Cell));
    char *type = calloc(3, sizeof(char));
    double x, y;
    int index1, index2, index3, index4;
    int id;
    type[0] = 60;
    int itter = 0;


    if (!infile) {
        fprintf(stderr, "Something went wrong when reading %s", filename);
        exit(2);
    }
    while (fscanf(infile, "%c %d", type, &id) != EOF)
    {
        if(type[0] == 'p'){
            fscanf(infile, "%lg %lg", &x, &y);
            add_point(id, x, y, point_list, num_points, &max_points);
            fscanf(infile, "%*c");
        }else if(type[0] == 'e'){
            fscanf(infile, "%s %d %d", type, &index1, &index2);
            add_edge(id, type[0], index1, index2, edge_list, num_edges, &max_edges);
            fscanf(infile, "%*c");
        }else if(type[0] == 'c'){
            fscanf(infile, "%d %d %d %d", &index1, &index2, &index3, &index4);
            add_cell(id, index1, index2, index3, index4, cell_list, num_cells, &max_cells);
            fscanf(infile, "%*c");
        }else{
            //error
            fprintf(stderr, "num points: %d num edges: %d num cells: %d\n", *num_points, *num_edges, *num_cells);
            fprintf(stderr, "itter: %d type: %d id: %d\n", itter, (int)(type[0]), id);
            fprintf(stderr, "Error reading input. Exiting.\n");
            exit(2);
        }
        itter++;
        //add2grid(type1[0], type2[0], type3[0], x1, y1, x2, y2, x3, y3, edge_list, num_edges, &max_edges, map, num_cells, &max_cells);
    }

    fclose(infile);
	free(type);
}

int hasSource(struct Cell aCell, struct Edge *edge_list){
    if(edge_list[aCell.edge1].type == 's'){
        return aCell.edge1;
    }
    if(edge_list[aCell.edge2].type == 's'){
        return aCell.edge2;
    }
    if(edge_list[aCell.edge3].type == 's'){
        return aCell.edge3;
    }
    if(edge_list[aCell.edge4].type == 's'){
        return aCell.edge4;
    }
    return -1;
}

int hasWall(struct Cell aCell, struct Edge *edge_list){
    if(edge_list[aCell.edge1].type == 'w'){
        return aCell.edge1;
    }
    if(edge_list[aCell.edge2].type == 'w'){
        return aCell.edge2;
    }
    if(edge_list[aCell.edge3].type == 'w'){
        return aCell.edge3;
    }
    if(edge_list[aCell.edge4].type == 'w'){
        return aCell.edge4;
    }
    return -1;
}

double findLength(struct Point *point_list, int num_points, double y_pt){
    double min_x = DBL_MAX;
    double max_x = DBL_MIN;
    for(int i = 0; i < num_points; i++){
        if(point_list[i].y == y_pt){
            if(point_list[i].x < min_x){
                min_x = point_list[i].x;
            }
            if(point_list[i].x > max_x){
                max_x = point_list[i].x;
            }
        }
    }
    return max_x - min_x;
}

double findHeight(struct Point *point_list, int num_points, double x_pt){
    double min_y = DBL_MAX;
    double max_y = DBL_MIN;
    for(int i = 0; i < num_points; i++){
        if(point_list[i].x == x_pt){
            if(point_list[i].y < min_y){
                min_y = point_list[i].y;
            }
            if(point_list[i].y > max_y){
                max_y = point_list[i].y;
            }
        }
    }
    return max_y - min_y;
}

int main ( int argc, char **argv ){
    validate_inputs(argc, argv);
    clock_t start, end;
    double cpu_time_used;

    start = clock();

    int num_cells, num_edges, num_points;
    struct Cell *cell_list;
    struct Edge *edge_list;
    struct Point *point_list;
    double v_in, p_in, p_out, p_grad, mu;

    num_cells = 0;
    num_edges = 0;
    num_points = 0;
    cell_list = 0;
    edge_list = 0;
    point_list = 0;

    read_inputs(argv[1], &cell_list, &num_cells, &edge_list, &num_edges, &point_list, &num_points);

    v_in = atoi(argv[2]);
    p_in = atoi(argv[3]);
    p_out = atoi(argv[4]);
    p_grad = p_in - p_out;
    mu = 5;

    double *p_field_mag = malloc(num_points*sizeof(double));
    double *v_field_mag = malloc(num_points*sizeof(double));
    double *p_field_dir = malloc(num_points*sizeof(double));
    double *v_field_dir = malloc(num_points*sizeof(double));
    char *p_field_valid = calloc(num_points, sizeof(char));
    char *v_field_valid = calloc(num_points, sizeof(char));

    /*
    //inital conditions
    for(int i = 0; i < num_edges; i++){
        if(edge_list[i].type == 's'){
            int point1 = edge_list[i].point1;
            int point2 = edge_list[i].point2;
            p_field_mag[point1] = p_in;
            p_field_dir[point1] = 0.0;
            v_field_mag[point1] = v_in;
            v_field_dir[point1] = 0.0;
            p_field_valid[point1] = 1;
            v_field_valid[point1] = 1;
            p_field_mag[point2] = p_in;
            p_field_dir[point2] = 0.0;
            v_field_mag[point2] = v_in;
            v_field_dir[point2] = 0.0;
            p_field_valid[point2] = 1;
            v_field_valid[point2] = 1;
        }
        //no slip condition
        if(edge_list[i].type == 'w'){
            int point1 = edge_list[i].point1;
            int point2 = edge_list[i].point2;
            v_field_mag[point1] = 0.0;
            v_field_dir[point1] = 0.0;
            v_field_valid[point1] = 1;
            v_field_mag[point2] = 0.0;
            v_field_dir[point2] = 0.0;
            v_field_valid[point2] = 1;
        }
    }
    */

    //because this is steady state parallel flow b/t plates, these are constant
    double height = findHeight(point_list, num_points, 0.0);
    double length = findLength(point_list, num_points, 0.0);
    for(int i = 0; i < num_points; i++){
        p_field_mag[i] = p_in - p_grad*((point_list[i].x)/length);
        p_field_dir[i] = 0.0;
        v_field_mag[i] = p_grad/(2*mu)*point_list[i].y*(height-point_list[i].y);
        v_field_dir[i] = 0.0;
        p_field_valid[i] = 1;
        v_field_valid[i] = 1;
    }

    FILE *outfile;
    outfile = fopen("./results.txt", "w");

    if (!outfile) {
        fprintf(stderr, "Something went wrong when opening outfile\n");
        exit(3);
    }

    fprintf(outfile, "%s\t%s\t%s\t%s\n", "x", "y", "pressure", "velocity");
    for(int i = 0; i < num_points; i++){
        fprintf(outfile, "%g\t%g\t%g\t\t%g\n", point_list[i].x, point_list[i].y, p_field_mag[i], v_field_mag[i]);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    fclose(outfile);
    fprintf(stdout, "Took %g seconds.\n",cpu_time_used);
}
