/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <iostream>
#include <math.h>
#include "mex.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iterator>
#include <random>
#include <algorithm>
#include <queue>
#include <stack>
#include <limits>

using namespace std;

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTCONNECT  1
#define RRTSTAR     2
#define PRM         3

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10

typedef struct {
  int X1, Y1;
  int X2, Y2;
  int Increment;
  int UsingYIndex;
  int DeltaX, DeltaY;
  int DTerm;
  int IncrE, IncrNE;
  int XIndex, YIndex;
  int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size)
{
    double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize));
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params)
{
  params->UsingYIndex = 0;

  if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
    (params->UsingYIndex)++;

  if (params->UsingYIndex)
    {
      params->Y1=p1x;
      params->X1=p1y;
      params->Y2=p2x;
      params->X2=p2y;
    }
  else
    {
      params->X1=p1x;
      params->Y1=p1y;
      params->X2=p2x;
      params->Y2=p2y;
    }

   if ((p2x - p1x) * (p2y - p1y) < 0)
    {
      params->Flipped = 1;
      params->Y1 = -params->Y1;
      params->Y2 = -params->Y2;
    }
  else
    params->Flipped = 0;

  if (params->X2 > params->X1)
    params->Increment = 1;
  else
    params->Increment = -1;

  params->DeltaX=params->X2-params->X1;
  params->DeltaY=params->Y2-params->Y1;

  params->IncrE=2*params->DeltaY*params->Increment;
  params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
  params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

  params->XIndex = params->X1;
  params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y)
{
  if (params->UsingYIndex)
    {
      *y = params->XIndex;
      *x = params->YIndex;
      if (params->Flipped)
        *x = -*x;
    }
  else
    {
      *x = params->XIndex;
      *y = params->YIndex;
      if (params->Flipped)
        *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params)
{
  if (params->XIndex == params->X2)
    {
      return 0;
    }
  params->XIndex += params->Increment;
  if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
    params->DTerm += params->IncrE;
  else
    {
      params->DTerm += params->IncrNE;
      params->YIndex += params->Increment;
    }
  return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
		   int x_size,
 		   int y_size)

{
	bresenham_param_t params;
	int nX, nY; 
    short unsigned int nX0, nY0, nX1, nY1;

    //printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
    
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

    //printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
            return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
		   int x_size, int y_size)
{
    double x0,y0,x1,y1;
    int i;
    
 	//iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
    y1 = 0;
	for(i = 0; i < numofDOFs; i++)
	{
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
				return 0;
	}    
    return 1;
}

//======================================================================================================================

struct configuration{
//    int n_dof;
    vector<double> angles;

    friend bool operator== (const configuration &c1, const configuration &c2);
    friend bool operator!= (const configuration &c1, const configuration &c2);

    configuration(){}

    configuration(vector<double> config_angles):
    angles(config_angles){
    }

    configuration(int num_dof,double* config_angles)
    {
        angles.reserve(num_dof);
        for (double *it = config_angles; it != config_angles + num_dof; ++it) {
            angles.push_back(*it);
        }
    }

    configuration(const configuration &c2) {angles=c2.angles;}

    void print_config() const
    {
        for(size_t i=0;i<angles.size();i++)
            cout<<angles[i]<<"\t";
        cout<<endl;
    }

    double get_distance(const configuration &c2) const
    {
        assert(angles.size()==c2.angles.size());
        double dist = 0;
        for(size_t i=0;i<angles.size();i++)
        {
            dist+= std::min(pow((angles[i]-c2.angles[i]),2),pow((angles[i]-(c2.angles[i]-2*PI)),2));
        }
        //Square root can be put, but not needed as we are just comparing
        return dist;
    }

};

bool operator== (const configuration &c1, const configuration &c2)
{
    return (c1.angles == c2.angles);
}

bool operator!= (const configuration &c1, const configuration &c2)
{
    return !(c1== c2);
}

//======================================================================================================================

struct Node{
    configuration c;
    vector<int> neighbors;

    Node(){}

    Node (configuration state_config):
    c(state_config)
    {}

    void print_Node() const {
        c.print_config();
        cout<<"Neighbors are: "<<endl;
        for(size_t i=0;i<neighbors.size();i++)
        {
            cout<<neighbors[i]<<"\t";
        }
        cout<<endl;
    }
};

//======================================================================================================================

struct RRT_Node{
    configuration c;
    int parent;

    RRT_Node(){}

    RRT_Node (configuration state_config):
            c(state_config),
            parent(-1)
    {}

    RRT_Node (configuration state_config, int node_parent):
            c(state_config),
            parent(node_parent)
    {}

    void print_RRT_Node() const {
        c.print_config();
        cout<<"Parent is: "<<parent<<endl;
    }
};

//======================================================================================================================

struct RRT_Star_Node{
    configuration c;
    int parent;
    double gcost;

    RRT_Star_Node(){}

    RRT_Star_Node (configuration state_config):
            c(state_config),
            parent(-1),
            gcost(INT_MAX)
    {}

    RRT_Star_Node (configuration state_config, int node_parent,double new_gcost):
            c(state_config),
            parent(node_parent),
            gcost(new_gcost)
    {}

    void print_RRT_Star_Node() const {
        c.print_config();
        cout<<"Parent is: "<<parent<<endl;
        cout<<"g_cost is: "<<gcost<<endl;
    }
};

//======================================================================================================================

int interpolation_based_plan(   double*	map,
                                const int &x_size,
                                const int &y_size,
                                const double* armstart_anglesV_rad,
                                const double* armgoal_anglesV_rad,
                                const int &numofDOFs,
                                double*** plan)
{
    //for now just do straight interpolation between start and goal checking for the validity of samples

    double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
            distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));
    if(numofsamples < 2){
        printf("the arm is already at the goal\n");
        return 0;
    }
    *plan = (double**) malloc(numofsamples*sizeof(double*));
    int firstinvalidconf = 1;
    for (i = 0; i < numofsamples; i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double));
        for(j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i)/(numofsamples-1))*(armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
        }
        if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size) && firstinvalidconf)
        {
            firstinvalidconf = 1;
            printf("ERROR: Invalid arm configuration!!!\n");
        }
    }
    return numofsamples;
}

//======================================================================================================================

double * generate_random_config(const int &numofDOFs)
{
    const double MIN_ANGLE = 0;
    const double MAX_ANGLE = 2*PI;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::uniform_real_distribution<double> distribution(MIN_ANGLE,MAX_ANGLE);
    double *random_config = new double[numofDOFs];
    for(size_t i=0;i<numofDOFs;i++)
    {
        random_config[i] = distribution(generator);
    }
    return std::move(random_config);
}

//======================================================================================================================

struct less1{
    bool operator()(const pair<double,int> &a,const pair<double,int> &b) const{
        return a.first>b.first;
    }
};

//======================================================================================================================

vector<int> get_neighbors(const configuration &c,
                           unordered_map<int,Node> road_map,
                           const int &nearest_neighbors_to_consider)
{
    if(road_map.size()==1)
        return vector<int> {};

    vector<int> neighbors;
    if(road_map.size()<=nearest_neighbors_to_consider+1)
    {
        for(int i=0;i<road_map.size()-1;i++)
        {
            neighbors.push_back(i);
        }
        return std::move(neighbors);
    }

    vector<pair<double,int>> dist_index;
    for(int i=0;i<road_map.size()-1;i++)
    {
//        auto test_config = road_map[i].c;
        dist_index.emplace_back(make_pair(c.get_distance(road_map[i].c),i));
    }

    std::make_heap(dist_index.begin(), dist_index.end(), less1());
    int count = 0;
    while(count<nearest_neighbors_to_consider)
    {
        neighbors.push_back(dist_index.front().second);
        std::pop_heap(dist_index.begin(),dist_index.end(),less1());
        dist_index.pop_back();
        count++;
    }
    return std::move(neighbors);
}

//======================================================================================================================

bool is_connection_possible(  const configuration &start_config,
                              const configuration &end_config,
                              const int &numofDOFs,
                              double*	map,
                              const int &x_size,
                              const int &y_size )
{
    double distance = 0;

    for (int j = 0; j < numofDOFs; j++){
        if(distance < fabs(start_config.angles[j] - end_config.angles[j]))
            distance = fabs(start_config.angles[j] - end_config.angles[j]);
    }

    int numofsamples = (int)(distance/(PI/20));
    for (int i = 0; i < numofsamples; i++){
        double *intermediate_config = new double[numofDOFs];
        for(int j = 0; j < numofDOFs; j++){
            intermediate_config[j] = start_config.angles[j] + ((double)(i)/(numofsamples-1))*(end_config.angles[j] - start_config.angles[j]);
        }
        if(!IsValidArmConfiguration(intermediate_config, numofDOFs, map, x_size, y_size))
        {
            return false;
        }
    }
    return true;
}

//======================================================================================================================

unordered_map<int,Node> add_edges(unordered_map<int,Node> road_map,
                                  int sample_count,
                                  vector<int> neighbors,
                                  const int &numofDOFs,
                                  double*	map,
                                  const int &x_size,
                                  const int &y_size)
{
    for(auto neighbor:neighbors)
    {
        if(is_connection_possible(road_map[neighbor].c,road_map[sample_count].c,numofDOFs,map,x_size,y_size))
        {
            road_map[sample_count].neighbors.push_back(neighbor);
            road_map[neighbor].neighbors.push_back(sample_count);
        }
    }
    return std::move(road_map);
}
//======================================================================================================================

unordered_map<int,Node> build_road_map(double*	map,
                                        const int &x_size,
                                        const int &y_size,
                                        const int &numofDOFs,
                                        const int &num_samples_PRM)
{
    unordered_map<int,Node> road_map;
    int sample_count=0;
    int firstinvalidconf = 1;
    int NEAREST_NEIGHBORS_TO_CONSIDER = 3;
    while(sample_count<num_samples_PRM)
    {
        auto random_config = generate_random_config(numofDOFs);
        if(IsValidArmConfiguration(random_config, numofDOFs, map, x_size, y_size))
        {   configuration c(numofDOFs,random_config);
            road_map[sample_count] = Node(c);
            const auto neighbors = get_neighbors(c,road_map,NEAREST_NEIGHBORS_TO_CONSIDER);
//            for(size_t i=0;i<neighbors.size();i++)
//            {
//                cout<<"Neighbors_indices: "<<neighbors[i]<<"\t";
//            }
//            cout<<endl;
            road_map = add_edges(std::move(road_map),sample_count,neighbors,numofDOFs,map,x_size,y_size);
            sample_count++;
        }
    }

    /// Visualize Road_Map
//    for(const auto &elt:road_map)
//    {
//        cout<<"Node index: "<<elt.first<<endl;
//        elt.second.print_Node();
//    }

    return std::move(road_map);
}

//======================================================================================================================

bool add_start_and_goal_to_road_map(double*	map,
                                                       const int &x_size,
                                                       const int &y_size,
                                                       const int &numofDOFs,
                                                       unordered_map<int,Node> &road_map,
                                                       const configuration &start_config,
                                                       const configuration &goal_config)
{
    int road_map_length = road_map.size();
    int flag = 0;
    for(const auto &elt:road_map)
        {
            if(is_connection_possible(start_config,elt.second.c,numofDOFs,map,x_size,y_size))
            {
                road_map[road_map_length]  = Node {start_config};
                road_map[road_map_length].neighbors.push_back(elt.first);
                road_map[elt.first].neighbors.push_back(road_map_length);
                cout<<"Start added to Graph"<<endl;
                road_map_length++;
                flag=1;
                break;
            }
        }

    if(!flag)
        return false;

    flag = 0;
    for(const auto &elt:road_map)
    {
        if(is_connection_possible(goal_config,elt.second.c,numofDOFs,map,x_size,y_size))
        {
            road_map[road_map_length]  = Node {goal_config};
            road_map[road_map_length].neighbors.push_back(elt.first);
            road_map[elt.first].neighbors.push_back(road_map_length);
            cout<<"Goal added to Graph"<<endl;
            flag=1;
            break;
        }
    }

    if(!flag)
        return false;

    return true;
}

//======================================================================================================================

void dfs_util(unordered_map<int,Node> road_map,
              unordered_set<int> &visited,
              stack<int> &s,
              const int &start_index,
              const int &goal_index,
              int &flag)
{
    visited.insert(start_index);
    s.push(start_index);

    if(start_index==goal_index)
    {
        cout<<"Goal found"<<endl;
        flag=1;
        return;
    }

    for(size_t i=0;i<road_map[start_index].neighbors.size();i++)
    {
        if(!visited.count(road_map[start_index].neighbors[i]))
        {
            dfs_util(road_map,visited,s,road_map[start_index].neighbors[i],goal_index,flag);
            if(flag)
                return;
        }
    }
    s.pop();
}

//======================================================================================================================

void get_plan_from_stack(unordered_map<int,Node> road_map,
              double*** plan,
              stack<int> s,
              const int &numofDOFs)
{
    stack<int> reversed_stack;
    while(!s.empty())
    {
        reversed_stack.push(s.top());
        s.pop();
    }

    int i=0;
    *plan = (double**) malloc(reversed_stack.size()*sizeof(double*));
    while(!reversed_stack.empty())
    {
//        cout<<"Vertex is: "<<reversed_stack.top()<<endl;
        (*plan)[i] = (double *) malloc(numofDOFs * sizeof(double));
        for (int j = 0; j < numofDOFs; j++) {
            (*plan)[i][j] = road_map[reversed_stack.top()].c.angles[j];
        }
        reversed_stack.pop();
        i++;
    }
}

//======================================================================================================================

int search_road_map(unordered_map<int,Node> road_map,
                    double*** plan,
                    const int &numofDOFs)
{
    int start_index = road_map.size()-2;
    int goal_index = road_map.size()-1;
    unordered_set<int> visited;
    stack<int> s;
    int flag = 0;
    dfs_util(road_map,visited,s,start_index,goal_index,flag);
    get_plan_from_stack(road_map,plan,s,numofDOFs);
    return s.size();
}


//======================================================================================================================

int PRM_planner(double*	map,
        const int &x_size,
        const int &y_size,
        double* armstart_anglesV_rad,
        double* armgoal_anglesV_rad,
        const int &numofDOFs,
        double*** plan,
        int num_samples_PRM)
{
  /*
   Use an unordered_map<configuration,pair<vector<configuration>,int>>. The int here would sort of represent the connected component
   that configuration is a part of. The vector<configuration> is are the configuration which are it's edges. Graphs should be
   represented as undirected. Intialise the neighbpurhood as 1 when making the graph initially, and then assign the same if a
   neighbor is found else allot it a present_count+1 number of component.
   Now this graph can be searched by BFS,DFS, Dijkstra
   */

  const configuration start_config(numofDOFs,armstart_anglesV_rad);
  const configuration goal_config(numofDOFs,armgoal_anglesV_rad);
  int iteration_number = 1;
  while(true)
          {
              cout<<"Iteration Number: "<<iteration_number<<endl;
              auto road_map = build_road_map(map,x_size,y_size,numofDOFs,num_samples_PRM);
//              cout<<"Road_Map_length = "<<road_map.size()<<endl;
              auto got_connected = add_start_and_goal_to_road_map(map,x_size,y_size,numofDOFs,road_map,start_config,goal_config);
//              cout<<"Road_Map_length = "<<road_map.size()<<endl;

              auto plan_length = search_road_map(road_map,plan,numofDOFs);
              if(plan_length && got_connected)
              {
//                  cout<<"Plan_Length: "<<plan_length<<endl;
                  return plan_length;
              }

              num_samples_PRM *=2;
              cout<<"New number of samples are: "<<num_samples_PRM<<endl;
              iteration_number++;
          }
}

//======================================================================================================================

template <typename T>
vector<int> get_k_nearest_neighbors(const configuration &c,
                                    unordered_map<int,T> Tree,
                                    int nearest_neighbors_to_consider)
{
    vector<pair<double,int>> dist_index;
//    for(int i=0;i<Tree.size();i++)
//    {
//        dist_index.emplace_back(make_pair(c.get_distance(Tree[i].c),i));
//    }
    for(const auto &elt:Tree)
    {
        dist_index.emplace_back(make_pair(c.get_distance(elt.second.c),elt.first));
    }
    vector<int> neighbors;
    std::make_heap(dist_index.begin(), dist_index.end(), less1());
    int count = 0;
    while(count<nearest_neighbors_to_consider && !dist_index.empty())
    {
        neighbors.push_back(dist_index.front().second);
        std::pop_heap(dist_index.begin(),dist_index.end(),less1());
        dist_index.pop_back();
        count++;
    }
    return std::move(neighbors);
}

//======================================================================================================================

template <typename T>
configuration get_epsilon_connect_config(unordered_map<int,T> Tree,
                                        const configuration &start_config,
                                        configuration end_config,
                                        const int &numofDOFs,
                                        double*	map,
                                        const int &x_size,
                                        const int &y_size,
                                        const double &epsilon,
                                        const int &NUM_OF_SAMPLES)
{
    double distance = 0;

    for (int j = 0; j < numofDOFs; j++){
        if(distance < fabs(start_config.angles[j] - end_config.angles[j]))
            distance = fabs(start_config.angles[j] - end_config.angles[j]);
    }

    if(fabs(distance)>epsilon)
    {
        configuration new_end_configuration = end_config;
        for(int j = 0; j < numofDOFs; j++)
            new_end_configuration.angles[j] = start_config.angles[j] + epsilon*(1/distance)*(end_config.angles[j]-start_config.angles[j]);
        end_config = new_end_configuration;
    }

    auto previous_valid_configuration = start_config;
    for(int i=1;i<=NUM_OF_SAMPLES;i++)
    {
        double *intermediate_config = new double[numofDOFs];
        for(int j = 0; j < numofDOFs; j++){
            intermediate_config[j] = start_config.angles[j] + ((double)(i)/(NUM_OF_SAMPLES))*(end_config.angles[j] - start_config.angles[j]);
        }
        if(!IsValidArmConfiguration(intermediate_config, numofDOFs, map, x_size, y_size))
        {
            break;
        }
       previous_valid_configuration = configuration{numofDOFs,intermediate_config};
    }

    return previous_valid_configuration;

}

//======================================================================================================================

template <typename T>
bool epsilon_connect(const int &nearest_neighbor_index,
                     const configuration &random_config,
                     unordered_map<int,T> &Tree,
                     const double &epsilon,
                     const int &numofDOFs,
                     double*	map,
                     const int &x_size,
                     const int &y_size,
                     const int &sample_count,
                     const configuration &goal_config,
                     const int &discretization_factor=10)
{
    const auto q_new = get_epsilon_connect_config(Tree,Tree[nearest_neighbor_index].c,random_config,numofDOFs,map,x_size,y_size,epsilon,discretization_factor);
    if(q_new!=Tree[nearest_neighbor_index].c)
    {
        Tree[sample_count] = RRT_Node(q_new,nearest_neighbor_index);
//        if(random_config==goal_config)
//            cout<<"Epsilon connection made"<<endl;
        return true;
    }
//    if(random_config==goal_config)
//        cout<<"Failed to make epsilon connection"<<endl;
    return false;

}

//======================================================================================================================

template <typename T>
bool epsilon_connect_RRTstar(const int &nearest_neighbor_index,
                             const configuration &random_config,
                             unordered_map<int,T> &Tree,
                             const double &epsilon,
                             const int &numofDOFs,
                             double*	map,
                             const int &x_size,
                             const int &y_size,
                             const int &sample_count,
                             const configuration &goal_config,
                             const int &discretization_factor=10)
{
    /// This and the function above are entirely common. The only difference is the node type. Make common function.
    const auto q_new = get_epsilon_connect_config(Tree,Tree[nearest_neighbor_index].c,random_config,numofDOFs,map,x_size,y_size,epsilon,discretization_factor);
    if(q_new!=Tree[nearest_neighbor_index].c)
    {
        auto distance_bw_configs = q_new.get_distance(Tree[nearest_neighbor_index].c);
        Tree[sample_count] = RRT_Star_Node(q_new,nearest_neighbor_index,Tree[nearest_neighbor_index].gcost + std::move(distance_bw_configs));
        return true;
    }
    return false;

}

//======================================================================================================================

int change_path_vector_to_path(double***plan,
                          const vector<configuration> &path_vector,
                          const int &numofDOFs)
{
    cout<<"Vector size is" <<path_vector.size()<<endl;
    *plan = (double**) malloc(path_vector.size()*sizeof(double*));
    for (int i = 0; i < path_vector.size(); i++) {
        (*plan)[i] = (double *) malloc(numofDOFs * sizeof(double));
        for (int j = 0; j < numofDOFs; j++) {
            (*plan)[i][j] = path_vector[i].angles[j];
        }
    }
    return path_vector.size();
}
//======================================================================================================================

template <typename T>
vector<configuration> get_path_vector(
               unordered_map<int,T> Tree,
               const int &sample_count,
               const int &numofDOFs,
               const configuration &start_config)
{
    auto present_node = Tree[sample_count];
    vector<configuration> vec;
    while(present_node.c!=start_config)
    {
        vec.emplace_back(present_node.c);
        present_node = Tree[present_node.parent];
    }
    vec.emplace_back(present_node.c);
    cout<<"Path vector formed"<<endl;
    return std::move(vec);
}

//======================================================================================================================

int base_cases_RRT_family(double*	map,
                          const int &x_size,
                          const int &y_size,
                          double* armstart_anglesV_rad,
                          double* armgoal_anglesV_rad,
                          const int &numofDOFs,
                          double*** plan,
                          const configuration &start_config,
                          const configuration &goal_config)

{
    if(!IsValidArmConfiguration(armgoal_anglesV_rad, numofDOFs, map, x_size, y_size))
    {
        cout<<"Goal Position is invalid "<<endl;
        return 0;
    }

    if(!IsValidArmConfiguration(armstart_anglesV_rad, numofDOFs, map, x_size, y_size))
    {
        cout<<"Start Position is invalid "<<endl;
        return 0;
    }

    if(start_config==goal_config)
    {
        *plan = (double**) malloc(2*sizeof(double*));
        (*plan)[0] = (double *) malloc(numofDOFs * sizeof(double));
        (*plan)[0] = armstart_anglesV_rad;
        (*plan)[1] = (double *) malloc(numofDOFs * sizeof(double));
        (*plan)[1] = armgoal_anglesV_rad;
        return 2;
    }

    return -1;
}

//======================================================================================================================

bool is_new_sample_in_goal_region(const configuration &new_tree_config,
                                  const configuration &goal_config,
                                  const double &goal_region_threshold)
{
    assert(new_tree_config.angles.size()==goal_config.angles.size());
    for(size_t i=0;i<goal_config.angles.size();i++)
    {
        if(goal_config.angles[i] - goal_region_threshold > new_tree_config.angles[i] || goal_config.angles[i] + goal_region_threshold < new_tree_config.angles[i])
            return false;
    }
    return true;
}

//======================================================================================================================

template <typename T>
void rewire_RRT_star(const int &child_index,
                     const int &potential_parent_index,
                     unordered_map<int,T> &Tree,
                     double*	map,
                     const int &x_size,
                     const int &y_size,
                     const int &numofDOFs)
{
    const auto potential_parent_config = Tree[potential_parent_index].c;
    const auto child_config = Tree[child_index].c;
    if(is_connection_possible(potential_parent_config,child_config,numofDOFs,map,x_size,y_size))
    {
        const auto transition_cost = potential_parent_config.get_distance(child_config);
        if(Tree[potential_parent_index].gcost + transition_cost < Tree[child_index].gcost)
            {
                Tree[child_index].gcost = Tree[potential_parent_index].gcost + transition_cost;
                Tree[child_index].parent = potential_parent_index;
            }
    }
}

//======================================================================================================================

int RRT_planner(double*	map,
                const int &x_size,
                const int &y_size,
                double* armstart_anglesV_rad,
                double* armgoal_anglesV_rad,
                const int &numofDOFs,
                double*** plan,
                int num_samples_RRT)
{
    const double EPSILON = PI/25;
    const int NEAREST_NEIGHBORS_TO_CONSIDER = 1;
    double GOAL_BIAS = 0.01;
    const int &DISCRETIZATION_FACTOR = 10;  //Interpolation for collision checker
    const double GOAL_REGION_THRESHOLD = PI/36;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1

    const configuration start_config(numofDOFs,armstart_anglesV_rad);
    const configuration goal_config(numofDOFs,armgoal_anglesV_rad);

    auto base_case_result = base_cases_RRT_family(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan,start_config,goal_config);
    if(base_case_result!=-1)
        return base_case_result;

    unordered_map<int,RRT_Node> Tree;
    Tree[0] = RRT_Node(start_config);
    int sample_count = 1;
    bool is_goal_reached = false;
    while(sample_count<num_samples_RRT)
    {
        if(sample_count%1000==0)
            cout<<"Sample count is: "<<sample_count<<endl;

        if(sample_count>10000)
            GOAL_BIAS = 0.05;

        configuration random_config;
        if(dis(gen)<GOAL_BIAS)
        {
            random_config = goal_config;
        }
        else
        {
            auto random_state = generate_random_config(numofDOFs);
            random_config = configuration{numofDOFs,random_state};
        }
        const auto k_nearest_neighbors = get_k_nearest_neighbors(random_config,Tree,NEAREST_NEIGHBORS_TO_CONSIDER);
        const auto nearest_neighbor_index = k_nearest_neighbors[0];
        auto was_epsilon_connection_made = epsilon_connect(nearest_neighbor_index,random_config,Tree,EPSILON,numofDOFs,map,x_size,y_size,sample_count,goal_config,DISCRETIZATION_FACTOR);
        if(was_epsilon_connection_made)
        {
            if(is_new_sample_in_goal_region(Tree[sample_count].c,goal_config,GOAL_REGION_THRESHOLD))
            {
                cout<<"Connection to goal region made"<<endl;
                Tree[sample_count+1] = RRT_Node(goal_config,sample_count);  //Adding goal configuration to my tree
                is_goal_reached = true;
                sample_count++;     //Adding this for goal. So that backtrack starts at goal.
                break;
            }
            sample_count++;
        }
    }

//    /// Visualize Tree
//    for(const auto &elt:Tree)
//    {
//        cout<<"Node index: "<<elt.first<<endl;
//        elt.second.print_RRT_Node();
//    }

    if(Tree[sample_count].c == goal_config)
        cout<<"Goal has been correctly reached"<<endl;

    if(is_goal_reached)
    {
        auto path_vector = get_path_vector(Tree,sample_count,numofDOFs,start_config);
        reverse(path_vector.begin(), path_vector.end());
        return change_path_vector_to_path(plan,path_vector,numofDOFs);
    }

    cout<<"Could not reach goal. Sample more points"<<endl;
    return -1;
}

//======================================================================================================================

int RRT_connect(double*	map,
                const int &x_size,
                const int &y_size,
                double* armstart_anglesV_rad,
                double* armgoal_anglesV_rad,
                const int &numofDOFs,
                double*** plan,
                int num_samples_RRT_connect)
{
    const double EPSILON = PI/25;
    const double CONNECT_EPSILON = INT_MAX;
    const int NEAREST_NEIGHBORS_TO_CONSIDER = 1;
    const int &DISCRETIZATION_FACTOR = 20;  //Interpolation for collision checker

    const configuration start_config(numofDOFs,armstart_anglesV_rad);
    const configuration goal_config(numofDOFs,armgoal_anglesV_rad);

    auto base_case_result = base_cases_RRT_family(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan,start_config,goal_config);
    if(base_case_result!=-1)
        return base_case_result;

    unordered_map<int,RRT_Node> forward_Tree;
    unordered_map<int,RRT_Node> backward_Tree;
    forward_Tree[0] = RRT_Node(start_config);
    backward_Tree[0] = RRT_Node(goal_config);
    int sample_count = 1;
    bool is_goal_reached = false;
    while(sample_count<num_samples_RRT_connect)
    {
        configuration random_config;
        auto random_state = generate_random_config(numofDOFs);
        random_config = configuration{numofDOFs,random_state};
        if(sample_count%2!=0)   //Expand forward_tree
        {
            const auto k_nearest_neighbors = get_k_nearest_neighbors(random_config,forward_Tree,NEAREST_NEIGHBORS_TO_CONSIDER);
            const auto nearest_neighbor_index = k_nearest_neighbors[0];
            auto was_epsilon_connection_made = epsilon_connect(nearest_neighbor_index,random_config,forward_Tree,EPSILON,numofDOFs,map,x_size,y_size,sample_count,goal_config,DISCRETIZATION_FACTOR);
            if(was_epsilon_connection_made)
            {
                const auto k_nearest_neighbors = get_k_nearest_neighbors(forward_Tree[sample_count].c,backward_Tree,NEAREST_NEIGHBORS_TO_CONSIDER);
                const auto nearest_neighbor_index = k_nearest_neighbors[0];
                auto did_both_trees_connect = epsilon_connect(nearest_neighbor_index,forward_Tree[sample_count].c,backward_Tree,CONNECT_EPSILON,numofDOFs,map,x_size,y_size,sample_count,goal_config,DISCRETIZATION_FACTOR);
                if(did_both_trees_connect)
                {
                    cout<<"Both trees got connected A"<<endl;
                    is_goal_reached = true;
                    break;
                }
            }
        }

        else
        {
            const auto k_nearest_neighbors = get_k_nearest_neighbors(random_config,backward_Tree,NEAREST_NEIGHBORS_TO_CONSIDER);
            const auto nearest_neighbor_index = k_nearest_neighbors[0];
            auto was_epsilon_connection_made = epsilon_connect(nearest_neighbor_index,random_config,backward_Tree,EPSILON,numofDOFs,map,x_size,y_size,sample_count,goal_config,DISCRETIZATION_FACTOR);
            if(was_epsilon_connection_made)
            {
                const auto k_nearest_neighbors = get_k_nearest_neighbors(backward_Tree[sample_count].c,forward_Tree,NEAREST_NEIGHBORS_TO_CONSIDER);
                const auto nearest_neighbor_index = k_nearest_neighbors[0];
                auto did_both_trees_connect = epsilon_connect(nearest_neighbor_index,backward_Tree[sample_count].c,forward_Tree,CONNECT_EPSILON,numofDOFs,map,x_size,y_size,sample_count,goal_config,DISCRETIZATION_FACTOR);
                if(did_both_trees_connect)
                {
                    cout<<"Both trees got connected B"<<endl;
                    is_goal_reached = true;
                    break;
                }
            }
        }
        sample_count++;
    }
    if(is_goal_reached)
    {
        cout<<"Goal has been reached"<<endl;
        auto path_vector_A= get_path_vector(forward_Tree,sample_count,numofDOFs,start_config);
        std::reverse(path_vector_A.begin(),path_vector_A.end());
        path_vector_A.pop_back();
        auto path_vector_B= get_path_vector(backward_Tree,sample_count,numofDOFs,goal_config);
        for(auto elt:path_vector_B)
        {
            path_vector_A.push_back(elt);
        }
        return change_path_vector_to_path(plan,path_vector_A,numofDOFs);
    }

    cout<<"Could not reach goal. Sample more points"<<endl;
    return -1;
}

//======================================================================================================================

int RRT_star(double*	map,
                const int &x_size,
                const int &y_size,
                double* armstart_anglesV_rad,
                double* armgoal_anglesV_rad,
                const int &numofDOFs,
                double*** plan,
                int num_samples_RRT_star)
 {
    const double EPSILON = PI / 25;
    const int NEAREST_NEIGHBORS_TO_CONSIDER = 1;
    const int NEAREST_NEIGHBORS_FOR_REWIRING = 3;
    double GOAL_BIAS = 0.01;
    const int &DISCRETIZATION_FACTOR = 10;  //Interpolation for collision checker
    const double GOAL_REGION_THRESHOLD = PI/36;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1

    const configuration start_config(numofDOFs, armstart_anglesV_rad);
    const configuration goal_config(numofDOFs, armgoal_anglesV_rad);

    auto base_case_result = base_cases_RRT_family(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad,
                                                  numofDOFs, plan, start_config, goal_config);
    if (base_case_result != -1)
        return base_case_result;

    unordered_map<int, RRT_Star_Node> Tree;
    Tree[0] = RRT_Star_Node(start_config);
    Tree[0].gcost = 0;
    int sample_count = 1;
    bool is_goal_reached = false;
    while (sample_count < num_samples_RRT_star) {
        if (sample_count % 1000 == 0)
            cout << "Sample count is: " << sample_count << endl;

        if (sample_count > 10000)
            GOAL_BIAS = 0.05;

        configuration random_config;
        if (dis(gen) < GOAL_BIAS) {
            random_config = goal_config;

        } else {
            auto random_state = generate_random_config(numofDOFs);
            random_config = configuration{numofDOFs, random_state};
        }

        const auto k_nearest_neighbors = get_k_nearest_neighbors(random_config, Tree, NEAREST_NEIGHBORS_TO_CONSIDER);
        const auto nearest_neighbor_index = k_nearest_neighbors[0];
        auto was_epsilon_connection_made = epsilon_connect_RRTstar(nearest_neighbor_index, random_config, Tree, EPSILON,
                                                           numofDOFs, map, x_size, y_size, sample_count, goal_config,
                                                           DISCRETIZATION_FACTOR);

        if(was_epsilon_connection_made)
        {
            auto k_nearest_neighbors = get_k_nearest_neighbors(Tree[sample_count].c, Tree, NEAREST_NEIGHBORS_FOR_REWIRING+1);
            k_nearest_neighbors.erase(k_nearest_neighbors.begin()); /// Because the closest element will be the vertex itself.

            for(int i=0;i<k_nearest_neighbors.size();i++)
                rewire_RRT_star(sample_count,k_nearest_neighbors[i],Tree,map,x_size,y_size,numofDOFs);
            for(int i=0;i<k_nearest_neighbors.size();i++)
                rewire_RRT_star(k_nearest_neighbors[i],sample_count,Tree,map,x_size,y_size,numofDOFs);

            if(is_new_sample_in_goal_region(Tree[sample_count].c,goal_config,GOAL_REGION_THRESHOLD))
            {
                cout<<"Connection to goal region made"<<endl;
                auto distance_bw_configs = goal_config.get_distance(Tree[sample_count].c);
                Tree[sample_count+1] = RRT_Star_Node(goal_config,sample_count,Tree[sample_count].gcost + std::move(distance_bw_configs));  //Adding goal configuration to my tree
                is_goal_reached = true;
                sample_count++;     //Adding this for goal. So that backtrack starts at goal.
                break;
            }
        }

        sample_count++;
    }

     if(Tree[sample_count].c == goal_config)
         cout<<"Goal has been correctly reached"<<endl;

     if(is_goal_reached)
     {
         auto path_vector = get_path_vector(Tree,sample_count,numofDOFs,start_config);
         reverse(path_vector.begin(), path_vector.end());
         return change_path_vector_to_path(plan,path_vector,numofDOFs);
     }

     cout<<"Could not reach goal. Sample more points"<<endl;
     return -1;
}

//======================================================================================================================

static void planner(
		   double*	map,
		   int x_size,
 		   int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
	   int numofDOFs,
	   double*** plan,
	   int* planlength)
{
	//no plan by default
	*plan = NULL;
	*planlength = 0;
    int PRM_NUM_SAMPLES = 10;
    int RRT_NUM_SAMPLES = 60000;
    int RRT_CONNECT_NUM_SAMPLES = 10000;
    int RRT_STAR_NUM_SAMPLES = 100000;
    //auto num_samples = PRM_planner(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan,PRM_NUM_SAMPLES);
    auto num_samples = RRT_planner(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan,RRT_NUM_SAMPLES);
    //auto num_samples = RRT_connect(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan,RRT_NUM_SAMPLES);
    //auto num_samples = RRT_star(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan,RRT_STAR_NUM_SAMPLES);
    //num_samples = interpolation_based_plan(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan);
    *planlength = num_samples;
    return;
}

//======================================================================================================================

//prhs contains input parameters (3): 
//1st is matrix with all the obstacles
//2nd is a row vector of start angles for the arm 
//3nd is a row vector of goal angles for the arm 
//plhs should contain output parameters (2): 
//1st is a 2D matrix plan when each plan[i][j] is the value of jth angle at the ith step of the plan
//(there are D DoF of the arm (that is, D angles). So, j can take values from 0 to D-1
//2nd is planlength (int)
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[])
     
{ 
    
    /* Check for proper number of arguments */    
    if (nrhs != 4) { 
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                "Four input arguments required."); 
    } else if (nlhs != 2) {
	    mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required."); 
    } 
        
    /* get the dimensions of the map and the map matrix itself*/     
    int x_size = (int) mxGetM(MAP_IN);
    int y_size = (int) mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);
    
    /* get the start and goal angles*/     
    int numofDOFs = (int) (MAX(mxGetM(ARMSTART_IN), mxGetN(ARMSTART_IN)));
    if(numofDOFs <= 1){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "it should be at least 2");         
    }
    double* armstart_anglesV_rad = mxGetPr(ARMSTART_IN);
    if (numofDOFs != MAX(mxGetM(ARMGOAL_IN), mxGetN(ARMGOAL_IN))){
        	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "numofDOFs in startangles is different from goalangles");         
    }
    double* armgoal_anglesV_rad = mxGetPr(ARMGOAL_IN);
 
    //get the planner id
    int planner_id = (int)*mxGetPr(PLANNER_ID_IN);
    if(planner_id < 0 || planner_id > 3){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidplanner_id",
                "planner id should be between 0 and 3 inclusive");         
    }
    
    //call the planner
    double** plan = NULL;
    int planlength = 0;
    
    //you can may be call the corresponding planner function here
    //if (planner_id == RRT)
    //{
    //    plannerRRT(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    //}
    
    //dummy planner which only computes interpolated path
    planner(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength); 
    
    printf("planner returned plan of length=%d\n", planlength); 
    
    /* Create return values */
    if(planlength > 0)
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)planlength, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);        
        //copy the values
        int i,j;
        for(i = 0; i < planlength; i++)
        {
            for (j = 0; j < numofDOFs; j++)
            {
                plan_out[j*planlength + i] = plan[i][j];
            }
        }
    }
    else
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);
        //copy the values
        int j;
        for(j = 0; j < numofDOFs; j++)
        {
                plan_out[j] = armstart_anglesV_rad[j];
        }     
    }
    PLANLENGTH_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)1, mxINT8_CLASS, mxREAL); 
    int* planlength_out = (int*) mxGetPr(PLANLENGTH_OUT);
    *planlength_out = planlength;

    
    return;
    
}





