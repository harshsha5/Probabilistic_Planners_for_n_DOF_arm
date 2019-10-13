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
//#include <boost/functional/hash.hpp>

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
            dist+= pow((angles[i]-c2.angles[i]),2);
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

//struct config_hasher{
//    long double operator()(vector<double> const& vec) const
//    {
//        long double seed = vec.size();
//        for(auto& i : vec) {
//            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//        }
//        return seed;
//    }
//};

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
    const double MIN_ANGLE = -PI;
    const double MAX_ANGLE = PI;
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

bool is_connection_possible(unordered_map<int,Node> road_map,
                                  const configuration &start_config,
                                  const configuration &end_config,
                                  const int &numofDOFs,
                                  double*	map,
                                  const int &x_size,
                                  const int &y_size)
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
        if(is_connection_possible(road_map,road_map[neighbor].c,road_map[sample_count].c,numofDOFs,map,x_size,y_size))
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
            if(is_connection_possible(road_map,start_config,elt.second.c,numofDOFs,map,x_size,y_size))
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
        if(is_connection_possible(road_map,goal_config,elt.second.c,numofDOFs,map,x_size,y_size))
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
              const int &goal_index)
{
    visited.insert(start_index);
    s.push(start_index);

    if(start_index==goal_index)
        return;

    for(size_t i=0;i<road_map[start_index].neighbors.size();i++)
    {
        if(!visited.count(road_map[start_index].neighbors[i]))
            dfs_util(road_map,visited,s,road_map[start_index].neighbors[i],goal_index);
    }
    s.pop();
}

//======================================================================================================================

void get_plan_from_stack(unordered_map<int,Node> road_map,
              double*** plan,
              stack<int> s,
              const int &numofDOFs)
{
    s.pop();        //Removing the goal. As the plan only needs to include the start till goal-1
    stack<int> reversed_stack;
    while(!s.empty())
    {
        reversed_stack.push(s.top());
        s.pop();
    }

    int i=0;
    while(!reversed_stack.empty())
    {
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
    dfs_util(road_map,visited,s,start_index,goal_index);
    cout<<"Stack length: "<<s.size()<<endl;
    get_plan_from_stack(road_map,plan,s,numofDOFs);
    return s.size()-1;
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
  while(iteration_number<5)
          {
              cout<<"Iteration Number: "<<iteration_number<<endl;
              auto road_map = build_road_map(map,x_size,y_size,numofDOFs,num_samples_PRM);
              auto got_connected = add_start_and_goal_to_road_map(map,x_size,y_size,numofDOFs,road_map,start_config,goal_config);

              if(got_connected){
                  auto plan_length = search_road_map(road_map,plan,numofDOFs);
                  cout<<"Plan_Length: "<<plan_length<<endl;
                  return plan_length;
              }
              num_samples_PRM *=2;
              cout<<"New number of samples are: "<<num_samples_PRM<<endl;
              iteration_number++;
          }
  cout<<"PRM didn't find route"<<endl;
  return 5;
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
    auto num_samples = PRM_planner(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan,PRM_NUM_SAMPLES);
    //auto num_samples = interpolation_based_plan(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan);
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





