#include <iostream>
#include "Forest.h"
#include "point.h"
#include "geometryUtil.h"
#include <cmath>
#include <sstream>
using namespace std;

Forest::Forest(Configuration* config, Circuit* circuit):circuit(circuit), config(config){
    geometry::Point* initialState = new geometry::Point( circuit->getDimension(), circuit->getInitialState() ) ;
    kd = kd_create(circuit->getDimension());
    Node* root = new Node( initialState, true ) ;
    addRoot( root, forward ) ;
}

Forest::~Forest(){
    for(int i=0;i<nodes.size();i++)
        delete nodes[i];
     kd_free(kd);
}

Node* Forest::find_nearest_point(geometry::Point* q_sample){
    set = kd_nearest(kd, q_sample->getData());
    if(kd_res_size(set)>0){
        return (Node*) kd_res_item_data (set);
    }else{
        cout << "[error]" << endl ;
        return NULL;
    }
    //kd_res_free(set);
 
    /*
    int nearest_node= 0 ;
    //initializing nearest node with node 0 (root)
    
    double minimum_distance = nodes[nearest_node]->distance(q_sample, circuit);
    // O(n) loop to find nearest node, this loop can be optimized to O(log n) using KD-tree.
    for(int i=0;i<nodes.size();i++){
        double distance = nodes[i]->distance(q_sample, circuit);
        if(distance<minimum_distance){
            minimum_distance = distance;
            nearest_node= i;
        }
    }
    //cout << "nearest node is " << nearest_node << "  with distance = " << minimum_distance << endl ;
    return nodes[nearest_node];
     */
}


pair<geometry::Point*, Type> Forest::find_optimum_trajectory(geometry::Point* q_sample, Node* q_near){
    int shoots;                 config->getParameter("edu.uiuc.crhc.core.sampling.shoots", &shoots);
    double min_distance = 9999; config->getParameter("edu.uiuc.crhc.core.general.a-really-big-number", &min_distance);
    Type type ;
    geometry::Point* q_new = new geometry::Point(circuit->getDimension());
    bool mode;
    double coupling_random_v = random(0,1);
    cout << q_near->getNodeType() << endl ;
    if( (q_near->getNodeType() == forward) ||
       ((q_near->getNodeType()== coupling) && (coupling_random_v < 0.5))){
        mode = true;    //forward integration
        type=forward;
    }else{
        mode = false;   //backward integration
        type=backward;
    }

    if(shoots>1) cout << "there is a bug here. check" << endl ;
    for(int i=0;i<shoots; i++){
        double* tmp;
        bool isNan=false;
        tmp = circuit->integration(q_near->getState(), config, mode, q_near->getTime());
        //Set the new point to maximum value if the trace goes outside the bounding box
        for(int i=0; i<circuit->getDimension();i++){
            if(isnan(tmp[i])){
                isNan = true;
                continue;
            }
            //if(tmp[i]>circuit->getMax(i)) tmp[i]=circuit->getMax(i);
            //if(tmp[i]<circuit->getMin(i)) tmp[i]=circuit->getMin(i);
        }
        
        double distance = 0;
        for(int i=0;i<circuit->getDimension();i++){
            double di = tmp[i] - q_near->getData(i) ;
            double di_normalized = di / (circuit->getMax(i)-circuit->getMin(i));
            distance +=  di_normalized*di_normalized;
        }
        distance = sqrt(distance);
       
        //if(distance < min_distance ){
        //    min_distance = distance;
            q_new->setData(tmp);
        //}
        delete tmp;
    }
    bool invalidData=false;
    for(int i=0;i<circuit->getDimension();i++){
        if(isnan(q_new->getData(i))){
            invalidData=true;
        }
        if(q_new->getData(i) > 10 || q_new->getData(i) < -10)
            invalidData=true;
    }
    
    if(invalidData) cout  << "   i***** " << endl ;
    if(invalidData) return make_pair(q_new, invalid);
    else return make_pair(q_new, type);
}

void Forest::collision(Node* q_new, Node* q_near){
    for(int i=0;i<getSize();i++){
        if(!nodes[i]->isRoot()){
            if( geometry::intersection(q_new->getPoint(), q_near->getPoint(), nodes[i]->getPoint(), nodes[i]->getParent()->getPoint() ) ){
                if(!(q_new == nodes[i] || q_new==nodes[i]->getParent() || q_near==nodes[i] || q_near==nodes[i]->getParent()) ){
                    if( ((nodes[i]->getTreeType()==forward)&&(q_new->getTreeType()==coupling)) ||
                       ((nodes[i]->getTreeType()==coupling) && (q_new->getTreeType()==forward ))){
                        //Collision occured, do something.
                    }
                    
                }else{
                    //do something else
                }
            }
        }
    }
}

Node* Forest::addNode(geometry::Point* q_new, Node* parent, Type type){
    Node* new_node = new Node( q_new, false /*false means that this new node is NOT root */ );
    double dt=0;
    if(type==forward){
        config->getParameter("edu.uiuc.crhc.core.simulation.step", &dt);
    }else{
        config->getParameter("edu.uiuc.crhc.core.simulation.backwardstep", &dt);
        dt=-dt;
    }
    new_node->setTime( parent->getTime() + dt);
    new_node->setNodeType(type);
    new_node->setTreeType(parent->getTreeType());
    new_node->setParent(parent);
    nodes.push_back(new_node);
    kd_insert(kd, new_node->getState(), new_node);
    //collision(new_node, parent);
    return new_node;
}

// in case that the distance between the sampled node and nearest node is too much,
// we plant a new rrt instead of growing the nearest one.
Node* Forest::addRoot(geometry::Point* q, Type type){
    Node* new_node = new Node( q, true );
    new_node->setTime(0);
    new_node->setNodeType(type);
    new_node->setTreeType(type);
    nodes.push_back(new_node);
    kd_insert(kd, new_node->getState(), new_node);
    return new_node;
}

Node* Forest::addRoot(Node* new_node, Type type){
    new_node->setTime(0);
    new_node->setNodeType(type);
    new_node->setTreeType(type);
    nodes.push_back(new_node);
    kd_insert(kd, new_node->getState(), new_node);
    return new_node;
}


void Forest::load(string fileName){
    cout << "implement me! [Forest::load]" << endl ;
}

void Forest::save(string fileName){
    cout << "implement me! [Forest::save]" << endl ;
}

string Forest::toString(){
    cout << "implement me! [Forest::toString]" << endl ;
    return "";
}

int Forest::getDimension(){
    return d;
}

string Forest::draw(){
    return draw(0, 1);
}

string Forest::draw(int d0, int d1){
    
    stringstream str ;
    
    double xmin=9999, xmax=-99999, ymin=9999, ymax=-9999;
     for(int i=0;i<nodes.size();i++){
        if(nodes.at(i)->getState()[d0] > xmax) xmax=nodes.at(i)->getState()[d0];
        if(nodes.at(i)->getState()[d0] < xmin) xmin=nodes.at(i)->getState()[d0];
        if(nodes.at(i)->getState()[d1] > ymax) ymax=nodes.at(i)->getState()[d1];
        if(nodes.at(i)->getState()[d1] < ymin) ymin=nodes.at(i)->getState()[d1];
    }
    xmin = -7;
    xmax = 8;
    ymin = -0.2;
    ymax= 0.2;
    str << "plot [ " << xmin << ":" << xmax << "][" << ymin << ":" << ymax << "] 0 with linespoints lt \"white\" pt 0.01" ;
    str  << " title \"" << " " << "\"  \n";
    //cmdstr << "set ylabel " << "y(" << index << " )"  << " \n";
    //cmdstr << "set xlabel " << "time" << " \n";

    //cout << "drawing rrt " << nodes.size() << endl ;
    for(int i=0;i<nodes.size();i++){
        if( !nodes.at(i)->isRoot() ){
            double iToX = nodes.at(i)->getState()[d0] ;
            double iToY = nodes.at(i)->getState()[d1] ;
            double iFromX = nodes.at(i)->getParent()->getState()[d0] ;
            double iFromY = nodes.at(i)->getParent()->getState()[d1] ;
           // iFromX = iToX+0.001 ;
           // iFromY = iToY+0.001 ;
            //For debugging purposes: to check if we have an unusually long traces
            if(sqrt( ((iFromX-iToX)*(iFromX-iToX)) + ((iFromY-iToY)*(iFromY-iToY))) > 0.05){
                cout << "oops" << endl ;
                cout << iFromX << " " << iFromY << " to " << iToX << " " << iToY << endl ;
                cout << nodes.at(i)->isRoot() << " " << nodes.at(i)->getParent()->isRoot() << endl ;
            }
            if(nodes.at(i)->getNodeType()==forward){
                str << " set arrow from " << iFromX << "," << iFromY <<  "   to     " << iToX << "," << iToY  << "  nohead  lc rgb \"blue\" lw 2 \n" ;
            }else{
                 cout << "Hooray!=========================================" << endl ;
                cout << " set arrow from " << iFromX << "," << iFromY <<  "   to     " << iToX << "," << iToY  << "  nohead  lc rgb \"red\" lw 2 \n" ;

                str << " set arrow from " << iFromX << "," << iFromY <<  "   to     " << iToX << "," << iToY  << "  nohead  lc rgb \"red\" lw 2 \n" ;
            }
            //cout << str.str() << endl ;
        }
    }
    
    return str.str();
}


// Draw the signal index w.r.t. time
string Forest::draw(int index){
    stringstream str ;
    
    //scaling
    double min=9999, max=-99999, tmax=0;
    for(int i=0;i<nodes.size();i++){
        if(nodes.at(i)->getState()[index] > max) max=nodes.at(i)->getState()[index];
        if(nodes.at(i)->getState()[index] < min) min=nodes.at(i)->getState()[index];
        if(nodes.at(i)->getTime()>tmax) tmax = nodes.at(i)->getTime();
    }
    if( tmax==0 ){
        cout << "error: this is NOT a time-augmented RRT, so time-plotting is not allowed." << endl ;
        return "";
    }
    
    str << "plot [ " << 0 << ":" << tmax << "][" << min << ":" << max << "] 0 with linespoints lt \"white\" pt 0.01" ;
    str  << " title \"" << " " << "\"  \n";
    //cmdstr << "set ylabel " << "y(" << index << " )"  << " \n";
    //cmdstr << "set xlabel " << "time" << " \n";
    
    //cout << "drawing rrt " << nodes.size() << endl ;
    for(int i=0;i<nodes.size();i++){
        if( !nodes.at(i)->isRoot() ){
            double y_parent = nodes.at(i)->getParent()->getState()[index];
            double t_parent = nodes.at(i)->getParent()->getTime();
            double y_node = nodes.at(i)->getState()[index] ;
            double t_node = nodes.at(i)->getTime();
            if(nodes.at(i)->getNodeType()==forward){
                str << " set arrow from " << t_parent << "," << y_parent <<  "   to     " << t_node << "," << y_node  << "  nohead  lc rgb \"blue\" lw 2 \n" ;
            }else{
                cout << "Hooray!=========================================" << endl ;
                str << " set arrow from " << t_parent << "," << y_parent <<  "   to     " << t_node << "," << y_node  << "  nohead  lc rgb \"red\" lw 2 \n" ;
            }
            //cout << str.str() << endl ;
        }
    }
    
    return str.str();
}



int Forest::getSize(){
    return (int)(nodes.size());
}

Node* Forest::getNode(int i){
    return nodes.at(i);
}