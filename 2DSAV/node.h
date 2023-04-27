#ifndef NODE_H
#define NODE_H

#include "vec3.h"
#include <math.h>
#include <iostream>

class Node{
 public:
  double x;
  double y;
  double z;
  double dx;
  double dy;
  double dz;

  Node(){
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
    this->dx = 0.0;
    this->dy = 0.0;
    this->dz = 0.0;
  }

  Node(double x, double y, double z, double dx, double dy, double dz){
    this->x = x;
    this->y = y;
    this->z = z;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
  }

  Node(vec3 pos, vec3 dir){
    this->x = pos.x;
    this->y = pos.y;
    this->z = pos.z;
    this->dx = dir.x;
    this->dy = dir.y;
    this->dz = dir.z;
  }

  double getLength(Node node1, Node node2){
    return sqrt(pow(node1.x-node2.x,2)+pow(node1.y-node2.y,2)
		+pow(node1.z-node2.z,2));
  }

  void printInfo(){
    std::cout << "The node is at (" 
	      << x << ", " << y << ", " << z <<")" << std::endl;
    std::cout << "The direction is (" 
	      << dx << ", " << dy << ", " << dz <<")" << std::endl;
  }
};

#endif
