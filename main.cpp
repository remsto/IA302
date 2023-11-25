#include "ibex.h"
#include <iostream>
#include <vector>

using namespace std;
using namespace ibex;

struct Constraint{
  unsigned int id1;
  unsigned int id2;
  Interval distance;
};

vector<Constraint> readConstraintFromFile(string filename){
  ifstream file(filename.c_str());
  if(!file){
    cerr << "Error: cannot open file " << filename << endl;
    exit(1);
  }
  string line;
  vector<Constraint> constraints;
  while(getline(file, line)){
    // the first number is id1, the second number is id2, the third and fourth numbers are the lower and upper bounds of the distance
    unsigned int id1, id2;
    double lb, ub;
    istringstream iss(line);
    iss >> id1 >> id2 >> lb >> ub;
    cout << id1 << " " << id2 << " " << lb << " " << ub << endl;
    Constraint c;
    c.id1 = id1;
    c.id2 = id2;
    c.distance = Interval(-lb, ub); // Weird here: negative for lower bound 
    constraints.push_back(c);
  }
  file.close();
  return constraints;
}

bool checkTriangleInequality(vector<Constraint> &constraints){
  // Find the number of points
  unsigned int n = 0;
  for (unsigned int i=0; i<constraints.size(); i++){
    Constraint c = constraints[i];
    if (c.id1 > n){
      n = c.id1;
    }
    if (c.id2 > n){
      n = c.id2;
    }
  }
  vector<Constraint> neighbors[n]; // neighbors[i] is the list of constraints that involve point i
  for (unsigned int i=0; i<constraints.size(); i++){
    Constraint c = constraints[i];
    neighbors[c.id1-1].push_back(c);
    neighbors[c.id2-1].push_back(c);
  }
  // Find the triangles
  for (unsigned int i=0; i<n; i++){
    for (Constraint neighbor: neighbors[i]){
      unsigned int id_neighbor = (neighbor.id1-1 == i) ? neighbor.id2-1 : neighbor.id1-1;
      for (Constraint neighbor2: neighbors[id_neighbor]){
        unsigned int id_neighbor2 = (neighbor2.id1-1 == id_neighbor) ? neighbor2.id2-1 : neighbor2.id1-1;
        if (id_neighbor2 == i) // Avoid the same point
          continue;
        bool triangle_exists = false;
        for (Constraint neighbor3: neighbors[id_neighbor2]){
          unsigned int id_neighbor3 = (neighbor3.id1-1 == id_neighbor2) ? neighbor3.id2-1 : neighbor3.id1-1;
          if (id_neighbor3 == i){
            // Check the triangle inequality
            Interval d1 = neighbor.distance;
            Interval d2 = neighbor2.distance;
            Interval d3 = neighbor3.distance;
            if (d1.lb() > d2.ub()+d3.ub() || d2.lb() > d1.ub()+d3.ub() || d3.lb() > d1.ub()+d2.ub())
              return false;
            triangle_exists = true;
          }
        }
        if (!triangle_exists){
          // Add the triangle constraint
          Constraint triangle_constraint;
          triangle_constraint.id1 = i+1;
          triangle_constraint.id2 = id_neighbor2+1;
          triangle_constraint.distance = Interval(0, neighbor.distance.ub()+neighbor2.distance.ub());
          constraints.push_back(triangle_constraint);
        }
        }
      }
    }
  return true;
}
    

bool ConstraintIn(IntervalVector input_intervals, vector<Constraint> dist_constraints){
  Variable x1, y1, z1, x2, y2, z2;
  Function dist(x1, y1, z1, x2, y2, z2, sqrt(sqr(x1-x2)+sqr(y1-y2)+sqr(z1-z2)),"dist");
  for (unsigned int i=0; i<dist_constraints.size(); i++){
    Constraint c = dist_constraints[i];
    Interval xa = input_intervals[c.id1*3];
    Interval ya = input_intervals[c.id1*3+1];
    Interval za = input_intervals[c.id1*3+2];
    Interval xb = input_intervals[c.id2*3];
    Interval yb = input_intervals[c.id2*3+1];
    Interval zb = input_intervals[c.id2*3+2];
    Interval d = dist.eval(IntervalVector({xa, ya, za, xb, yb, zb}));
    if (!d.is_subset(c.distance)){
      return false;
    }
  }
  return true;
}

bool ConstraintOut(IntervalVector input_intervals, vector<Constraint> dist_constraints){
  Variable x1, y1, z1, x2, y2, z2;
  Function dist(x1, y1, z1, x2, y2, z2, sqrt(sqr(x1-x2)+sqr(y1-y2)+sqr(z1-z2)),"dist");
  for (unsigned int i=0; i<dist_constraints.size(); i++){
    Constraint c = dist_constraints[i];
    Interval xa = input_intervals[c.id1*3];
    Interval ya = input_intervals[c.id1*3+1];
    Interval za = input_intervals[c.id1*3+2];
    Interval xb = input_intervals[c.id2*3];
    Interval yb = input_intervals[c.id2*3+1];
    Interval zb = input_intervals[c.id2*3+2];
    Interval d = dist.eval(IntervalVector({xa, ya, za, xb, yb, zb}));
    if (d.is_disjoint(c.distance)){
      return true;
    }
  }
  return false;
}

void branchAndContract(IntervalVector start_interval, vector<Constraint> constraints, double tau=1e-4){
  stack<IntervalVector> stack_base;
  stack_base.push(start_interval);
  stack<IntervalVector> stack_acc;
  stack<IntervalVector> stack_rej;
  stack<IntervalVector> stack_unc;
  // Create distance function
  Variable x1, y1, z1, x2, y2, z2;
  Function dist(x1, y1, z1, x2, y2, z2, sqrt(sqr(x1-x2)+sqr(y1-y2)+sqr(z1-z2)),"dist");
  // Create the contractor
  Array<Ctc> contractors(0);
  for (auto constraint: constraints){
    contractors.add(*new CtcFwdBwd(*new NumConstraint(x1, y1, z1, x2, y2, z2, dist(x1, y1, z1, x2, y2, z2)=constraint.distance)));
  }
  CtcUnion ctc(contractors);
  cout << "Contractor created" << endl;
  while (stack_base.size()){
    IntervalVector x = stack_base.top();
    stack_base.pop();
    ctc.contract(x);
    if (ConstraintIn(x, constraints)){
      stack_acc.push(x);
    } else if (ConstraintOut(x, constraints)){
      stack_rej.push(x);
    } else if (x.max_diam() > tau){
      // Bisect x into two subintervals
      LargestFirst bb(0);
      pair<IntervalVector, IntervalVector> bisectResult = bb.bisect(x);
      stack_base.push(bisectResult.first);
      stack_base.push(bisectResult.second);
    }
    else {
      stack_unc.push(x);
    }
  }
  // Write the results to files
  ofstream file("results.txt");
  file << "Accepted intervals:" << endl;
  while (stack_acc.size()){
    IntervalVector x = stack_acc.top();
    stack_acc.pop();
    file << x << endl;
  }
  file << "Rejected intervals:" << endl;
  while (stack_rej.size()){
    IntervalVector x = stack_rej.top();
    stack_rej.pop();
    file << x << endl;
  }
  file << "Uncertain intervals:" << endl;
  while (stack_unc.size()){
    IntervalVector x = stack_unc.top();
    stack_unc.pop();
    file << x << endl;
  }
  file.close();
}

int main(int argc, char** argv) {
  vector<Constraint> constraints = readConstraintFromFile("dgsol-1.3/data/data_set_1/graph.01.data");
  cout << "Read " << constraints.size() << " constraints" << endl;
  cout << "Checking triangle inequality..." << endl;
  if (!checkTriangleInequality(constraints)){
    cout << "Triangle inequality is violated" << endl;
    exit(1);
  }
  cout << "Triangle inequality is satisfied" << endl;
  cout << "New constraints size: " << constraints.size() << endl;
  // Print all the constraints
  // for (unsigned int i=0; i<constraints.size(); i++){
  //   Constraint c = constraints[i];
  //   cout << c.id1 << " " << c.id2 << " " << c.distance << endl;
  // }
  // vector<Constraint> constraints = {Constraint({0, 1, Interval(1.2, 1.3)}), Constraint({0, 2, Interval(2.4, 2.5)}), Constraint({1, 2, Interval(1.6, 1.7)})};
  // double TAU = 0.01;
  // IntervalVector start_sol_interval = IntervalVector({Interval(0).inflate(0.001), Interval(0).inflate(0.001), Interval(1.25).inflate(0.001), Interval(0).inflate(0.001), Interval(1.976).inflate(0.001), Interval(1.448).inflate(0.001)});
  // IntervalVector start_interval = IntervalVector({Interval(0, 0.1), Interval(0, 0.1), Interval(1.225, 1.275), Interval(0, 0.1), Interval(-10, 10), Interval(-10, 10)});
  // IntervalVector big_start_interval = IntervalVector({Interval(0, 0.1), Interval(0, 0.1), Interval(-10, 10), Interval(-10, 10), Interval(-10, 10), Interval(-10, 10)});
  // branchAndContract(start_sol_interval, constraints, TAU);
}
