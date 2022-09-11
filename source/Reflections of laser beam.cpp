/* 
Author: Aman Patel
Last Modified Date: 9/24/2021

Description: Takes in an initial reflections horizontal location along the AB segment
             of the mirror triangle, and calculates the number of reflections until the vector
             exits the system.
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>

double dot(std::pair<double, double> p1, std::pair<double, double> p2);
std::pair<double, double> turnVector(std::pair<double, double> p1, std::pair<double, double> p2);
std::pair<double, double> mulVector(double constant, std::pair<double, double> p2);
std::pair<double, double> findIntersection(std::pair<double, double> p1, std::pair<double, double> p2);
std::pair<double, double> sumVector(std::pair<double, double> p1, std::pair<double, double> p2);
std::pair<double, double> findNormal(std::pair<double, double> p1);

std::vector<std::pair<double,double >> Normals;
std::vector<int> lineSegment;  // ab = 0 , bc = 1, ca = 2
std::pair<double, double> A;
std::pair<double, double> B;
std::pair<double, double> C;
const double tolerance = 0.00000001;

std::vector<double> getT(std::pair<double,double> PriviousIntersect, std::pair<double,double> R){
    std::vector<double> T;
    //AB
    double t0 = ((10* sqrt(3) - PriviousIntersect.second)/R.second);
    T.push_back(t0);
    //BC
    double t1 = ((PriviousIntersect.second - sqrt(3)*PriviousIntersect.first)/(sqrt(3)*R.first - R.second));
    T.push_back(t1);
    double t2 = -((PriviousIntersect.second + sqrt(3)*PriviousIntersect.first)/(sqrt(3)*R.first + R.second));
    //AC
    T.push_back(t2);
    return T;
}

std::pair<double, double> intersection(std::pair<double,double> intersect, std::pair<double, double> reflector, std::vector<double> T){
    std::vector<std::pair<double,double >> vector;
    vector.emplace_back((intersect.first + T[0] * reflector.first), (intersect.second + T[0] * reflector.second));
    vector.emplace_back((intersect.first + T[1] * reflector.first), (intersect.second + T[1] * reflector.second));
    vector.emplace_back((intersect.first + T[2] * reflector.first), (intersect.second + T[2] * reflector.second));
    std::pair<double,double> point = vector[0];
    // determining if the intesection points are within tolerance
    for(int i = 0 ; i< vector.size(); i++){
        if((vector[i].first <= (intersect.first+tolerance) && vector[i].first >= (intersect.first- tolerance)) && (vector[i].second <= (intersect.second+ tolerance) && vector[i].second >= (intersect.second- tolerance))){
            continue;
        }else{
            if((vector[i].first >= -10 && vector[i].first <= 10) && (vector[i].second <= (10 * sqrt(3)) && vector[i].second >= 0)){
                return vector[i];
            }
        }
    }
    return std::make_pair(0,0);
}
// main function taking in x as an argument representing the initial beam's x location along the AB Wall segment
int main(int argc, char* argv[]) {

    //std::cout << "Please enter the starting angle:" << std::endl;
    //std::cin >> x;
    if (argc < 1) {
        std::cout << "Invalid amount of inputs" << std::endl;
    }
    const unsigned long x = std::stof(argv[1]);

    std::pair<double,double> intersect  = {x, 10 * sqrt(3)};
    std::pair<double, double> PrevIntersect(0,0);
    std::pair<double, double> V(0,0);
    std::pair<double, double> N(0,0);
    std::pair<double, double> R(0,0);
    std::pair<double, double> newLinePoint(0,0);
    lineSegment.push_back(-1);
    lineSegment.push_back(0);
    lineSegment.push_back(0);
    A = std::make_pair(-10,10 * sqrt(3));
    B = std::make_pair(10,10 * sqrt(3));
    C = std::make_pair(0,0);
    Normals.push_back(std::make_pair(0,-1));
    Normals.push_back(std::make_pair(-(sqrt(3)/2.0),0.5));
    Normals.push_back(std::make_pair((sqrt(3)/2.0),0.5));
    unsigned long total = 0;
    std::vector<double> T;
    do{
        // Calculating reflection vector using given formula
        total++;
        V = turnVector(intersect,PrevIntersect);
        N = findNormal(intersect);
        R = mulVector(-2*dot(N,V), N);
        R = sumVector(V, R);
        PrevIntersect = intersect;
        T = getT(PrevIntersect, R);
        newLinePoint = intersection(PrevIntersect, R, T);
        intersect = findIntersection(newLinePoint,PrevIntersect);
    }while(intersect.second>0.01);
    //Outputting total to text file output3.txt
    std::ofstream out("output3.txt");
    out << total;
    //printing total in console
    std::cout<< total<< std::endl;
    return 0;
}

std::pair<double, double> mulVector(double constant,std::pair<double, double> p2){
    return std::make_pair((constant* p2.first),(constant*p2.second));
}
std::pair<double, double> sumVector(std::pair<double, double> p1,std::pair<double, double> p2){
    return std::make_pair((p1.first + p2.first),(p1.second+p2.second));
}
std::pair<double, double> turnVector(std::pair<double, double> p1,std::pair<double, double> p2){
    return std::make_pair((p1.first - p2.first),(p1.second - p2.second));
}
double dot(std::pair<double, double> p1,std::pair<double, double> p2 ){
    return (p1.first*p2.first) + (p1.second*p2.second);
}

// find normal vector for previous intersection vertex
std::pair<double, double> findNormal(std::pair<double, double> p1){
    if(((p1.second >= ((10 * sqrt(3))- tolerance) && (p1.second <= ((10 * sqrt(3))+ tolerance) )&& p1.second >=0))){
        return Normals[0];
    }else if(p1.first >0){
        return Normals[1];
    }else{
        return Normals[2];
    }
}

// p1&p2 are a new point for the line segment and the current intersection point respectivley
std::pair<double, double> findIntersection(std::pair<double, double> p1,std::pair<double, double> p2)
{
    int currentLine;
    for(int i = 0; i<lineSegment.size(); i++){
        if(lineSegment[i] == -1){
            currentLine = i;
        }
    }
    double pointA = p1.second - p2.second;
    double pointB = p2.first - p1.first;
    double pointC = pointA * (p2.first) + pointB * (p2.second);
    double pointD = A.second - C.second;
    double pointE = C.first - A.first;
    double pointF = pointD*(C.first)+ pointE*(C.second);
    double pointG = C.second - B.second;
    double pointH = B.first - C.first;
    double pointI = pointG * (B.first) + pointH * (B.second);
    double pointJ = B.second - A.second;
    double pointK = A.first - B.first;
    double pointL = pointJ * (A.first) + pointK * (A.second);

    switch(currentLine) {
        case 0: {
            // we are reflecting of the AB segment
            // so new intersection can be from BC or AC
            // find the intersection point for of the line
            double determinant1 = pointA * pointE - pointD * pointB; // for CA
            double determinant2 = pointA * pointH - pointG * pointB; // for CB
            /// calculating for CA line segment
            if (determinant1 != 0) {
                double x = (pointE * pointC - pointB * pointF) / determinant1;
                double y = (pointA * pointF - pointD * pointC) / determinant1;
                if ((x >-10.1 && x<10.1) && (y <= (10 * sqrt(3)) && y >=0)) {
                    /// now we are reflecting of the vertex  AC
                    lineSegment[2] = -1;
                    /// dont forget to make AB to zero
                    lineSegment[currentLine] = 0;
                    return std::make_pair(x, y);
                }
            }
            /// for CB line intersection
            if (determinant2 != 0) {
                double x = (pointH * pointC - pointB * pointI) / determinant2;
                double y = (pointA * pointI - pointG * pointC) / determinant2;
                if ((x >=-10 && x<=10) && (y <= (10 * sqrt(3)) && y >=0))  {
                    ///now we are reflecting of the vertex CB
                    lineSegment[1] = -1;
                    lineSegment[currentLine] = 0;
                    return std::make_pair(x, y);
                }
            }
            break;
        }
        //  this is if we are reflecting of bc
        case 1:{
            double determinant1 = pointA * pointK - pointJ * pointB; // for AB
            double determinant2 = pointA * pointE - pointD * pointB; // for AC
            /// calculating for AB
            if (determinant1 != 0) {
                double x = (pointK * pointC - pointB * pointL) / determinant1;
                double y = (pointA * pointL - pointJ * pointC) / determinant1;
                if ((x >= -10 && x<=10) && ((y >= ((10 * sqrt(3))- tolerance) && (y <= ((10 * sqrt(3))+ tolerance) )&& y >=0)))  {
                    /// we are reflecting of the line AB
                    lineSegment[0] = -1;
                    lineSegment[currentLine] = 0;
                    return std::make_pair(x, y);
                }
            }
            ///calculating for AC
            if (determinant2 != 0) {
                double x = (pointE * pointC - pointB * pointF) / determinant2;
                double y = (pointA * pointF - pointD * pointC) / determinant2;
                if((x >= -10 && x<=0) && (y <= (10 * sqrt(3)) && y >=0))  {
                    /// we are reflecting of the line CB
                    lineSegment[2] = -1;
                    lineSegment[currentLine] = 0;
                    return std::make_pair(x, y);
                }
            }
            break;
        }
        // in case of AC for reflective surface
        case 2:{
            double determinant1 = pointA * pointH - pointG * pointB; // for BC
            double determinant2 = pointA * pointK - pointJ * pointB; // for AB
            ///calculating for bc
            if (determinant1 != 0) {
                double x = (pointH * pointC - pointB * pointI) / determinant1;
                double y = (pointA * pointI - pointG * pointC) / determinant1;
                if ((x >=-10 && x<=10) && (y <= (10 * sqrt(3)) && y >=0))  {
                    ///  we are now reflecting of the line AC
                    lineSegment[1] = -1;
                    lineSegment[currentLine] = 0;
                    return std::make_pair(x, y);
                }
            }
            /// calculating for AB intersection
            if (determinant2 != 0) {
                double x = (pointK * pointC - pointB * pointL) / determinant2;
                double y = (pointA * pointL - pointJ * pointC) / determinant2;
                if ((x >= -10 && x<=10) && ((y >= ((10 * sqrt(3))- tolerance) && (y <= ((10 * sqrt(3))+ tolerance) )&& y >=0)))  {
                    /// we are now reflecting of the line AB
                    lineSegment[0] = -1;
                    lineSegment[currentLine] = 0;
                    return std::make_pair(x, y);
                }
            }
            break;
        }
    }
    return std::make_pair(0, 0);
}
