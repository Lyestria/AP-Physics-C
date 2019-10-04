#include<bits/stdc++.h>
#include "Vector2D.cpp"
using namespace std;
const double eps=1e-9;
double m1=0.02,m2=0.02;
double r=0.05;
double I1=0.4*m1*r*r,I2=0.4*m2*r*r;
double u=0.5;
Vector2D p1,p2,p;
double L1,L2,L;
double E1,E2,E;
double w1,w2;
Vector2D s1(0,0.01),s2(0.1,0);
Vector2D v1(0.1,0),v2(-0.2,0);
double k=100;
const double dt=0.0001;
void step(){
    Vector2D a1(0,0),a2(0,0);
    double alpha1=0,alpha2=0;
    double dis=(s1-s2).mag();
    double PE1=0,PE2=0;
    if(dis<2*r){
        double x=(2*r-dis)*0.5;
        double Ff=u*k*x;
        a1=(s1-s2).norm()*k*x/m1;
        a2=(s2-s1).norm()*k*x/m2;
        PE1=PE2=0.5*k*x*x;
        Vector2D perp=(s1-s2).norm().perp();
        double ev1=v1.dot(perp)+r*w1;
        double ev2=v2.dot(perp)-r*w2;
        double v12=ev1-ev2;
        if(v12>eps){
            a1=a1+perp*Ff/m1;
            a2=a2-perp*Ff/m2;
            alpha1=(s1-s2).norm().cross(perp*Ff)*r/I1;
            alpha2=(s2-s1).norm().cross(perp*-Ff)*r/I2;
        }
        if(v12<-eps){
            a1=a1-perp*Ff/m1;
            a2=a2+perp*Ff/m2;
            alpha1=(s1-s2).norm().cross(perp*-Ff)*r/I1;
            alpha2=(s2-s1).norm().cross(perp*Ff)*r/I2;
        }
    }
    p1=v1*m1,p2=v2*m2,p=p1+p2;
    L1=s1.cross(v1)*m1+I1*w1;
    L2=s2.cross(v2)*m2+I2*w2;
    L=L1+L2;
    E1=0.5*m1*v1.dot(v1)+0.5*I1*w1*w1+PE1;
    E2=0.5*m2*v2.dot(v2)+0.5*I2*w2*w2+PE2;
    E=E1+E2;
    v1=v1+a1*dt;
    v2=v2+a2*dt;
    s1=s1+v1*dt;
    s2=s2+v2*dt;
    w1+=alpha1*dt;
    w2+=alpha2*dt;
}
int main(){
    freopen("//Users//kevinwan//Desktop//sim1.csv","w",stdout);
    for(int i=0;i<10000;i++){
        step();
        printf("%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",s1.x,s1.y,s2.x,s2.y,v1.x,v1.y,v2.x,v2.y,p1.x,p1.y,p2.x,p2.y,L1,L2,L,E1,E2,E);
    }
}
