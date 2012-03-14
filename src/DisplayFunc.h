#ifndef DISPLAYFUNC_H
#define DISPLAYFUNC_H

#include <iostream>

// Define macros
#define output(x) std::cout << (x) << std::endl ;
#define output2(x,y) std::cout << (x) << (y) << std::endl ;
#define output3(x,y,z) std::cout << "x = "<< (x) << " y = "<< (y) << " z = "<< (z) << std::endl ;

int intlen(int start) {
    int end = 0;
    while(start > 0) {
        start = start/10;
        end++;
    }
    return end;
}

void display(int count){
    int len,j;
    cout<<count;
    len = intlen(count);
    for(j = 1; j<=len; j++){
        cout<<"\b";
    }
}
template <typename T>
void displayPoint(T p, const char* prefix = ""){
    cout<<prefix<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;
}

void displayVect(geometry::Vector v, const char* prefix = ""){
    cout<<prefix<<v.getX()<<" "<<v.getY()<<" "<<v.getZ()<<endl;
}

void displayStdVect(std::vector <double> v, const char* prefix = ""){
    cout<<"size = "<<v.size()<<endl;
    for(unsigned int i = 0; i < v.size(); i++){
        cout<<i<<" = "<<v[i]<<endl;
    }
}

void displayVector(itk::VariableLengthVector <double> v){
    cout<<"size = "<<v.GetSize()<<endl;
    for(unsigned int i = 0; i < v.GetSize(); i++){
        cout<<i<<" = "<<v.GetElement(i)<<endl;
    }
}

void displayVector(std::vector <double> v){
    cout<<"size = "<<v.size()<<endl;
    for(unsigned int i = 0; i < v.size(); i++){
        cout<<i<<" = "<<v[i]<<endl;
    }
}

void displayFiberSegment(FiberSegment f){
    cout<<"segment.p1 = "<<f.p1<<", segment.p2 = "<<f.p2<<endl;
}

#endif
