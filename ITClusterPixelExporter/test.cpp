


#include <iostream>
#include <map>


int main(int argc, char const *argv[])
{
    
    std::map<int,int> m;


    auto l = m.find(5);
 
    auto a = *l;

    return 0;
}
