#include <iostream>
#include <vector>
#include <algorithm>4?

struct node {
        int id;
        std::vector<int> neighbs;
};

class graph {
    public:
        std::vector<node> nodes;
        
        int Num_Nodes() {
            return nodes.size();
        }


};

int main() {
    node a;
    node b;
    a.id = 0;
    b.id = 1;
    a.neighbs.push_back (b.id);
    b.neighbs.push_back (a.id);
    for (int i: b.neighbs)
        std::cout << i << '\n';
    graph G1;
    G1.nodes.push_back (a); 
    G1.nodes.push_back (b);
    //int i = G1.Has_Edge(G1.nodes); 
    //std::cout << i << "\n";
}