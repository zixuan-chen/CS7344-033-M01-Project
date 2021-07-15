#include <omp.h>
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <ctime>
using namespace std;

#define NODE_NUM 1024000
#define MAX_LINK 5
#define ITER 5
#define THTEADS 100

class Node{
    private:
        int id;
        double page_rank;
        vector<int> edges;
    public:
    Node(int _id = 0, double _weight = 1.0){
        id = _id;
        page_rank = _weight;
    }
    Node(const Node & n){
        id = n.id;
        page_rank = n.page_rank;
        edges.assign(n.edges.begin(), n.edges.end());
    }

    void linkTo(int nid){
        if (nid == id) return;
        if (find(edges.begin(), edges.end(), nid) != edges.end())
            return;
        edges.push_back(nid);
    }
    int getDegree(){return edges.size();}
    int getId(){return id;}
    int getPageRank(){return page_rank;}
    vector<int> & getAllLinks(){return edges;}
    void setId(int _id){id = _id;}
    void setPageRank(double _pgrk){page_rank = _pgrk;}

};

class Graph {
    public:
        vector<Node> nodes;
        Graph(){nodes.clear();}
        Graph(int n){
            nodes.resize(n);
            for (int i = 0; i < n; i++)
                nodes[i].setId(i);
        }
        int getNodeCnt(){
            return nodes.size();
        }
        int getDegree(int id){
            if(id >= 0 && id < nodes.size())
                return nodes[id].getDegree();
            else{
                cout << "invalid id in Graph::getDegree" << endl;
                exit(-1);
            }
        }
        void link(int n1, int n2){
            if (n1 == n2)
                return;
            nodes[n1].linkTo(n2);
            nodes[n2].linkTo(n1);
        }
        Node & getNode(int id){
            if(id >= 0 && id < nodes.size())
                return nodes[id];
            else{
                cout << "invalid id in Graph::getNode" << endl;
                exit(-1);
            }
        }
        void setPageRank(int id, double value){
            if(id >= 0 && id < nodes.size())
                nodes[id].setPageRank(value);
        }

        
};

class SubGraph: public Graph{
    map<int, int> lidOf; // map from global index -> local index
    vector<int> gidOf;
    map<int, vector<int>> outLinksOfNode;
    vector<double> page_rank_old;

    public:
    void addNewNode(const Node &n){
        nodes.push_back(n);
    }
    vector<int> & getAllLinks(int gid){
        return nodes[lidOf[gid]].getAllLinks();
    }
    int getDegree(int gid){
        return nodes[lidOf[gid]].getDegree();
    }
    double getPageRank(int gid){
        return nodes[lidOf[gid]].getPageRank();
    }
    void setPageRank(int gid, double value){
        nodes[lidOf[gid]].setPageRank(value);
    }
    void initialize(){
        gidOf.resize(nodes.size());
        page_rank_old.resize(nodes.size());
        for(int lid = 0; lid < nodes.size(); ++lid){
            lidOf[nodes[lid].getId()] = lid; // construct a map from global index to local index
            gidOf[lid] = nodes[lid].getId(); // construct a map from local index to global index
        }
        
        for(int lid = 0; lid < nodes.size(); ++lid){
            for(int gid : nodes[lid].getAllLinks()){
                if(lidOf.find(gid) == lidOf.end()){
                    outLinksOfNode[gidOf[lid]].push_back(gid);
                }
            }
        }
    }
    void updatePageRankLocal(){
        double page_rank_sum ;
        // lid: local index in sub graph 
        // gid: global index in father graph
        for (int lid = 0; lid < nodes.size(); ++lid){
            page_rank_old[lid] = nodes[lid].getPageRank();
        }
        for (int lid = 0; lid < nodes.size(); ++lid){
            page_rank_sum = 0;
            for(int gid : nodes[lid].getAllLinks()){
                
                if(lidOf.find(gid) != lidOf.end())
                    page_rank_sum += nodes[lidOf[gid]].getPageRank();
            }
            nodes[lid].setPageRank(page_rank_sum / nodes[lid].getDegree());
        }
    }
    vector<int> & getAllNodeGlobalIndex(){            
        return gidOf;
    }
    map<int, vector<int>> & getAllLinksOutMap(){
        return outLinksOfNode;
    }
    double getOldPageRank(int gid){
        return page_rank_old[lidOf[gid]];
    }
};


// void graphPartition(Graph &g, int num, vector<SubGraph> &sgs){
//     int numNodes = g.getNodeCnt();
//     int limit = numNodes / num + 1;
//     vector<int> searched(g.getNodeCnt(), -1);
//     vector<int> degrees(g.getNodeCnt(), 0);
//     list<int> buffer;
//     int root, maxDegree, cnt;
//     int i, j;

//     for (i = 0; i < numNodes; ++i)
//         degrees[i] = g.getDegree(i);

//     cnt = 0;
//     for (i = 0; i < num; ++i){
//         for (j = 0; j < limit && j+cnt < numNodes; ++j){
//             sgs[i].addNewNode(g.getNode(j+cnt));
//         }
//         cnt += j;
//         sgs[i].initialize();
//     }
// }

void graphPartition(Graph &g, int num, vector<SubGraph> &sgs){
    int numNodes = g.getNodeCnt();
    int limit = numNodes / num + 1;
    vector<int> searched(g.getNodeCnt(), -1);
    vector<int> degrees(g.getNodeCnt(), 0);
    list<int> buffer;
    int root, maxDegree, cnt;
    int i, j;

    for (i = 0; i < numNodes; ++i)
        degrees[i] = g.getDegree(i);

    for (i = 0; i < num; ++i){
        root = 0;

        // find root node with highest degree
        for(j = 0; j < numNodes; ++j){
            if (searched[j] ==-1 && degrees[j] >= degrees[root]){
                root = j;
            }
        }    

        // start from root, DFS find the rest nodes
        buffer.clear();
        buffer.push_back(root);
        cnt = 0;
        while(cnt <= limit && !buffer.empty()){
            root = buffer.front();
            sgs[i].addNewNode(g.getNode(root));
            ++cnt;
            searched[root] = i;
            buffer.pop_front();
            for(auto gid : g.getNode(root).getAllLinks()){
                if (searched[gid] == -1){
                    buffer.push_back(gid);
                }
            }
        }
        sgs[i].initialize();
    }
}

void graphMerge(Graph &g, int num, vector<SubGraph> &sgs){
    for(int i = 0; i < num; ++i){
        for(int gid : sgs[i].getAllNodeGlobalIndex()){
            g.setPageRank(gid, sgs[i].getPageRank(gid));
        }
    }
}

void syncGraph(vector<SubGraph> &sgs, int num){
    // construct a map from global index -> index of subgraph
    map<int, int> mapGid2SubG;
    vector<int> gids;
    double page_rank;
    int sgid, gid;

    for (int i = 0; i < num; ++i){
        gids = sgs[i].getAllNodeGlobalIndex();
        for(auto _gid: gids)
            mapGid2SubG[_gid] = i;
    }
    
    for (int i = 0; i < num; ++i){
        page_rank = 0;
        for (auto outlink: sgs[i].getAllLinksOutMap()){
            gid = outlink.first;
            for(int gid2 : outlink.second){
                sgid = mapGid2SubG[gid2];
                page_rank += sgs[sgid].getOldPageRank(gid2);
            }
            page_rank = page_rank / sgs[i].getDegree(gid) + sgs[i].getPageRank(gid);
            sgs[i].setPageRank(gid, page_rank);
        }
    }
    
}

void PageRank(Graph &g, int iterations, int threads){
    vector<SubGraph> sgs(threads);
    graphPartition(g, threads, sgs);
    // partition graphs into multiple sub-graphs, each process
    // manages one part syncGraph(sgs, threads);
    for (int t = 0; t < iterations; t++){
        cout << "iter :  " << t << " running" << endl;
        #pragma omp parallel
        {
            #pragma omp for
            for (int tid = 0; tid < threads; tid++){
                // do pagerank on the local sub-Graph
                sgs[tid].updatePageRankLocal();
            }
            syncGraph(sgs, threads);    //synchronize sub-grpahs
        }
    }
    graphMerge(g, threads, sgs);
}

int main(int argc, char const *argv[])
{
    printf("Initiallzing graph of node=%d...\n", NODE_NUM);
    Graph g(NODE_NUM);
    int numlink;
    int rand_id;
    clock_t startTime,endTime;

    printf("generate edge count of different nodes ranges from 1 to 10...\n");
    for (int id = 0; id < NODE_NUM; ++id)
    {
        numlink = (rand() % (MAX_LINK))+ 1;
        for(int j = 0; j < numlink; ++j)
        {
            rand_id = (rand() %(NODE_NUM));
            g.link(id, rand_id);
        }
    }
    printf("executing page rank\n");
    startTime = clock();
    PageRank(g, ITER, THTEADS);
    endTime = clock();
    printf("Total elapsed time: %10.6f\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);
    return 0;
}
