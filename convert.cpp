#include <set>
#include <map>
#include <cmath>
#include <array>
#include <queue>
#include <time.h>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include <cstdio>
#include <algorithm>
#include <unordered_map>
#define LL long long
using namespace std;
bool ifnum(char x){return x>='0'&&x<='9';}
bool ifupchr(char x){return x>='A'&&x<='Z';}
bool iflochr(char x){return x>='a'&&x<='z';}
bool ifchr(char x){return x=='#'||x==' '||x=='-'||x=='.';}
struct Rin
{
    char c;
    char gc(FILE *f){c=EOF;fscanf(f,"%c",&c);return c;}
    void rin(int &x,FILE *f)
    {
        x=-1;
        for(c=gc(f);c>'9'||c<'0';c=gc(f))if(c==EOF)return;x=0;
        for(;c>='0'&&c<='9';c=gc(f))x=(x<<1)+(x<<3)+(c^'0');
    }
    void rin(LL &x,FILE *f)
    {
        x=-1;
        for(c=gc(f);c>'9'||c<'0';c=gc(f))if(c==EOF)return;x=0;
        for(;c>='0'&&c<='9';c=gc(f))x=(x<<1)+(x<<3)+(c^'0');
    }
    void rin(char &x,FILE *f)
    {
        for(c=gc(f);!ifnum(c)&&!iflochr(c)&&!ifupchr(c)&&!ifchr(c);c=gc(f))if(c==EOF)return;
        x=c;
    }
    void rin(string &x,FILE *f)
    {
        x.clear();
        for(c=gc(f);!ifnum(c)&&!iflochr(c)&&!ifupchr(c)&&!ifchr(c);c=gc(f))if(c==EOF)return;
        for(;ifnum(c)||iflochr(c)||ifupchr(c)||ifchr(c);c=gc(f))x.push_back(c);
    }
}rin;
void jh(int &x,int &y){if(x^y)x^=y^=x^=y;return;}
void jh(LL &x,LL &y){if(x^y)x^=y^=x^=y;return;}
int min(int x,int y){return x<y?x:y;}
int max(int x,int y){return x>y?x:y;}
LL min(LL x,LL y){return x<y?x:y;}
LL max(LL x,LL y){return x>y?x:y;}



mt19937 engine(chrono::_V2::steady_clock::now().time_since_epoch().count());
inline int rand(int l,int r){return uniform_int_distribution<int>(l,r)(engine);}

// #define Directed
const double eps=1e-9;
const double revolution=2e-9;
bool sgn(double v){return (v>eps)-(v<-eps);}
struct Graph{
    struct edge
    {
        int u,v,w;
        double delQ;
    };
    double cost_time;
    int n,m;
    vector<int>V;
    vector<array<int,3> >E;
    vector<array<int,3> >E2;//正反方向都存
    vector<array<int,3> >Inc_E;
    vector<vector<pair<int,int> > >eg;
    vector<vector<pair<int,int> > >eg_b;

    vector<int>nd_deg_in;
    vector<int>nd_deg_out;
    vector<int>nd_deg_tot;
    vector<int>cm_deg_in;
    vector<int>cm_deg_out;
    vector<int>cm_deg_tot;

    int community_cnt;
    vector<int>bel;//归属社群编号
    vector<vector<int> >group;
    double norm=double(1.0);

    map<int,int>renum;//离散化
    vector<int>backnum;
    
    vector<map<int,int> >k_in;
    priority_queue<pair<double,int>,vector<pair<double,int> >,less<pair<double,int> > >q;
    void clear()
    {
        n=0;m=0;community_cnt=0;
        V.clear();E.clear();E2.clear();
        eg.clear();
        nd_deg_in.clear();nd_deg_out.clear();nd_deg_tot.clear();
        cm_deg_in.clear();cm_deg_out.clear();cm_deg_tot.clear();
        bel.clear();group.clear();renum.clear();backnum.clear();
        k_in.clear();
        for(;!q.empty();q.pop());
        return;
    }
    void Read(FILE *fin,FILE *fout)
    {
        int u,v,w;
        while(true)
        {
            w=1;
            rin.rin(u,fin);
            rin.rin(v,fin);
            if(u==-1)break;

            // printf("%d %d\n",u,v);
            u=(renum.find(u)==renum.end()?(renum[u]=++n):renum[u]);
            v=(renum.find(v)==renum.end()?(renum[v]=++n):renum[v]);

            E.push_back({u,v,w});m+=w;
            // printf("%d %d\n",u,v);
            fprintf(fout,"%d %d\n",u,v);
        }
        backnum.reserve(n+1);
        nd_deg_in.reserve(n+1);
        nd_deg_out.reserve(n+1);
        nd_deg_tot.reserve(n+1);
        for(auto [u,num]:renum)backnum[num]=u;
        
        map<int,int>qpt;qpt.clear();
        vector<pair<int,int> >ept;ept.clear();
        for(int i=0;i<=n;i++)
        {
            k_in.push_back(qpt);
            eg.push_back(ept);
            eg_b.push_back(ept);
            nd_deg_in[i]=nd_deg_out[i]=nd_deg_tot[i]=0;
        }
        for(auto [u,v,w]:E)
        {
            eg[u].push_back({v,w});
            eg_b[v].push_back({u,w});
            nd_deg_in[v]+=w;
            nd_deg_out[u]+=w;
            nd_deg_tot[v]+=w;
            nd_deg_tot[u]+=w;
            E2.push_back({u,v,w});
            E2.push_back({v,u,w});
        }
        return;
    }
};

class Louvain
{
    public:
        Graph G[1];
        void Load_G(string path,string path2)
        {
            puts("Start loading G");
            FILE *pub=fopen(path.c_str(),"r");
            FILE *out=fopen(path2.c_str(),"w");
            G[0].clear();
            G[0].Read(pub,out);
            puts("End loading G");
            fclose(pub);
            fclose(out);
            return;
        }
};

// void facebook()
// {
//     Louvain louvain;
//     louvain.Load_G("data/facebook_combined.txt","data");
//     return;
// }
void twitter()
{
    Louvain louvain;
    louvain.Load_G("data/twitter_combined.txt","data/twitter_renumber.txt");
    return;
}
void gplus()
{
    Louvain louvain;
    louvain.Load_G("data/gplus_combined.txt","data/gplus_renumber.txt");
    return;
}
// void Epinions()
// {
//     Louvain louvain;
//     louvain.Load_G("data/soc-Epinions1.txt");
//     louvain.Work("output/Epinions_sb.txt");
// }
// void LiveJournal()
// {
//     Louvain louvain;
//     louvain.Load_G("data/soc-LiveJournal1.txt");
//     louvain.Work("output/LiveJournal_sb.txt");
// }
int main()
{
    twitter();
    gplus();
    return 0;
}