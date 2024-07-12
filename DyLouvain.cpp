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
    void print()
    {
        printf("n:%d\n",n);
        printf("V:");for(auto u:V)printf("%d ",u);puts("");
        printf("E:");for(auto [u,v,w]:E)printf("[%d,%d,%d]",u,v,w);puts("");
        printf("eg:\n");
        for(int u=1;u<=n;u++)
        {
            printf("u:%d",u);
            for(auto [v,w]:eg[u])printf("[%d,%d]",v,w);puts("");
        }
        return;
    }
    void Read(FILE *f)
    {
        int u,v,w;
        while(true)
        {
            w=1;
            rin.rin(u,f);
            rin.rin(v,f);
            if(u==-1)break;

            u=(renum.find(u)==renum.end()?(renum[u]=++n):renum[u]);
            v=(renum.find(v)==renum.end()?(renum[v]=++n):renum[v]);

            E.push_back({u,v,w});m+=w;
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
    void Read_Inc(FILE *f)
    {
        int u,v,w;
        while(true)
        {
            w=1;
            rin.rin(u,f);
            rin.rin(v,f);
            if(u==-1)break;
            //暂时默认所有点都出现过
            u=(renum.find(u)==renum.end()?(renum[u]=++n):renum[u]);
            v=(renum.find(v)==renum.end()?(renum[v]=++n):renum[v]);
            Inc_E.push_back({u,v,w});m+=w;
        }
        for(auto [u,v,w]:Inc_E)
        {
            eg[u].push_back({v,w});
            eg_b[v].push_back({u,w});
            nd_deg_in[v]+=w;nd_deg_out[u]+=w;
            nd_deg_tot[v]+=w;nd_deg_tot[u]+=w;
            E2.push_back({u,v,w});
            E2.push_back({v,u,w});
        }
        return;
    }
    void insert(int u,int p)//u塞入社群p中
    {
        int lst=bel[u];
        if(lst!=0)
        {
            cm_deg_in[lst]-=nd_deg_in[u];
            cm_deg_out[lst]-=nd_deg_out[u];
            cm_deg_tot[lst]-=nd_deg_tot[u];
            for(auto [v,w]:eg[u])k_in[v][lst]-=w;
            for(auto [v,w]:eg_b[u])k_in[v][lst]-=w;
        }
        bel[u]=p;
        cm_deg_in[p]+=nd_deg_in[u];
        cm_deg_out[p]+=nd_deg_out[u];
        cm_deg_tot[p]+=nd_deg_tot[u];
        for(auto [v,w]:eg[u])k_in[v][p]+=w;
        for(auto [v,w]:eg_b[u])k_in[v][p]+=w;
        return;
    }
    double delQ(int u,int p,int k_in_1,int k_in_2)
    {
        double ans=k_in_1-k_in_2;
        ans+=(double(nd_deg_tot[u]))*(cm_deg_tot[bel[u]]-cm_deg_tot[p]-nd_deg_tot[u])/(m>>1);
        ans/=m;
        return ans;
    }
    void delQ_pr(int u,int p,int k_in_1,int k_in_2)
    {
        printf("del u:%d p:%d k_in_1:%d k_in_2:%d tot1:%d tot2:%d ki:%d\n",u,p,k_in_1,k_in_2,cm_deg_tot[p],cm_deg_tot[bel[u]],nd_deg_tot[u]);
    }
    double Modu()
    {
        map<int,int>Q;Q.clear();
        double s1=0,s2=0;
        for(auto [u,v,w]:E)
        {
            if(bel[u]==bel[v])
            {
                s1+=w;
            }
            // Q[bel[u]]+=w;
            // Q[bel[v]]+=w;
        }
        
        for(int i=1;i<=n;i++)
        {
            s2+=cm_deg_tot[i]*cm_deg_tot[i];
        }
        // printf("Modu s1:%.2lf s2:%.2lf m:%d\n",s1,s2,m);
        return (s1-s2/m)/m;
    }
    
    void rechecktop()
    {
        for(;true;)
        {
            auto [dQ,num]=q.top();
            auto [u,v,w]=E2[num];
            double nowQ=delQ(u,bel[v],k_in[u][bel[v]],k_in[u][bel[u]]);
            if(sgn(nowQ-dQ)!=0)
            {
                q.pop();
                q.push({nowQ,num});
                continue;
            }
            return;
        }
    }
    void Init_bel()
    {
        // int T=30;
        int T=3;
        vector<int>num;
        for(int i=0;i<=n;i++)num.push_back(i);
        for(;T-->0;)
        {
            shuffle(num.begin()+1,num.end(),engine);
            //maybe random shuffle
            for(int u=1;u<=n;u++)
            {
                int col=bel[u],cnt=0;
                // shuffle(eg[u].begin(),eg[u].end(),engine);
                // shuffle(eg_b[u].begin(),eg_b[u].end(),engine);
                unordered_map<int,int>mp;mp.clear();
                for(auto [v,w]:eg[u])if((mp[bel[v]]+=w)>cnt)cnt=mp[bel[v]],col=bel[v];
                for(auto [v,w]:eg_b[u])if((mp[bel[v]]+=w)>cnt)cnt=mp[bel[v]],col=bel[v];
                bel[u]=col;
            }
        }
        for(int u=1;u<=n;u++)
        {
            int p=bel[u];bel[u]=0;
            insert(u,p);
        }
        return;
    }
    void Combine()
    {
        puts("Start Combining");
        puts("End Combining");
        return;
    }
    void unfold()
    {
        printf("n:%d\n",n);
        double start=clock();
        bel.reserve(n+1);
        cm_deg_in.reserve(n+1);
        cm_deg_out.reserve(n+1);
        cm_deg_tot.reserve(n+1);
        for(int i=1;i<=n;i++)
        {
            bel[i]=i;
            cm_deg_in[i]=cm_deg_out[i]=cm_deg_tot[i]=0;
            // insert(i,i);
        }
        Init_bel();
        for(int i=0;i<E2.size();i++)
        {
            auto [u,v,w]=E2[i];
            q.push({delQ(u,bel[v],k_in[u][bel[v]],k_in[u][bel[u]]),i});
        }
        double lst_Modu=Modu();
        for(int rd=1;true;rd++)
        {
            rechecktop();
            auto [dQ,num]=q.top();
            if(dQ<revolution)break;
            auto [u,v,w]=E2[num];
            // printf("%d -> %d dQ:%.12lf\n",u,bel[v],dQ);
            insert(u,bel[v]);
        }
        puts("Start renum the community");
        renum.clear();//给社群离散化标号
        for(int i=1;i<=n;i++)
        {
            if(renum.find(bel[i])==renum.end())
            {
                renum[bel[i]]=++community_cnt;
            }
            bel[i]=renum[bel[i]];
        }
        group.reserve(community_cnt+1);
        vector<int>ept;ept.clear();
        for(int i=0;i<=community_cnt;i++)group.push_back(ept);
        for(int i=1;i<=n;i++)
        {
            group[bel[i]].push_back(i);
        }
        puts("End Unfolding");
        Combine();
        double end=clock();
        double t=double(end-start)/CLOCKS_PER_SEC;
        printf("time total :%.6lf\n",t);
        cost_time=t;
        return;
    }
    void print_community(FILE *f)
    {
        fprintf(f,"Cost Time  : %.10lf\n",cost_time);
        fprintf(f,"Modularity : %.10lf\n",Modu());
        for(int i=1,cntt=0;i<=community_cnt;i++,cntt=0)
        {
            fprintf(f,"Community %d:\n",i);
            for(auto u:group[i])
            {
                fprintf(f,"%5d ",u);
                if((++cntt)==20)
                {
                    cntt=0;
                    fprintf(f,"\n");
                }
            }
            fprintf(f,"\n");
        }
    }
};

class Louvain
{
    public:
        Graph G[1];
        void Load_G(string path)
        {
            puts("Start loading G");
            FILE *pub=fopen(path.c_str(),"r");
            G[0].clear();G[0].Read(pub);
            puts("End loading G");
            fclose(pub);
            return;
        }
        void Work(string path)
        {
            FILE *output=fopen(path.c_str(),"w");
            G[0].unfold();
            G[0].print_community(output);
            fclose(output);
        }
};

void facebook()
{
    Louvain louvain;
    louvain.Load_G("data/facebook_combined.txt");
    louvain.Work("output/facebook_sb.txt");
    // scp "D:\code\lab\DynaMo\data\facebook_combined.txt" zjj@211.87.232.180:/home/zjj/code/Dyno/data
    // scp zjj@211.87.232.180:/home/zjj/code/Dyno/facebook_sb.txt "D:\code\lab\DynaMo\data"
    return;
}
void twitter()
{
    Louvain louvain;
    louvain.Load_G("data/twitter_combined.txt");
    louvain.Work("output/twitter_sb.txt");
    // scp zjj@211.87.232.180:/home/zjj/code/Dyno/twitter_sb.txt "D:\code\lab\DynaMo\data"
    return;
}
void gplus()
{
    Louvain louvain;
    louvain.Load_G("data/gplus_combined.txt");
    louvain.Work("output/gplus_sb.txt");
    // scp zjj@211.87.232.180:/home/zjj/code/Dyno/gplus_sb.txt "D:\code\lab\DynaMo\data"
    return;
}
void Epinions()
{
    Louvain louvain;
    louvain.Load_G("data/soc-Epinions1.txt");
    louvain.Work("output/Epinions_sb.txt");
}
void LiveJournal()
{
    Louvain louvain;
    louvain.Load_G("data/soc-LiveJournal1.txt");
    louvain.Work("output/LiveJournal_sb.txt");
}
int main()
{
    facebook();
    twitter();
    gplus();
    Epinions();
    LiveJournal();
    return 0;
}
/*
https://sites.google.com/site/findcommunities/

tar -xzf gplus.tar.gz
gzip -d gplus_combined.txt.gz
g++ DyLouvain.cpp -o DyLv -std=c++17
scp "D:\code\lab\DynaMo\DyLouvain.cpp" zjj@211.87.232.180:/home/zjj/code/Dyno

# .tar.gz解压
tar -zxvf jdk-17_linux-aarch64_bin.tar.gz
# .tar.gz压缩
tar -czvf archive.tar.gz /path/to/directory

# .gz解压
gzip -d file.gz
# .gz压缩
gzip filename

# zip解压
unzip filename.zip
# zip压缩
zip archive.zip /path/to/file
*/