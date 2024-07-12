#include <set>
#include <map>
#include <cmath>
#include <array>
#include <queue>
#include <time.h>
#include <vector>
#include <string>
#include <cstdio>
#include <stdlib.h>
#include <algorithm>
#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
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
        x=0;
        for(c=gc(f);c>'9'||c<'0';c=gc(f))if(c==EOF)return;
        for(;c>='0'&&c<='9';c=gc(f))x=(x<<1)+(x<<3)+(c^'0');
    }
    void rin(LL &x,FILE *f)
    {
        x=0;
        for(c=gc(f);c>'9'||c<'0';c=gc(f))if(c==EOF)return;
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

string name;
string st1="curl https://snap.stanford.edu/data/";
string st2=" --output ";
string tst="cd ..";
int main()
{
    FILE *f=fopen("download.in","r");
    for(rin.rin(name,f);name.size()>0;rin.rin(name,f))
    {
        system((st1+name+st2+name).c_str());
        // system(tst.c_str());
    }
}
/*
gcc download.cpp -o download

*/