/*
  Problem: https://codeforces.com/gym/101673
  Explication:
*/

#include <bits/stdc++.h>
using namespace std;

#define ll long long int
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()

const int Inf = 1e9;
const ll INF = 1e18;
const int maxn = 1e5;

using ld = long double;
const ld eps = 1e-9, pi = acos(-1);
bool geq(ld a, ld b){return a - b >= -eps;}
bool leq(ld a, ld b){return b - a >= -eps;}
bool ge(ld a, ld b){return a - b > eps;}
bool le(ld a, ld b){return b - a > eps;}
bool eq(ld a, ld b){return abs(a - b) <= eps;}
bool neq(ld a, ld b){return abs(a - b) > -eps;}

struct point{
    ld x, y;
    point(): x(0), y(0){}
    point(ld x, ld y): x(x), y(y){}
    
    point operator+(const point & p) const{return point(x + p.x, y + p.y);}
    point operator-(const point & p) const{return point(x - p.x, y - p.y);}
    point operator*(const ld & p) const{return point(x * p, y * p);}
    point operator/(const ld & p) const{return point(x / p, y / p);}

    point operator+=(const point & p){*this = *this + p; return *this;}
    point operator-=(const point & p){*this = *this - p; return *this;}
    point operator*=(const ld & p){*this = *this * p; return *this;}
    point operator/=(const ld & p){*this = *this / p; return *this;}

    point rotate(const ld & a) const{return point(x * cos(a) - y * sin(a), x * sin(a) + y * cos(a));}
    point perp() const{return point(-1 * y, x);}
    ld dot(const point & p) const{return x * p.x + y * p.y;}
    ld cross(const point & p) const{return x * p.y - y * p.x;}
    ld norm() const{return x * x + y * y;}
    ld length() const{return sqrtl(x * x + y * y);}
    point unit() const{return (*this) / length();}

    bool operator==(const point & p) const{return eq(x, p.x) && eq(y, p.y);}
    bool operator!=(const point & p) const{return !(*this == p);}
    bool operator<(const point & p) const{return le(x, p.x) || (eq(x, p.x) && le(y, p.y));}
    bool operator>(const point & p) const{return ge(x, p.x) || (eq(x, p.x) && ge(y, p.y));}

};

ostream &operator<<(ostream &os, const point & p){return os << "("<< p.x<<", "<<p.y<<")";}

vector<vector<point>> tangents(const point & c1, ld r1, const point & c2, ld r2, bool inner){
	//returns a vector of segments or a single point
	if(inner) r2 = -r2;
	point d = c2 - c1;
	ld dr = r1 - r2, d2 = d.norm(), h2 = d2 - dr*dr;
	if(eq(d2, 0) || le(h2, 0)) return {};
	point v = d*dr/d2;
	if(eq(h2, 0)) return {{c1 + v*r1}};
	else{
		point u = d.perp()*sqrt(h2)/d2;
		return {{c1 + (v - u)*r1, c2 + (v - u)*r2}, {c1 + (v + u)*r1, c2 + (v + u)*r2}};
	}
}

bool cmp(pair<point, pair<point, ld>> a, pair<point, pair<point, ld>> b){
    return a.fi < b.fi;
}

bool cmp1(pair<point, ld> a, pair<point, ld> b){
    return a.se > b.se;
}

vector<pair<point, pair<point, ld>>> convexHull(vector<pair<point, pair<point, ld>>> P){
        int n = P.size();
        sort(all(P), cmp);
        vector<pair<point, pair<point, ld>>> L, U;
        for(int i = 0; i < n; i++){
            while(L.size() >= 2 && leq((L[L.size() - 2].fi - P[i].fi).cross(L[L.size() - 1].fi - P[i].fi), 0)){
                L.pop_back();
            }
            L.pb(P[i]);
        }
        for(int i = n - 1; i >= 0; i--){
            while(U.size() >= 2 && leq((U[U.size() - 2].fi - P[i].fi).cross(U[U.size() - 1].fi - P[i].fi), 0)){
                U.pop_back();
            }
            U.pb(P[i]);
        }
        L.pop_back(); U.pop_back();
        L.insert(L.end(), U.begin(), U.end());
        return L;
    }

ld segmentCircular(point c, point a, point b){
    ld r = (c - a).length();
    point v1 = a - c, v2 = b - c;
    ld tet = acosl(v1.dot(v2) / v1.length() / v2.length());
    point rt = v1.rotate(tet);
    if(!(abs(rt.x - v2.x) <= .5 && abs(rt.y - v2.y) <= .5)) tet = pi * 2 - tet;
    return (a - c).length() * tet;
}

vector<pair<point, ld>> clear(vector<pair<point, ld>> P){
    int n = P.size();
    int buck[n + 1] = {};
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i == j) continue;
            if((P[i].fi - P[j].fi).length() + P[i].se <= P[j].se + eps) buck[i] = 1;
        }
    }
    vector<pair<point, ld>> res;
    for(int i = 0; i < n; i++){
        if(!buck[i]) res.pb(P[i]);
    }
    return res;
}

int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int n; cin>>n;
    if(!n){
        cout<<0; return 0;
    }
    ld t5 = 0;
    vector<pair<point, ld>> v(n);
    for(int i = 0; i < n; i++){
        int a, b, c; cin>>a>>b>>c;
        v[i] = {point(a, b), c + 10};
        t5 = max(t5, (ld) c + 10);
    }
    if(n == 1){
        cout<<setprecision(20)<<(ld)v[0].se * pi * 2; return 0;
    }
    v = clear(v);
    n = v.size();
    if(n <= 1){
        cout<<setprecision(25)<<t5 * pi * 2; return 0;
    }
    sort(all(v), cmp1);
    vector<pair<point, pair<point, ld>>> res;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            bool f = false;
            auto u = tangents(v[i].fi, v[i].se, v[j].fi, v[j].se, f);
            res.pb({u[0][0], {v[i].fi, v[i].se}});
            res.pb({u[0][1], {v[j].fi, v[j].se}});
            res.pb({u[1][0], {v[i].fi, v[i].se}});
            res.pb({u[1][1], {v[j].fi, v[j].se}});
        }
    }
    if(res.size() <= 1){
        cout<<setprecision(25)<<t5 * pi * 2; return 0;
    }
    res = convexHull(res);
    n = res.size();
    vector<pair<point, pair<point, ld>>> at;
    for(int i = 0; i < n; i++){
        bool f = true;
        for(int j = 0; j < v.size(); j++){
            if((res[i].fi - v[j].fi).length() < v[j].se - eps) f = false;
        }
        if(f) at.pb(res[i]);
    }
    n = at.size();
    ld ans = 0;
    for(int i = 0; i < n; i++){
        int j = i + 1; if(j == n) j = 0;
        if(at[i].se.fi == at[j].se.fi) ans += segmentCircular(at[i].se.fi, at[i].fi, at[j].fi);
        else ans += (at[i].fi - at[j].fi).length();
    }
    cout<<setprecision(25)<<ans;
}
