/*
  Problem: https://qoj.ac/contest/442/problem/1196
  Explication:
*/

#include <bits/stdc++.h>
using namespace std;

// Holi c:

#define ll long long int
//#define ld long double
#define ii __int128
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()

const int Inf = 1e9;
const ll mod =  1e9 + 7;
const ll INF = 4e18;
const int maxn = 2e5 + 5;

using ld = long double;
const ld eps = 1e-12, inf = numeric_limits<ld>::max(), pi = acos(-1);
// For use with integers, just set eps=0 and everything remains the same
bool geq(ld a, ld b){return a-b >= -eps;}     //a >= b
bool leq(ld a, ld b){return b-a >= -eps;}     //a <= b
bool ge(ld a, ld b){return a-b > eps;}        //a > b
bool le(ld a, ld b){return b-a > eps;}        //a < b
bool eq(ld a, ld b){return abs(a-b) <= eps;}  //a == b
bool neq(ld a, ld b){return abs(a-b) > eps;}  //a != b

struct point{
	ld x, y;
	point(): x(0), y(0){}
	point(ld x, ld y): x(x), y(y){}

	point operator+(const point & p) const{return point(x + p.x, y + p.y);}
	point operator-(const point & p) const{return point(x - p.x, y - p.y);}
	point operator*(const ld & k) const{return point(x * k, y * k);}
	point operator/(const ld & k) const{return point(x / k, y / k);}

	point operator+=(const point & p){*this = *this + p; return *this;}
	point operator-=(const point & p){*this = *this - p; return *this;}
	point operator*=(const ld & p){*this = *this * p; return *this;}
	point operator/=(const ld & p){*this = *this / p; return *this;}

	point rotate(const ld & a) const{return point(x*cos(a) - y*sin(a), x*sin(a) + y*cos(a));}
	point perp() const{return point(-y, x);}
	ld ang() const{
		ld a = atan2l(y, x); a += le(a, 0) ? 2*pi : 0; return a;
	}
	ld dot(const point & p) const{return x * p.x + y * p.y;}
	ld cross(const point & p) const{return x * p.y - y * p.x;}
	ld norm() const{return x * x + y * y;}
	ld length() const{return sqrtl(x * x + y * y);}
	point unit() const{return (*this) / length();}

	bool operator==(const point & p) const{return eq(x, p.x) && eq(y, p.y);}
	bool operator!=(const point & p) const{return !(*this == p);}
	bool operator<(const point & p) const{return le(x, p.x) || (eq(x, p.x) && le(y, p.y));}
	bool operator>(const point & p) const{return ge(x, p.x) || (eq(x, p.x) && ge(y, p.y));}
	bool half(const point & p) const{return le(p.cross(*this), 0) || (eq(p.cross(*this), 0) && le(p.dot(*this), 0));}
};

ostream &operator<<(ostream &os, const point & p){return os << "(" << p.x << ", " << p.y << ")";}

ld f1(ld a, ld b, ld c, ld x){
    return a * x * x + b * x + c;
}

ld f(ld a, ld b, ld c, ld l, ld r){
    for(int i = 0; i < 75; i++){
        ld m1 = l + (r - l) / 3, m2 = r - (r - l) / 3;
        ld w1 = f1(a, b, c, m1), w2 = f1(a, b, c, m2);
        if(w1 > w2) l = m1;
        else r = m2;
    }
    return f1(a, b, c, l);
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	int n; cin>>n; vector<pair<point, point>> v(n);
	for(int i = 0; i < n; i++){ int a, b, c, d; cin>>a>>b>>c>>d; v[i] = {point(a, b), point(c, d)}; }
	if(n == 1){ cout<<0; return 0; }
	vector<pair<ld, ld>> ix, iy;
    //Velocidad, posicion
    for(auto e : v) ix.pb({e.se.x, e.fi.x}); for(auto e : v) iy.pb({e.se.y, e.fi.y}); sort(all(ix)); sort(all(iy));
	//velocidad, inicio, tiempo
    vector<pair<ld, pair<ld, ld>>> minx, maxx, miny, maxy;
    //Computar todos los intervalos de tiempo ¿Mismas velocidades? -> mayor maxx, mayor maxy, menor minx, menor miny
    //maxx
    for(int i = 0; i < n; i++){
        if(!maxx.size()){ maxx.pb({ix[i].fi, {ix[i].se, 0}}); continue; }
        while(maxx.size() && (((maxx[maxx.size() - 1].se.fi - (maxx[maxx.size() - 1].fi * maxx[maxx.size() - 1].se.se) - ix[i].se) / (ix[i].fi - maxx[maxx.size() - 1].fi) <= maxx[maxx.size() - 1].se.se)
         || (maxx[maxx.size() - 1].fi == ix[i].fi && maxx[maxx.size() - 1].se.fi - maxx[maxx.size() - 1].fi * maxx[maxx.size() - 1].se.se <= ix[i].se + eps))){
            maxx.pop_back();
        }
        ld at = 0; if(maxx.size()) at = (maxx[maxx.size() - 1].se.fi - (maxx[maxx.size() - 1].fi * maxx[maxx.size() - 1].se.se) - ix[i].se) / (ix[i].fi - maxx[maxx.size() - 1].fi);
        maxx.pb({ix[i].fi, {ix[i].se + ix[i].fi * at, at}});
    } reverse(all(ix));
    //minx
    for(int i = 0; i < n; i++){
        if(!minx.size()){ minx.pb({ix[i].fi, {ix[i].se, 0}}); continue; }
        while(minx.size() && (((minx[minx.size() - 1].se.fi - (minx[minx.size() - 1].fi * minx[minx.size() - 1].se.se) - ix[i].se) / (ix[i].fi - minx[minx.size() - 1].fi) <= minx[minx.size() - 1].se.se)
         || (minx[minx.size() - 1].fi == ix[i].fi && minx[minx.size() - 1].se.fi - minx[minx.size() - 1].fi * minx[minx.size() - 1].se.se >= ix[i].se - eps))){
            minx.pop_back();
        }
        ld at = 0; if(minx.size()) at = (minx[minx.size() - 1].se.fi - (minx[minx.size() - 1].fi * minx[minx.size() - 1].se.se) - ix[i].se) / (ix[i].fi - minx[minx.size() - 1].fi);
        minx.pb({ix[i].fi, {ix[i].se + ix[i].fi * at, at}});        
    }
    //maxy
    for(int i = 0; i < n; i++){
        if(!maxy.size()){ maxy.pb({iy[i].fi, {iy[i].se, 0}}); continue; }
        while(maxy.size() && (((maxy[maxy.size() - 1].se.fi - (maxy[maxy.size() - 1].fi * maxy[maxy.size() - 1].se.se) - iy[i].se) / (iy[i].fi - maxy[maxy.size() - 1].fi) <= maxy[maxy.size() - 1].se.se)
         || (maxy[maxy.size() - 1].fi == iy[i].fi && maxy[maxy.size() - 1].se.fi - maxy[maxy.size() - 1].fi * maxy[maxy.size() - 1].se.se <= iy[i].se + eps))){
            maxy.pop_back();
        }
        ld at = 0; if(maxy.size()) at = (maxy[maxy.size() - 1].se.fi - (maxy[maxy.size() - 1].fi * maxy[maxy.size() - 1].se.se) - iy[i].se) / (iy[i].fi - maxy[maxy.size() - 1].fi);
        maxy.pb({iy[i].fi, {iy[i].se + iy[i].fi * at, at}});        
    } reverse(all(iy));
    //miny
    for(int i = 0; i < n; i++){
        if(!miny.size()){ miny.pb({iy[i].fi, {iy[i].se, 0}}); continue; }
        while(miny.size() && (((miny[miny.size() - 1].se.fi - (miny[miny.size() - 1].fi * miny[miny.size() - 1].se.se) - iy[i].se) / (iy[i].fi - miny[miny.size() - 1].fi) <= miny[miny.size() - 1].se.se)
         || (miny[miny.size() - 1].fi == iy[i].fi && miny[miny.size() - 1].se.fi - miny[miny.size() - 1].fi * miny[miny.size() - 1].se.se >= iy[i].se - eps))){
            miny.pop_back();
        }
        ld at = 0; if(miny.size()) at = (miny[miny.size() - 1].se.fi - (miny[miny.size() - 1].fi * miny[miny.size() - 1].se.se) - iy[i].se) / (iy[i].fi - miny[miny.size() - 1].fi);
        miny.pb({iy[i].fi, {iy[i].se + iy[i].fi * at, at}});          
    }
    minx.pb({minx[minx.size() - 1].fi, {minx[minx.size() - 1].se.fi + minx[minx.size() - 1].fi * Inf, Inf}}); 
    maxx.pb({maxx[maxx.size() - 1].fi, {maxx[maxx.size() - 1].se.fi + maxx[maxx.size() - 1].fi * Inf, Inf}}); 
    miny.pb({miny[miny.size() - 1].fi, {miny[miny.size() - 1].se.fi + miny[miny.size() - 1].fi * Inf, Inf}}); 
    maxy.pb({maxy[maxy.size() - 1].fi, {maxy[maxy.size() - 1].se.fi + maxy[maxy.size() - 1].fi * Inf, Inf}});

    vector<pair<ld, pair<ld, pair<ld, ld>>>> vx, vy;
    int i = 0, j = 0;
    //vx
    while(i < minx.size() - 1 && j < maxx.size() - 1){
        ld ti = minx[i].se.se, tf = minx[i + 1].se.se, ti1 = maxx[j].se.se, tf1 = maxx[j + 1].se.se;
        ld tr = max(ti, ti1);
        ld v = minx[i].fi, v1 = maxx[j].fi;
        ld d = minx[i].se.fi + (tr - ti) * v, d1 = maxx[j].se.fi + (tr - ti1) * v1; 
        ld re = v1 + v * -1; ld er = d1 + d * -1;
        vx.pb({re, {er, {max(ti, ti1), min(tf, tf1)}}});
        if(tf < tf1) i++;
        else j++;
    }
    i = j = 0;
    //vy
    while(i < miny.size() - 1 && j < maxy.size() - 1){
        ld ti = miny[i].se.se, tf = miny[i + 1].se.se, ti1 = maxy[j].se.se, tf1 = maxy[j + 1].se.se;
        ld tr = max(ti, ti1);
        ld v = miny[i].fi, v1 = maxy[j].fi;
        ld d = miny[i].se.fi + (tr - ti) * v, d1 = maxy[j].se.fi + (tr - ti1) * v1; 
        ld re = v1 + v * -1; ld er = d1 + d * -1;
        vy.pb({re, {er, {max(ti, ti1), min(tf, tf1)}}});
        if(tf < tf1) i++;
        else j++;
    }
 
    ld ans = INF; i = j = 0;
    while(i < vx.size() && j < vy.size()){
        ld a = vx[i].fi, b = vx[i].se.fi, ti = vx[i].se.se.fi, tf = vx[i].se.se.se;
        ld a1 = vy[j].fi, b1 = vy[j].se.fi, ti1 = vy[j].se.se.fi, tf1 = vy[j].se.se.se;
        b += a * (max(ti, ti1) - ti); b1 += a1 * (max(ti, ti1) - ti1);
        ans = min(ans, f(a * a1, b * a1 + b1 * a, b * b1, 0, min(tf, tf1) - max(ti, ti1)));
        if(tf < tf1) i++;
        else j++;
    }
    cout<<setprecision(25)<<ans;
}
