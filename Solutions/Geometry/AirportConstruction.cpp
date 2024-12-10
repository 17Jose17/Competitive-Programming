/*
  Problem: https://codeforces.com/gym/101471
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

using ld = double;
const ld eps = 1e-9, inf = numeric_limits<ld>::max(), pi = acos(-1);
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

int sgn(ld x){
	if(ge(x, 0)) return 1;
	if(le(x, 0)) return -1;
	return 0;
}
 
bool pointInLine(const point & a, const point & v, const point & p){
	//line a+tv, point p
	return eq((p - a).cross(v), 0);
}
 
bool pointInSegment(const point & a, const point & b, const point & p){
	//segment ab, point p
	return pointInLine(a, b - a, p) && leq((a - p).dot(b - p), 0);
}
 
int intersectLineSegmentInfo(const point & a, const point & v, const point & c, const point & d){
	//line a+tv, segment cd
	point v2 = d - c;
	ld det = v.cross(v2);
	if(eq(det, 0)){
		if(eq((c - a).cross(v), 0)){
			return -1; //infinity points
		}else{
			return 0; //no point
		}
	}else{
		return sgn(v.cross(c - a)) != sgn(v.cross(d - a)); //1: single point, 0: no point
	}
}
 
point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
	//lines a1+tv1, a2+tv2
	//assuming that they intersect
	ld det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
}
 
ld dist(point & a, point & b){
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}
 
bool pointInPerimeter(const vector<point> & P, const point & p){
	int n = P.size();
	for(int i = 0; i < n; i++){
	    int l = i + 1; if(l == n) l = 0;
		if(pointInSegment(P[i], P[l], p)){
			return true;
		}
	}
	return false;
}
 
bool crossesRay(const point & a, const point & b, const point & p){
	return (geq(b.y, p.y) - geq(a.y, p.y)) * sgn((a - p).cross(b - p)) > 0;
}
 
int pointInPolygon(const vector<point> & P, const point & p){
	if(pointInPerimeter(P, p)){
		return 1; //point in the perimeter
	}
	int n = P.size();
	int rays = 0;
	for(int i = 0; i < n; i++){
	    int l = i + 1; if(l == n) l = 0;
		rays += crossesRay(P[i], P[l], p);
	}
	return rays & 1; //0: point outside, 1: point inside
}

bool f(int l, int r, vector<point> P){
    int n = P.size(), i = (l + 1) % (int)P.size();
    vector<point> side1, side2;
    while(i != r){
        side1.pb(P[i]); i = (i + 1) % n;
    }
    while(i != l){
        side2.pb(P[i]); i = (i + 1) % n;
    }
    int m = side1.size(), w = side2.size();
    vector<int> s1, s2;
    for(int j = 0; j < m; j++){
        s1.pb(sgn((P[r] - P[l]).cross(side1[j] - P[l])));
    }
    for(int j = 0; j < w; j++){
        s2.pb(sgn((P[r] - P[l]).cross(side2[j] - P[l])));
    }
    sort(all(s1)); sort(all(s2));
    if((m && abs(s1[0] - s1[m - 1]) > 1) || (w && abs(s2[0] - s2[w - 1]) > 1)) return false;
    int a = 2, b = 2;
    for(auto e : s1) if(e != 0) a = e; for(auto e : s2) if(e != 0) b = e;
    if(a == 2 || b == 2) return true;
    if(a == b) return false;
    return true;
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	int n; cin>>n; vector<point> v(n);
	for(int i = 0; i < n; i++){
	    int a, b; cin>>a>>b; v[i] = point(a, b);
	}
	ld ans = 0;
	for(int i = 0; i < n; i++){
	    for(int j = i + 1; j < n; j++){
	        vector<point> inter;
	        for(int l = 0; l < n; l++){
	            int l1 = l + 1; if(l1 == n) l1 = 0;
	            auto u = intersectLineSegmentInfo(v[i], v[j] - v[i], v[l], v[l1]);
	            if(u == 1) inter.pb(intersectLines(v[i], v[j] - v[i], v[l], v[l1] - v[l]));
	            else if(u == -1) inter.pb(v[l]), inter.pb(v[l1]);
	        }
	        sort(all(inter)); int m = inter.size();
	        int in = 0, fi = 0;
	        for(int l = 0; l < m - 1; l++){
	            if(inter[l] == inter[l + 1]) continue;
	            if(pointInPolygon(v, (inter[l] + inter[l + 1]) / 2) == 0){
	                ans = max(ans, (inter[in] - inter[fi]).length());
	                in = fi = l + 1;
	            }else fi = l + 1;
	        }
	        ans = max(ans, (inter[in] - inter[fi]).length());
	    }
	}
	cout<<setprecision(25)<<ans;
}
