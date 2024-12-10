/*
  Problem: https://codeforces.com/contest/62/problem/C
  Explication:
*/

#include <bits/stdc++.h>
using namespace std;
 
// Holi c:
 
#define ll long long int
//#define ld long double
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()
 
const int Inf = 1e9;
const ll mod = 1e9+7;
const ll INF = 1e15;
 
using ld = long double;
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

struct line{
	point a, v;
	line(): a(), v(){}
	line(const point& a, const point& v): a(a), v(v){}
};
 
istream &operator>>(istream &is, point & p){return is >> p.x >> p.y;}
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

int intersectSegmentsInfo(const point & a, const point & b, const point & c, const point & d){
	//segment ab, segment cd
	point v1 = b - a, v2 = d - c;
	int t = sgn(v1.cross(c - a)), u = sgn(v1.cross(d - a));
	if(t == u){
		if(t == 0){
			if(pointInSegment(a, b, c) || pointInSegment(a, b, d) || pointInSegment(c, d, a) || pointInSegment(c, d, b)){
				return -1; //infinity points
			}else{
				return 0; //no point
			}
		}else{
			return 0; //no point
		}
	}else{
		return sgn(v2.cross(a - c)) != sgn(v2.cross(b - c)); //1: single point, 0: no point
	}
}

point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
	//lines a1+tv1, a2+tv2
	//assuming that they intersect
	ld det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
}

bool lexCompare(const point & a, const point & b){
    if(neq(a.x, b.x))
        return a.x < b.x;
    return a.y < b.y;
}

char segmentType(line seg, point oth){
    if(eq(seg.a.x, seg.v.x))
        return 0;
    if(!lexCompare(seg.a, seg.v))
        swap(seg.a, seg.v);
    return (seg.v - seg.a).cross(oth - seg.a) > 0 ? 1 : -1;
}

bool pointInPerimeter(const vector<point> & P, const point & p){
	int n = P.size();
	for(int i = 0; i < n; i++){
		if(pointInSegment(P[i], P[(i + 1) % n], p)){
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
		return 0; //point in the perimeter
	}
	int n = P.size();
	int rays = 0;
	for(int i = 0; i < n; i++){
		rays += crossesRay(P[i], P[(i + 1) % n], p);
	}
	return rays & 1; //0: point outside, 1: point inside
}

ld perimeterUnionTriangles(vector<vector<point>> & P){
    int n = P.size();
    vector<line> segments(3 * n);
    for(int i = 0; i < n; i++){
        point a = P[i][0], b = P[i][1], c = P[i][2];
        segments[i * 3 + 0] = (line(a, b));
        segments[i * 3 + 1] = (line(b, c));
        segments[i * 3 + 2] = (line(c, a));
    }
    vector<vector<line>> LSegs(3 * n);
    for(int i = 0; i < 3 * n; i++){
        LSegs[i].emplace_back(segments[i]);
        for(int j = 0; j < 3 * n; j++){
            if(i == j) continue;
            int f = 0;
            int m = LSegs[i].size();
            for(int l = 0; l < m; l++){
                int t = intersectSegmentsInfo(LSegs[i][l].a, LSegs[i][l].v, segments[j].a, segments[j].v);
                if(t == 1){
                    auto u = intersectLines(LSegs[i][l].a, LSegs[i][l].v - LSegs[i][l].a, segments[j].a, segments[j].v - segments[j].a);
                    if(u == LSegs[i][l].a || u == LSegs[i][l].v || u == segments[j].a || u == segments[j].v) continue;
                    LSegs[i].emplace_back((line(u, LSegs[i][l].v)));
                    LSegs[i][l] = (line(LSegs[i][l].a, u));
                    break;
                }
            }
        }
    }
    ld ans = 0;
    for(int i = 0; i < 3 * n; i++){
        for(int j = 0; j < LSegs[i].size(); j++){
            int f = 0;
            for(int l = 0; l < n; l++){
                if(pointInPolygon(P[l], (LSegs[i][j].a + LSegs[i][j].v) / 2))
                    f = 1;
            }
            if(!f) ans += (LSegs[i][j].a - LSegs[i][j].v).length();
        }
    }
    return ans;
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	int n; cin>>n;
	vector<vector<point>> vTs(n);
	ld res = 0;
	for(int i = 0; i < n; i++){
	    int a, b, c, d, e, f, g, h; cin>>a>>b>>c>>d>>e>>f;
	    point t1(a, b), t2(c, d), t3(e, f);
	    vTs[i].pb(t1); vTs[i].pb(t2); vTs[i].pb(t3);
	}
	ld ans = perimeterUnionTriangles(vTs);
	cout<<setprecision(20)<<ans;
}
