/*
  Problem: https://acm.timus.ru/problem.aspx?space=1&num=1662
  Explication:
*/

#include <bits/stdc++.h>
using namespace std;
 
// Holi c:
 
#define ll long long int
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()
 
const int Inf = 1e9;
const ll mod = 1e9+7;
const ll INF = 1e18;
 
using ld = long double;
const ld eps = 1e-9, inf = numeric_limits<ld>::max(), pi = acos(-1);
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

istream &operator>>(istream &is, point & p){return is >> p.x >> p.y;}
ostream &operator<<(ostream &os, const point & p){return os << "(" << p.x << ", " << p.y << ")";}

int sgn(ld x){
	if(ge(x, 0)) return 1;
	if(le(x, 0)) return -1;
	return 0;
}

void polarSort(vector<point> & P, const point & o, const point & v){
	//sort points in P around o, taking the direction of v as first angle
	sort(P.begin(), P.end(), [&](const point & a, const point & b){
		return point((a - o).half(v), 0) < point((b - o).half(v), (a - o).cross(b - o));
	});
}

vector<point> intersectLineCircle(const point & a, const point & v, const point & c, ld r){
	//line a+tv, circle with center c and radius r
	ld h2 = r*r - v.cross(c - a) * v.cross(c - a) / v.norm();
	point p = a + v * v.dot(c - a) / v.norm();
	if(eq(h2, 0)) return {p}; //line tangent to circle
	else if(le(h2, 0)) return {}; //no intersection
	else{
		point u = v.unit() * sqrt(h2);
		return {p - u, p + u}; //two points of intersection (chord)
	}
}

vector<point> intersectionCircles(const point & c1, ld r1, const point & c2, ld r2){
	point d = c2 - c1;
	ld d2 = d.norm();
	if(eq(d2, 0)) return {};
	ld pd = (d2 + r1*r1 - r2*r2) / 2;
	ld h2 = r1*r1 - pd*pd/d2;
	point p = c1 + d*pd/d2;
	if(eq(h2, 0)) return {p};
	else if(le(h2, 0)) return {};
	else{
		point u = d.perp() * sqrt(h2/d2);
		return {p - u, p + u};
	}
}

bool pointInLine(const point & a, const point & v, const point & p){
	return eq((p - a).cross(v), 0);
}

bool pointInSegment(const point & a, const point & b, const point & p){
	return pointInLine(a, b - a, p) && leq((a - p).dot(b - p), 0);
}

ld distancePointLine(const point & a, const point & v, const point & p){
	return abs(v.cross(p - a)) / v.length();
}

point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
	ld det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
}

ld distancePointSegment(const point & a, const point & b, const point & p){
    if(eq((b - a).cross(p - a), 0)) return min((a - p).length(), (b - p).length());
    auto it = intersectLines(a, b - a, p, (b - a).perp());
    if(pointInSegment(a, b, it)) return distancePointLine(a, b - a, p);
    return min((a - p).length(), (b - p).length());
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
		return 1; //point in the perimeter
	}
	int n = P.size();
	int rays = 0;
	for(int i = 0; i < n; i++){
		rays += crossesRay(P[i], P[(i + 1) % n], p);
	}
	return rays & 1; //0: point outside, 1: point inside
}

ld distancePointPolygon(vector<point> P, point p){
    int n = P.size();
    ld ans = INF;
    for(int i = 0; i < n; i++){
        ans = min(ans, distancePointSegment(P[i], P[(i + 1) % n], p));
    }
    return ans;
}

pair<bool, point> solve1(vector<point> P, point c, point a, point b, ld r1, ld r2){
    int n = P.size();
    point ans;
    for(int i = 0; i < n; i++) if(ge((a - P[i]).length(), r2) || ge((b - P[i]).length(), r2)) return {false, ans};
    ld l = 0, r = acosl((a - c).dot(b - c) / (a - c).length() / (b - c).length());
    point d = a - c;
    for(int i = 0; i < 100; i++){
        ld m1 = l + (r - l) / 3, m2 = r - (r - l) / 3;
        point p1 = c + d.rotate(m1), p2 = c + d.rotate(m2);
        ld f1 = distancePointPolygon(P, p1), f2 = distancePointPolygon(P, p2);
        if(f1 < f2) l = m1;
        else r = m2;
    }
    point p0 = c + d.rotate(l);
    if(geq(distancePointPolygon(P, p0), r1)) return {true, p0};
    return {false, ans};
}

pair<bool, point> solve(vector<point> P){
    int n = P.size(); point ans;
    ld r1, r2; cin>>r1>>r2;
    vector<vector<point>> C(n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i == j) continue;
            auto u = intersectionCircles(P[i], r2, P[j], r2);
            if(u.size() == 0) return {false, ans};
            else if(u.size() == 1){
                point at = u[0];
                for(int l = 0; l < n; l++) if(ge((P[l] - at).length(), r1)) return {false, ans};
                for(int l = 0; l < n; l++) if(le(distancePointSegment(P[l], P[(l + 1) % n], at), r1)) return {false, ans};
                return {true, at};
            }else{
                C[i].pb(u[0]); C[i].pb(u[1]);
                C[j].pb(u[0]); C[j].pb(u[0]);
            }
        }
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            auto u = intersectLineCircle(P[j], P[(j + 1) % n] - P[j], P[i], r2);
            for(auto e : u) C[i].pb(e);
        }
    }
    for(int i = 0; i < n; i++) polarSort(C[i], P[i], point(1, 0));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < C[i].size(); j++){
            point p0 = C[i][j], p1 = C[i][(j + 1) % (int) C[i].size()];
            ld te = acosl((p0 - P[i]).dot(p1 - P[i]) / (p0 - P[i]).length() / (p1 - P[i]).length());
            point p2 = p0 - P[i]; p2.rotate(te); p2 += P[i];
            if(pointInPolygon(P, p2) == 1) continue;
            auto rep = solve1(P, P[i], p0, p1, r1, r2);
            if(rep.fi == true) return rep;
        }
    }
    return {false, ans};
}

int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int n; cin>>n;
    vector<point> P(n);
    for(int i = 0; i < n; i++){
        int a, b; cin>>a>>b;
        P[i] = point(a, b);
    }
    auto u = solve(P);
    if(u.fi == false){
        cout<<"No solution";
    }else{
        cout<<setprecision(25)<<u.se.x<<" "<<u.se.y;
    }
}
