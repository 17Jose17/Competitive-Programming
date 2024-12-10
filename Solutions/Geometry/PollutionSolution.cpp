/*
  Problem: https://codeforces.com/gym/101208
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
const ll mod = 1e9 + 7;
const ll INF = 1e18;

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

vector<vector<point>> triangulate(vector<point> & P){
    int n = P.size();
    vector<int> next(n);
    for(int i = 0; i < n - 1; i++) next[i] = i + 1;
    auto is_ear = [&](int i, int j, int k){
        if(sgn((P[j] - P[i]).cross(P[k] - P[i])) <= 0) return false;
        for(int l = next[k]; l != i; l = next[l])
            if(sgn((P[i] - P[l]).cross(P[j] - P[l])) >= 0 &&
             sgn((P[j] - P[l]).cross(P[k] - P[l])) >= 0 &&
             sgn((P[k] - P[l]).cross(P[i] - P[l])) >= 0) return false;
        return true;
    };
    vector<vector<point>> res;
    for(int i = 0; next[next[i]] != i;){
        if(is_ear(i, next[i], next[next[i]])){
            res.pb({P[i], P[next[i]], P[next[next[i]]]});
            next[i] = next[next[i]];
        }else i = next[i];
    }
    return res;
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

vector<point> cutPolygon(const vector<point> & P, const point & a, const point & v){
	//returns the part of the convex polygon P on the left side of line a+tv
	int n = P.size();
	vector<point> lhs;
	for(int i = 0; i < n; ++i){
		if(geq(v.cross(P[i] - a), 0)){
			lhs.push_back(P[i]);
		}
		if(intersectLineSegmentInfo(a, v, P[i], P[(i+1)%n]) == 1){
			point p = intersectLines(a, v, P[i], P[(i+1)%n] - P[i]);
			if(p != P[i] && p != P[(i+1)%n]){
				lhs.push_back(p);
			}
		}
	}
	return lhs;
}

ld area(vector<point> & P){
    int n = P.size();
    ld ans = 0;
    for(int i = 0; i < n; i++){
        ans += P[i].cross(P[(i + 1) % n]);
    }
    return ans / 2;
}

bool pointInLine(const point & a, const point & v, const point & p){
	//line a+tv, point p
	return eq((p - a).cross(v), 0);
}

bool pointInSegment(const point & a, const point & b, const point & p){
	//segment ab, point p
	return pointInLine(a, b - a, p) && leq((a - p).dot(b - p), 0);
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

vector<point> intersectSegmentCircle(const point & a, const point & b, const point & c, ld r){
	//segment ab, circle with center c and radius r
	vector<point> P = intersectLineCircle(a, b - a, c, r), ans;
	for(const point & p : P){
		if(pointInSegment(a, b, p)) ans.push_back(p);
	}
	return ans;
}

ld calc(vector<point> P, ld r){
    int n = P.size();
    vector<point> L;
    for(int i = 0; i < n; i++){
        auto u = intersectSegmentCircle(P[i], P[(i + 1) % n], point(0, 0), r);
        for(auto e : u) L.pb(e);
    }
    if(L.size() == 1) return 0;
    polarSort(L, point(0, 0), point(1, 0));
    ld ans = 0;
    for(int i = 0; i < L.size(); i += 2){
        if(le(L[i].cross(L[i + 1]), 0)) swap(L[i], L[i + 1]);
        P = cutPolygon(P, L[i], L[i + 1] - L[i]);
        ld tet = acosl(L[i].dot(L[i + 1]) / L[i].length() / L[i + 1].length());
        ans += r * r * tet / 2 - abs(L[i].cross(L[i + 1])) / 2;
    }
    ans += abs(area(P));
    return ans;
}

int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int n; ld r; cin>>n>>r;
    vector<point> P1(n);
    for(int i = 0; i < n; i++){
        int a, b; cin>>a>>b;
        P1[i] = point(a, b);
    }
    vector<point> P;
    for(int i = 0; i < n; i++){
        if(neq((P1[(i + 1) % n] - P1[i]).cross(P1[(i - 1 + n) % n] - P1[i]), 0)) P.pb(P1[i]);
    }
    auto T = triangulate(P);
    int k = T.size();
    for(int i = 0; i < k; i++){
        if(le(area(T[i]), 0)) reverse(all(T[i]));
    }
    point p(0, 0);
    k = T.size();
    ld ans = 0;
    for(int i = 0; i < k; i++){
        bool fla = false, flb = false;
        int t = T[i].size();
        for(int j = 0; j < T[i].size(); j++){
            if(le(T[i][j].length(), r)) fla = true;
            auto u = intersectSegmentCircle(T[i][j], T[i][(j + 1) % t], p, r);
            if(u.size()) flb = true;
        }
        if(flb) ans += calc(T[i], r);
        else if(fla) ans += abs(area(T[i]));
    }
    cout<<setprecision(25)<<ans;
}
