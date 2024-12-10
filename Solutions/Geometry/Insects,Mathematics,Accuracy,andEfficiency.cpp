/*
  Problem: https://codeforces.com/gym/105053/problem/I
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
const ll INF = 1e18;
const int maxn = 2e5 + 5;

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
		ld a = atan2l(y, x); a += le(a, 0) ? pi * 2 : 0;
        return a;
	}
    ld ang(const point & p) const{ return abs(ang() - p.ang()); }
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

void polarSort(vector<point> & P, const point & o, const point & v){
	//sort points in P around o, taking the direction of v as first angle
	sort(P.begin(), P.end(), [&](const point & a, const point & b){
		return point((a - o).half(v), 0) < point((b - o).half(v), (a - o).cross(b - o));
	});
}

vector<point> convexHull(vector<point> P){
	sort(P.begin(), P.end());
	vector<point> L, U;
	for(int i = 0; i < P.size(); i++){
		while(L.size() >= 2 && leq((L[L.size() - 2] - P[i]).cross(L[L.size() - 1] - P[i]), 0)){
			L.pop_back();
		}
		L.push_back(P[i]);
	}
	for(int i = P.size() - 1; i >= 0; i--){
		while(U.size() >= 2 && leq((U[U.size() - 2] - P[i]).cross(U[U.size() - 1] - P[i]), 0)){
			U.pop_back();
		}
		U.push_back(P[i]);
	}
	L.pop_back();
	U.pop_back();
	L.insert(L.end(), U.begin(), U.end());
	return L;
}

ld area(vector<point> & P){
	int n = P.size();
	ld ans = 0;
	for(int i = 0; i < n; i++){
		ans += P[i].cross(P[(i + 1) % n]);
	}
	return abs(ans / 2);
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

ld f1(point p0, point p1, point p2){
    return abs((p1 - p0).cross(p2 - p0)) / 2;
}

ld f(point p0, point p1, point t1, point t2){
    ld ans = 0, l = 0, r = 0;
    
    r = acosl(p1.dot(p0) / p1.length() / p0.length());
    
    point rt = p1.rotate(r);
    if(!(abs(rt.x - p0.x) <= .5 && abs(rt.y - p0.y) <= .5)) r = 2 * pi - r;
    
    for(int i = 0; i < 100; i++){
        ld m1 = l + (r - l) / 3, m2 = r - (r - l) / 3;
        ld w1 = f1(t1, t2, p1.rotate(m1)), w2 = f1(t1, t2, p1.rotate(m2));
        if(w1 < w2) l = m1;
        else r = m2;
    }
    return f1(t1, t2, p1.rotate(l));
}

ld calc(point p0, point p1, int t0, int t1, vector<point> & v){
    ld ans = f(p0, p1, v[t0], v[t1]);
    return ans;
}

ld solveto2(point p0, point p1, point p2, point p3){
    return max(f(p0, p1, p2, p3), f(p1, p0, p3, p2));
}


int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	int n; ld r; cin>>n>>r; r += eps;
	vector<point> v(n);
	for(int i = 0; i < n; i++){ int a, b; cin>>a>>b; v[i] = point(a, b); }
	if(n == 1){ cout<<0<<'\n'; return 0; }
	sort(all(v));
	v = convexHull(v); n = v.size(); ld res = area(v); reverse(all(v));
	vector<point> side;
	for(int i = 0; i < v.size(); i++){
	    auto u = intersectLineCircle(v[i], v[(i + 1) % (int)v.size()] - v[i], point(0, 0), r);
	    if(n == 2 && i + 1 == v.size()) continue;
	    side.pb(u[0]); side.pb(u[1]);
	}
	polarSort(side, point(0, 0), point(-1, 0));
	reverse(all(side)); vector<point> u; u.pb(v[v.size() - 1]); for(int i = 0; i < v.size() - 1; i++) u.pb(v[i]); v = u;
	int j = 0, l = n / 2, m = side.size(); n = v.size(); ld ans = 0, zes = 0;
	if(n == 2){ cout<<setprecision(20)<<solveto2(side[0], side[1], v[0], v[1])<<'\n'; return 0; }
	for(int i = 0; i < m; i++){
	    while((v[l] - side[i]).cross(v[(l + 1) % n] - side[i]) < eps){
	        if(i) zes -= f1(v[j], v[l], v[(l + 1) % n]);
	        l++; if(l == n) l = 0;
	    }
	    if(!i) j = (l + 1) % n;
	    while((v[j] - side[(i + 1) % m]).cross(v[(j + 1) % n] - side[(i + 1) % m]) > -eps){
	        zes += f1(v[l], v[j], v[(j + 1) % n]);
	        j++; if(j == n) j = 0;
	    }
	    if(abs((v[j] - side[(i + 1) % m]).cross(v[(j - 1 + n) % n] - side[(i + 1) % m])) <= eps){
	        j--;
	        if(j < 0) j += n;
	        zes -= f1(v[l], v[j], v[(j + 1) % n]);
	    }
	
	    ans = max(ans, calc(side[i], side[(i + 1) % m], l, j, v) - (zes));
	}
	cout<<setprecision(20)<<res + ans;
}
