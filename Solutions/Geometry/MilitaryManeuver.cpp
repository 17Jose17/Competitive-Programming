/*
  Problem: https://codeforces.com/gym/104869/problem/G
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
const ld eps = 1e-9, inf = 1e8, pi = acos(-1);
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

struct plane{
	point a, v;
	plane(): a(), v(){}
	plane(const point& a, const point& v): a(a), v(v){}

	point intersect(const plane& p) const{
		ld t = (p.a - a).cross(p.v) / v.cross(p.v);
		return a + v*t;
	}

	bool outside(const point& p) const{ // test if point p is strictly outside
		return le(v.cross(p - a), 0);
	}

	bool inside(const point& p) const{ // test if point p is inside or in the boundary
		return geq(v.cross(p - a), 0);
	}

	bool operator<(const plane& p) const{ // sort by angle
		auto lhs = make_tuple(v.half({1, 0}), ld(0), v.cross(p.a - a));
		auto rhs = make_tuple(p.v.half({1, 0}), v.cross(p.v), ld(0));
		return lhs < rhs;
	}

	bool operator==(const plane& p) const{ // paralell and same directions, not really equal
		return eq(v.cross(p.v), 0) && ge(v.dot(p.v), 0);
	}
};

vector<point> halfPlaneIntersection(vector<plane> planes){
	planes.push_back({{0, -inf}, {1, 0}});
	planes.push_back({{inf, 0}, {0, 1}});
	planes.push_back({{0, inf}, {-1, 0}});
	planes.push_back({{-inf, 0}, {0, -1}});
	sort(planes.begin(), planes.end());
	planes.erase(unique(planes.begin(), planes.end()), planes.end());
	deque<plane> ch;
	deque<point> poly;
	for(const plane& p : planes){
		while(ch.size() >= 2 && p.outside(poly.back())) ch.pop_back(), poly.pop_back();
		while(ch.size() >= 2 && p.outside(poly.front())) ch.pop_front(), poly.pop_front();
		if(p.v.half({1, 0}) && poly.empty()) return {};
		ch.push_back(p);
		if(ch.size() >= 2) poly.push_back(ch[ch.size()-2].intersect(ch[ch.size()-1]));
	}
	while(ch.size() >= 3 && ch.front().outside(poly.back())) ch.pop_back(), poly.pop_back();
	while(ch.size() >= 3 && ch.back().outside(poly.front())) ch.pop_front(), poly.pop_front();
	poly.push_back(ch.back().intersect(ch.front()));
	return vector<point>(poly.begin(), poly.end());
}

ld f(point a, point b, point c){
    ld w = abs((b - a).cross(c - a)) / 2, at = 1;
    if((b - c).length() <= eps || w <= eps) return 0;
    ld h = w * 2 / (b - c).length(), mx = max((b - a).length(), (c  - a).length());
    ld v = sqrt(mx * mx - h * h); ld u = v - (b - c).length();
    if((b - a).cross(c - a) <= -eps) at = -1;
    return h * (v - u) * (h * h * 3.0 + u * u + u * v + v * v) / 12.0 * at;
}

pair<point, point> bisector(point p0, point p1){
    point pm = (p0 + p1) / 2.0;
    point b = ((p0 - p1).perp());
    return {pm, b};
}

int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int ax, ay, bx, by; cin>>ax>>ay>>bx>>by;
    point i0(ax, ay), i1(bx, ay), i2(bx, by), i3(ax, by);
    int n; cin>>n;
    vector<point> v(n);
    for(int i = 0; i < n; i++){
        int a, b; cin>>a>>b; v[i] = point(a, b);
    }
    ld ans = 0;
    //Closest
    for(int i = 0; i < n; i++){
        vector<plane> vp;
        for(int j = 0; j < n; j++){
            if(i == j) continue;
            auto u = bisector(v[i], v[j]);
            vp.pb(plane(u.fi, u.se));
        }
        vp.pb(plane(i0, i1 - i0)); vp.pb(plane(i1, i2 - i1)); vp.pb(plane(i2, i3 - i2)); vp.pb(plane(i3, i0 - i3));
        auto u = halfPlaneIntersection(vp);
        for(int j = 0; j < u.size(); j++){
            int l = j + 1; if(l == u.size()) l = 0;
            ans -= f(v[i], u[j], u[l]);
        }
    }
    //Farthest
       for(int i = 0; i < n; i++){
        vector<plane> vp;
        for(int j = 0; j < n; j++){
            if(i == j) continue;
            auto u = bisector(v[i], v[j]);
            vp.pb(plane(u.fi, u.se * -1));
        }
        vp.pb(plane(i0, i1 - i0)); vp.pb(plane(i1, i2 - i1)); vp.pb(plane(i2, i3 - i2)); vp.pb(plane(i3, i0 - i3));
        auto u = halfPlaneIntersection(vp);
        for(int j = 0; j < u.size(); j++){
            int l = j + 1; if(l == u.size()) l = 0;
            ans += f(v[i], u[j], u[l]);
        }
    }
    cout<<setprecision(35)<<abs(ans) * pi / (ld)(bx - ax) / (ld)(by - ay);
}
