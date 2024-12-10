/*
  Problem: https://codeforces.com/gym/104196/problem/D
  Explication:
*/

#include <bits/stdc++.h>
using namespace std;

// Holi c:

#define ll long long int
#define dl double
#define ld long double
#define ff __float128
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()

const int Inf = 1e9;
const ll mod = 1e9+7;
const ll INF = 1e18;

const ld eps = 1e-15, inf = numeric_limits<ld>::max(), pi = acos(-1);
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
    ld ang() const{ ld a = atan2l(y, x); a += le(a, 0) ? 2*pi : 0; return a;
    }
    ld ang(const point &p) const{ return abs(ang() - p.ang());}
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

ld area(vector<point> & P){
	int n = P.size();
	ld ans = 0;
	for(int i = 0; i < n; i++){
		ans += P[i].cross(P[(i + 1) % n]);
	}
	return (ans / 2);
}

point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
	//lines a1+tv1, a2+tv2
	//assuming that they intersect
	ld det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
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
				return 0; //infinity points
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

ld distEuc(point p1, point p2){
    return sqrtl((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

bool check(point p1, point p2, point ctt, point cent){
    if(intersectSegmentsInfo(p1, p2, ctt, cent) == 1) return false;
    return true;
}

ld areaCirc(point p1, point p2, point ctt, point cent){
    ld aux = distEuc(ctt, p1);
    point a = p1 - ctt, b = p2 - ctt;
    ld prr = acosl(a.dot(b) / (a.length() * b.length()));
    if(abs(abs((a.dot(b)) / (a.length() * b.length())) - 1) <= eps) return aux * aux * pi / 2;
    if(check(p1, p2, ctt, cent)) return aux * aux * prr / 2 - abs(a.cross(b)) / 2;
    return aux * aux * pi - (aux * aux * prr / 2) + abs(a.cross(b)) / 2;
}

point midPoint(point p, point q){
    point nuv;
    nuv.x = p.x + (q.x - p.x) / 2;
    nuv.y = p.y + (q.y - p.y) / 2;
    return nuv;
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
    point cent; ld rad; cin>>cent>>rad;
    int n; cin>>n;
    vector<point> vp(n); vector<point> polyg; vector<point> res;
    for(int i = 0; i < n; i++) cin>>vp[i];
    for(int i = 0 ; i < n; i++){
        point rep = vp[i] - cent;
        rep = rep.unit() * (rad * rad / rep.length());
        rep += cent;
        polyg.pb(rep);
    }

    ld res_p = area(polyg);
    for(int i = 0; i < n; i++){
        
        point p = cent, q = polyg[i], r = polyg[(i + 1) % n];
        point midpq = midPoint(p, q), midpr = midPoint(p, r);
	      point lpq, lpr;
	      lpq.x = q.x - (q.y - p.y); lpq.y = q.y + (q.x - p.x); 
	      lpr.x = r.x - (r.y - p.y); lpr.y = r.y + (r.x - p.x);
	      lpq.x -= (q.x - midpq.x); lpq.y -= (q.y - midpq.y);
	      lpr.x -= (r.x - midpr.x); lpr.y -= (r.y - midpr.y);
        point ctt = intersectLines(midpq, lpq - midpq, midpr, lpr - midpr);
        
        if(((vp[i] - p).cross(vp[(i + 1) % n] - p)) == 0) continue;
        
        if((q - cent).cross(r - cent) > 0) res_p += areaCirc(q, r, ctt, p);
        else if((q - cent).cross(r - cent) < 0) res_p -= areaCirc(q, r, ctt, p);
    }
    cout<<setprecision(20)<<abs(res_p);
}
