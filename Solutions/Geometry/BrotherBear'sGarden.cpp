/*
  Problem: https://acm.timus.ru/problem.aspx?space=1&num=1681
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

const int Inf = 1e4 + 1;
const ll mod = 1e9+7;
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

ld area(vector<point> & P){
	int n = P.size();
	ld ans = 0;
	for(int i = 0; i < n; i++){
		ans += P[i].cross(P[(i + 1) % n]);
	}
	return abs(ans / 2);
}

point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
	//lines a1+tv1, a2+tv2
	//assuming that they intersect
	ld det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
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

vector<point> f(vector<point> P, vector<point> v, vector<vector<int>> L){
    vector<point> t; t.pb(P[0]); t.pb(P[2]);
	int n = L.size();
	for(int i = 0; i < n; i++){
		point a = v[i], b = v[L[i][0]], c = v[L[i][1]];
		bool fl = false;
		for(int j = 0; j < t.size(); j++) if(sgn((b - a).cross(t[j] - a)) == -1) fl = true;
		if(!fl) P = cutPolygon(P, a, b - a);
		else P = cutPolygon(P, a, a - b);
		fl = false;
		for(int j = 0; j < t.size(); j++) if(sgn((c - a).cross(t[j] - a)) == -1) fl = true;
		if(!fl) P = cutPolygon(P, a, c - a);
		else P = cutPolygon(P, a, a - c);
	}
	return P;
}

int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int n, k; cin>>n>>k;
    if(k == 1 || k == n - 1){
        cout<<1; return 0;
    }
    if(k * 2 == n){
        cout<<0; return 0;
    }
    vector<point> v(n);
    for(int i = 0; i < n; i++){
        int a, b; cin>>a>>b;
        v[i] = point(a, b);
    }
    vector<bool> flag(n, false);
    vector<vector<int>> L(n);
    for(int i = 0; i < n; i++){
        if(!flag[i]){
            int ant = i;
            while(!flag[ant]){
                flag[ant] = true;
                L[ant].pb((ant + k) % n);
				L[(ant + k) % n].pb(ant);
                ant = (ant + k) % n;
            }
        }
    }
	ld res = 0;
	for(int i = 0; i < n; i++){
		point a = v[i], b, c = v[(i + 1) % n], d;
		if((v[L[i][0]] - v[i]).cross(v[L[i][1]] - v[i]) < 0) b = v[L[i][1]];
		else b = v[L[i][0]];
		if((v[L[(i + 1) % n][0]] - v[(i + 1) % n]).cross(v[L[(i + 1) % n][1]] - v[(i + 1) % n]) < 0) d = v[L[(i + 1) % n][0]];
		else d = v[L[(i + 1) % n][1]];
		point e = intersectLines(a, b - a, c, d - c);
		vector<point> r; r.pb(a); r.pb(e); r.pb(c);
		auto u = f(r, v, L);
		res += area(u);
	}
	res /= area(v);
	cout<<setprecision(25)<<(ld) 1.0 - res;
}
