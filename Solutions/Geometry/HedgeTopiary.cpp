/*
  Problem: https://codeforces.com/gym/105446/problem/H
  Explication: https://codeforces.com/blog/entry/136940
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
 
using ld = double;
const ld eps = 1e-9, inf = numeric_limits<ld>::max(), pi = acos(-1);
bool geq(ld a, ld b){return a-b >= -eps;}
bool leq(ld a, ld b){return b-a >= -eps;}
bool ge(ld a, ld b){return a-b > eps;}
bool le(ld a, ld b){return b-a > eps;}
bool eq(ld a, ld b){return abs(a-b) <= eps;}
bool neq(ld a, ld b){return abs(a-b) > eps;}

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

point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
	ld det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
}

bool pointInLine(const point & a, const point & v, const point & p){
	return eq((p - a).cross(v), 0);
}

bool pointInSegment(const point & a, const point & b, const point & p){
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

pair<ld, ld> calc(point p0, point p1, point a, point b, point c){
    if(ge(p0.cross(p1), 0)) swap(p0, p1);
    vector<ld> v;
    if(neq(p0.cross(p1), 0)){
        int u = intersectLineSegmentInfo(point(0, 0), p0, a, b);
        int w = intersectLineSegmentInfo(point(0, 0), p1, a, b);    
        if(u == 1){
            auto it = intersectLines(point(0, 0), p0, a, b - a);
            point p2 = p0 * -1;
            if(ge((p2 - it).length(), (p0 - it).length())) v.pb(it.length() / p0.length());
        }
        if(w == 1){
            auto it = intersectLines(point(0, 0), p1, a, b - a);
            point p2 = p1 * -1;
            if(ge((p2 - it).length(), (p1 - it).length())) v.pb(it.length() / p1.length());
        }
        if(le(p0.cross(a), 0) && ge(p1.cross(a), 0)){
            auto it = intersectLines(point(0, 0), p0, a, p1 - p0);
            v.pb(it.length() / p0.length());
        }
        if(le(p0.cross(b), 0) && ge(p1.cross(b), 0)){
            auto it = intersectLines(point(0, 0), p0, b, p1 - p0);
            v.pb(it.length() / p0.length());
        }
    }else{
        int u = intersectLineSegmentInfo(point(0, 0), p0, a, b);
        if(u == 1){
            auto it = intersectLines(point(0, 0), p0, a, b - a);
            point p2 = p0 * -1;
            if(ge((p2 - it).length(), (p0 - it).length())){
                if(it == b){
                    if(le(b.cross(c), 0)){
                        v.pb(it.length() / p0.length());
                        v.pb(it.length() / p1.length());
                    }    
                }else if(it != a && it != b){
                    v.pb(it.length() / p0.length());
                    v.pb(it.length() / p1.length());
                }
            }
        }
    }
    if(v.size() < 2) return {1, 0};
    sort(all(v));
    return {v[0], v[1]};
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

bool checkSegment(vector<point> B, point a, point b, ld k){
    point p0 = a * k, p1 = b * k;
    int n = B.size();
    for(int i = 0; i < n; i++){
        int l = (i + 1) % n;
        if(intersectSegmentsInfo(B[i], B[l], p0, p1) == 1){
            if(pointInSegment(B[i], B[l], p0)) continue;
            if(pointInSegment(B[i], B[l], p1)) continue;
            if(pointInSegment(p0, p1, B[i])) continue;
            if(pointInSegment(p0, p1, B[l])) continue;
            return false;
        }
    }
    if(pointInPolygon(B, p0)) return true;
    return false;
}

ld area(vector<point> & P){
	int n = P.size();
	ld ans = 0;
	for(int i = 0; i < n; i++){
		ans += P[i].cross(P[(i + 1) % n]);
	}
	return (ans / 2);
}

int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int n; cin>>n;
    vector<point> A(n);
    for(int i = 0; i < n; i++){
        int a, b; cin>>a>>b; A[i] = point(a, b);
    }
    if(le(area(A), 0)) reverse(all(A));
    int m; cin>>m;
    vector<point> B(m);
    for(int i = 0; i < m; i++){
        int a, b; cin>>a>>b; B[i] = point(a, b);
    }
    if(le(area(B), 0)) reverse(all(B));
    vector<vector<pair<ld, ld>>> L(n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            auto u = calc(A[i], A[(i + 1) % n], B[j], B[(j + 1) % m], B[(j + 2) % m]);
            if(leq(u.fi, u.se)) L[i].pb(u);
        }
    }
    for(int i = 0; i < n; i++){
        if(neq(A[i].cross(A[(i + 1) % n]), 0)) continue;
        for(int j = 0; j < m; j++){
            point p2 = A[i] * -1;
            if(le((p2 - B[j]).length(), (A[i] - B[j]).length())) continue;
            int u = intersectLineSegmentInfo(A[i], A[(i + 1) % n] - A[i], B[j], B[(j + 1) % m]);
            if(u == 1){
                point c, d, e, f;
                if(le(A[i].length(), A[(i + 1) % n].length())) c = A[i], d = A[(i + 1) % n];
                else c = A[(i + 1) % n], d = A[i];
                if(le(B[j].length(), B[(j + 1) % m].length())) e = B[j], f = B[(j + 1) % m];
                else e = B[(j + 1) % m], f = B[j];
                L[i].pb({e.length() / c.length(), e.length() / c.length()});
                L[i].pb({f.length() / d.length(), f.length() / d.length()});
                L[i].pb({e.length() / d.length(), e.length() / d.length()});
                L[i].pb({f.length() / c.length(), f.length() / c.length()});
            }
        }
    }
    for(int i = 0; i < n; i++) sort(all(L[i]));
    for(int i = 0; i < n; i++) if(L[i].size()) L[i].pb({L[i][L[i].size() - 1].se, Inf});
    vector<vector<pair<ld, ld>>> W(n);
    for(int i = 0; i < n; i++){
        if(!L[i].size()) continue;
        if(L[i].size() == 1){
            W[i].pb(L[i][0]);
            continue;
        }
        ld ant = L[i][0].fi, sec = L[i][0].se;
        if(!checkSegment(B, A[i], A[(i + 1) % n], ant)){
            ant -= eps * 3;
        }
        for(int j = 1; j < L[i].size(); j++){
            if(le(L[i][j].fi, sec)) sec = max(sec, L[i][j].se);
            else if(eq(L[i][j].fi, sec)){
                if(checkSegment(B, A[i], A[(i + 1) % n], L[i][j].fi)){
                    W[i].pb({L[i][j].fi, L[i][j].fi});
                    W[i].pb({ant, sec});
                    ant = L[i][j].fi; sec = L[i][j].se;
                }else{
                    sec = max(sec, L[i][j].se);
                }
            }else{
            	if(!checkSegment(B, A[i], A[(i + 1) % n], ant)) ant -= eps * 3;
                if(checkSegment(B, A[i], A[(i + 1) % n], ant)) W[i].pb({ant, ant});
                W[i].pb({ant, sec});
                ant = L[i][j].fi; sec = L[i][j].se;
            }
        }
        if(!checkSegment(B, A[i], A[(i + 1) % n], ant)) ant -= eps * 3;
        if(checkSegment(B, A[i], A[(i + 1) % n], ant)) W[i].pb({ant, ant});
        W[i].pb({ant, sec});
    }
    vector<pair<ld, ld>> T;
    for(auto e : W) for(auto d : e) T.pb(d);
    sort(all(T));
    ld ans = T[0].fi;
    ld at = T[0].se;
    for(int i = 0; i < T.size() - 1; i++){
        if(at == Inf) break;
        if(geq(T[i].fi, at) && T[i].se != Inf) ans = max(ans, T[i].fi);
        at = max(at, T[i].se);
    }
    cout<<setprecision(25)<<ans;
}
