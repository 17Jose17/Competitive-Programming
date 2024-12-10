/*
  Problem: https://codeforces.com/contest/107/problem/E
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

ld areaUnionTriangles(vector<pair<point, pair<point, point>>> & P){
    int n = P.size();
    vector<line> segments(n * 3); vector<char> segmentsType(n * 3);
    for(int i = 0; i < n; i++){
        point a = P[i].fi, b = P[i].se.fi, c = P[i].se.se;
        segments[i * 3] = lexCompare(a, b) ? line(a, b) : line(b, a);
        segmentsType[i * 3] = segmentType(segments[i * 3], c);
        segments[i * 3 + 1] = lexCompare(b, c) ? line(b, c) : line(c, b);
        segmentsType[i * 3 + 1] = segmentType(segments[i * 3 + 1], a);
        segments[i * 3 + 2] = lexCompare(c, a) ? line(c, a) : line(a, c);
        segmentsType[i * 3 + 2] = segmentType(segments[i * 3 + 2], b);
    }
    vector<ld> k(n * 3), t(n * 3);
    for(int i = 0; i < n * 3; i++){
        if(segmentsType[i]){
            k[i] = (segments[i].v.y - segments[i].a.y) / (segments[i].v.x - segments[i].a.x);
            t[i] = segments[i].a.y - k[i] * segments[i].a.x;
        }
    }
    ld ans = 0;
    for(int i = 0; i < 3 * n; i++){
        if(!segmentsType[i]) continue;
        ld l = segments[i].a.x, r = segments[i].v.x;
        vector<pair<ld, int>> evs;
        for(int j = 0; j < 3 * n; j++){
            if(!segmentsType[j] || i == j) continue;
            ld l1 = segments[j].a.x, r1 = segments[j].v.x;
            if(geq(l1, r) || geq(l, r1)) continue;
            ld comml = max(l, l1), commr = min(r, r1);
            auto f = intersectSegmentsInfo(segments[i].a, segments[i].v, segments[j].a, segments[j].v);
            if(f == 0){
                ld yl1 = k[j] * comml + t[j], yl = k[i] * comml + t[i];
                if(le(yl1, yl) == (segmentsType[i] == 1)){
                    int evTy = -segmentsType[i] * segmentsType[j];
                    evs.emplace_back(comml, evTy);
                    evs.emplace_back(commr, -evTy);
                }
            }else if(f == 1){
                auto u = intersectLines(segments[i].a, segments[i].v - segments[i].a, segments[j].a, segments[j].v - segments[j].a);
                ld yl = k[i] * comml + t[i], yl1 = k[j] * comml + t[j];
                int evTy = -segmentsType[i] * segmentsType[j];
                if(le(yl1, yl) == (segmentsType[i] == 1)){
                    evs.emplace_back(comml, evTy);
                    evs.emplace_back(u.x, -evTy);
                }
                yl = k[i] * commr + t[i], yl1 = k[j] * commr + t[j];
                if(le(yl1, yl) == (segmentsType[i] == 1)){
                    evs.emplace_back(u.x, evTy);
                    evs.emplace_back(commr, -evTy);
                }
            }else{
                if(segmentsType[i] != segmentsType[j] || j > i){
                    evs.emplace_back(comml, -2);
                    evs.emplace_back(commr, 2);
                }
            }
        }
        evs.emplace_back(l, 0);
        sort(all(evs));
        int j = 0, balance = 0;
        while(j < (int)evs.size()){
            int ptr = j;
            while(ptr < evs.size() && eq(evs[j].fi, evs[ptr].fi)){
                balance += evs[ptr].se;
                ptr++;
            }
            if(!balance && !eq(evs[j].fi, r)){
                ld nextx = ptr == (int)evs.size() ? r : evs[ptr].fi;
                ans -= segmentsType[i] * (k[i] * (nextx + evs[j].fi) + 2 * t[i]) * (nextx - evs[j].fi);
            }
            j = ptr;
        }
    }
    return ans / 2;
}

ld area(vector<point> & P){
	int n = P.size();
	ld ans = 0;
	for(int i = 0; i < n; i++){
		ans += P[i].cross(P[(i + 1) % n]);
	}
	return abs(ans / 2);
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	int n; cin>>n;
	vector<pair<point, pair<point, point>>> vTs;
	ld res = 0;
	for(int i = 0; i < n; i++){
	    int a, b, c, d, e, f, g, h; cin>>a>>b>>c>>d>>e>>f>>g>>h;
	    point t1(a, b), t2(c, d), t3(e, f), t4(g, h);
	    vector<point> u, u1; u.pb(t1); u.pb(t2); u.pb(t3); u1.pb(t3); u1.pb(t4); u1.pb(t1);
	    res += area(u) + area(u1);
	    vTs.pb({t1, {t2, t3}}); vTs.pb({t3, {t4, t1}});
	}
	ld ans = areaUnionTriangles(vTs);
	cout<<setprecision(20)<<res / ans;
}
