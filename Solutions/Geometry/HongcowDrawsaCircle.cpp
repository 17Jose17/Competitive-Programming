/*
  Problem: https://codeforces.com/contest/744/problem/D
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

using ld = double;
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
	//Vamos a limpiar los planos identicos
	/*for(int i = 0; i < planes.size() - 1; i++){
	    if(abs(planes[i].v.cross(planes[i + 1].a - planes[i].a)) <= eps && abs(planes[i].v.cross(planes[i + 1].v)) <= eps) cout<<"Ayuda";
	}*/
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

plane bisector(point p0, point p1){
    point pm = (p0 + p1) / 2.0;
    point b = ((p1 - p0).perp());
    //if(pm.cross(p0 - b) < -eps) b = b * -1;
    return plane(pm, b);
}

vector<point> convexHull(vector<point> P){
    sort(P.begin(), P.end());
    vector<point> L, U;
    for(int i = 0; i < P.size(); i++){
        while(L.size() >= 2 && le((L[L.size() - 2] - P[i]).cross(L[L.size() - 1] - P[i]), 0)){
            L.pop_back();
        }
        L.push_back(P[i]);
    }
    for(int i = P.size() - 1; i >= 0; i--){
        while(U.size() >= 2 && le((U[U.size() - 2] - P[i]).cross(U[U.size() - 1] - P[i]), 0)){
            U.pop_back();
        }
        U.push_back(P[i]);
    }
    L.pop_back();
    U.pop_back();
    L.insert(L.end(), U.begin(), U.end());
    return L;
}

point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
	ld det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
}
 
int intersectLineSegmentInfo(const point & a, const point & v, const point & c, const point & d){
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

int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int n, m; cin>>n>>m;
    vector<point> blue(m), red(n);
    for(int i = 0; i < n; i++){
        int a, b; cin>>a>>b; red[i] = point(a, b);
    }
    for(int i = 0; i < m; i++){
        int a, b; cin>>a>>b; blue[i] = point(a, b);
    }
    vector<point> z; for(auto e : red) z.pb(e); for(auto e : blue) z.pb(e);
    z = convexHull(z);
    for(int i = 0; i < z.size(); i++){
        for(int j = 0; j < n; j++){
            if(z[i] == red[j]){
                cout<<-1; return 0;
            }
        }
    }
    ld ans = 0;
    random_device rd;
    mt19937 g(rd());
    shuffle(all(red), g);
    srand(time(0));
    for(int i = 0; i < n; i++){
        if((double)clock() / CLOCKS_PER_SEC > 5.0) break;
        vector<plane> vp;
        for(int j = 0; j < m; j++){
            vp.pb(bisector(red[i], blue[j]));
        }
        auto u = halfPlaneIntersection(vp); int k = u.size();
        vector<point> familyxD;
        for(int j = 0; j < m; j++){
            for(int l = 0; l < k; l++){
                int l1 = l + 1; if(l1 == k) l1 = 0;
                if(eq((red[i] - (u[l] + u[l1]) / 2).length(), (blue[j] - (u[l] + u[l1]) / 2).length())){
                    familyxD.pb(blue[j]); break;
                }
            }
        }
        int t = familyxD.size();
        vector<plane> q;
        for(int j = 0; j < t; j++){
            q.pb(bisector(red[i], familyxD[j]));
        }
        for(int j = 0; j < t; j++){
            auto w = q;
            for(int l = 0; l < t; l++){
                if(j == l) continue;
                w.pb(bisector(familyxD[j], familyxD[l]));
            }
            auto p = halfPlaneIntersection(w);
            for(int l = 0; l < p.size(); l++){
                ans = max(ans, (familyxD[j] - p[l]).length());
            }
        }
    }
    cout<<setprecision(25)<<ans;
}
