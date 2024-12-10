/*
  Problem: https://qoj.ac/contest/442/problem/1196
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

const int Inf = 1e5;
const ll mod = 1e9+7;
const ll INF = 1e16;

using ld = double;
const ld eps = 1e-9, inf = numeric_limits<ld>::max(), pi = acos(-1);
bool geq(ld a, ld b){return a-b >= -eps;}
bool leq(ld a, ld b){return b-a >= -eps;}
bool ge(ld a, ld b){return a-b > eps;}
bool le(ld a, ld b){return b-a > eps;}
bool eq(ld a, ld b){return abs(a-b) <= eps;}
bool neq(ld a, ld b){return abs(a-b) > eps;}

struct point {
    ld x, y;
    point(): x(0), y(0) {}
    point(ld x, ld y): x(x), y(y) {}

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
    ld ang() const {
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

void polarSort(vector<pair<point, int>> & P, const point & o, const point & v){
	//sort points in P around o, taking the direction of v as first angle
	sort(P.begin(), P.end(), [&](const pair<point, int> & a, const pair<point, int> & b){
		return point((a.fi - o).half(v), 0) < point((b.fi - o).half(v), (a.fi - o).cross(b.fi - o));
	});
}

bool pointInLine(const point & a, const point & v, const point & p){
    return eq((p - a).cross(v), 0);
}

bool pointInSegment(const point & a, const point & b, const point & p){
    return pointInLine(a, b - a, p) && leq((a - p).dot(b - p), 0);
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
            return -1; // infinity points
        } else {
            return 0; // no point
        }
    } else {
        return sgn(v.cross(c - a)) != sgn(v.cross(d - a)); // 1: single point, 0: no point
    }
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

int intersectSegmentsInfo(const point &a, const point &b, const point &c, const point &d) {
    point v1 = b - a, v2 = d - c;
    int t = sgn(v1.cross(c - a)), u = sgn(v1.cross(d - a));
    if (t == u) {
        if (t == 0) {
            if (pointInSegment(a, b, c) || pointInSegment(a, b, d) || pointInSegment(c, d, a) || pointInSegment(c, d, b)) {
                return -1;
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    } else {
        return sgn(v2.cross(a - c)) != sgn(v2.cross(b - c));
    }
}

pair<vector<pair<point, point>>, vector<point>> precFunPolygon(vector<point> P) {
    int n = P.size();
    vector<point> prov;
    vector<pair<point, point>> Lprov;
    for (int i = 0; i < n; i++) {
        if (geq((P[(i + 1) % n] - P[i]).cross(P[(i + 2) % n] - P[i]), 0)) {
            prov.pb(P[(i + 1) % n]); prov.pb(P[(i + 2) % n]);
            Lprov.pb({P[(i + 1) % n], P[(i + 2) % n]});
        } else {
            int l = -1;
            point at(Inf, Inf), seg;
            for (int j = 0; j < n; j++) {
                if (j == i || j == ((i + 1) % n) || ((j + 1) % n) == i || ((j + 1) % n) == ((i + 1) % n)) continue;
                auto u = intersectLineSegmentInfo(P[i], P[(i + 1) % n] - P[i], P[j], P[(j + 1) % n]);
                if (u == 1) {
                    auto v = intersectLines(P[i], P[(i + 1) % n] - P[i], P[j], P[(j + 1) % n] - P[j]);
                    if (le((P[i] - v).length(), (P[(i + 1) % n] - v).length())) continue;
                    if (v == P[(j + 1) % n]) continue;
                    if (v == P[j] && le((P[(i + 1) % n] - v).length(), (P[(i + 1) % n] - at).length())) {
                        if (ge((v - P[(i + 1) % n]).cross(P[(j - 1 + n) % n] - P[(i + 1) % n]), 0)) at = P[j], seg = P[(j - 1 + n) % n], l = j;
                        if (ge((v - P[(i + 1) % n]).cross(P[(j + 1) % n] - P[(i + 1) % n]), 0)) at = P[j], seg = P[(j + 1) % n], l = j;
                    } else if (le((P[(i + 1) % n] - v).length(), (P[(i + 1) % n] - at).length())) {
                        if (ge((v - P[(i + 1) % n]).cross(P[j] - P[(i + 1) % n]), 0)) at = v, seg = P[j];
                        if (ge((v - P[(i + 1) % n]).cross(P[(j + 1) % n] - P[(i + 1) % n]), 0)) at = v, seg = P[(j + 1) % n];
                    }
                }
            }
            prov.pb(P[(i + 1) % n]); prov.pb(at); if(l == -1) prov.pb(seg); else if(at != P[l]) prov.pb(seg);
            Lprov.pb({P[(i + 1) % n], at}); if(l == -1) Lprov.pb({at, seg}); else if(at != P[l]) Lprov.pb({at, seg});
        }
    }
    sort(all(prov));
    prov.erase(unique(all(prov)), prov.end());
    sort(all(Lprov));
    Lprov.erase(unique(all(Lprov)), Lprov.end());
    return {Lprov, prov};
}

pair<vector<vector<int>>, vector<point>> precFunPolygon1(vector<pair<point, point>> L, vector<point> P, vector<point> P0){
	int n = L.size();
	map<point, point> mp;
	map<point, vector<point>> mps;
	vector<pair<point, point>> Lprov;
	vector<point> prov;
	for(int i = 0; i < n; i++){
		point at = L[i].se;
		int l1 = -1;
		for(int j = 0; j < n; j++){
			if(L[i].fi == L[j].fi || L[i].fi == L[j].se || L[i].se == L[j].fi || L[i].se == L[j].se) continue;
			if(intersectSegmentsInfo(L[i].fi, L[i].se, L[j].fi, L[j].se) == 1){
				if(leq((L[j].se - L[j].fi).cross(L[i].fi - L[j].fi), 0)) continue;
				auto it = intersectLines(L[i].fi, L[i].se - L[i].fi, L[j].fi, L[j].se - L[j].fi);
				if(it == at) continue;
				if(le((it - L[i].fi).length(), (at - L[i].fi).length())) at = it, l1 = j;
			}
		}
		if(at != L[i].se){
			if(at != L[i].fi) Lprov.pb({L[i].fi, at});
			mps[L[l1].fi].pb(at);
			mp[L[i].fi] = at;	
	    	prov.pb(at); prov.pb(L[i].fi);
		}else{
		    mp[L[i].fi] = L[i].se;
		    prov.pb(L[i].fi); prov.pb(L[i].se);
			Lprov.pb({L[i].fi, L[i].se});
		}
	}
	for(auto e : mps){
		auto at = mp[e.fi];
		for(auto d : e.se){
			Lprov.pb({d, at});
		}
	}
	sort(all(prov));
	prov.erase(unique(all(prov)), prov.end());
	int k = prov.size();
	vector<vector<int>> Lf(k);
	for(int i = 0; i < Lprov.size(); i++){
		int i1 = lower_bound(all(prov), Lprov[i].fi) - prov.begin(), i2 = lower_bound(all(prov), Lprov[i].se) - prov.begin();
		if(i1 != i2) Lf[i1].pb(i2);
	}
	return {Lf, prov};
}

vector<point> funPolygon(vector<vector<int>> L, vector<point> P, int ini){
	int n = P.size();
	vector<int> res;
	vector<point> ans;
	vector<bool> fls(n, false);
	stack<int> q;
	q.push(ini);
	while(q.size()){
		int v = q.top();
		res.pb(v);
		q.pop();
		if(fls[v]){
			bool fl = false;
			for(int i = 0; i < res.size(); i++){
				if(res[i] == v) fl = true;
				if(fl) ans.pb(P[res[i]]);
			}
			break;
		}
		fls[v] = true;
		int u = L[v][0];
		q.push(u);
	}
	ans = convexHull(ans);
	return ans;
}

int main(){
	ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
	int n; cin>>n;
	vector<point> P(n);
	for(int i = 0; i < n; i++){
		int a, b; cin>>a>>b;
		P[i] = point(a, b);
	}
	vector<vector<int>> L0(n);
	auto vx = P;
	sort(all(vx));
	for(int i = 0; i < n; i++){
		L0[lower_bound(all(vx), P[i]) - vx.begin()].pb(lower_bound(all(vx), P[(i + 1) % n]) - vx.begin());
	}
	auto u = precFunPolygon(P);
	auto v = precFunPolygon1(u.fi, u.se, P);
	vector<vector<point>> Ps;
	for(int i = 0; i < v.se.size(); i++){
		auto t = funPolygon(v.fi, v.se, i);
		Ps.pb(t);
	}
	sort(all(Ps));
	Ps.erase(unique(all(Ps)), Ps.end());
	ld ans = 0;
	if(Ps.size()) ans = area(Ps[0]);
	if(Ps.size() > 1) ans = 0;
    cout<<setprecision(25)<<ans;
}
