#include <bits/stdc++.h>
using namespace std;
 
// Holi c:
 
#define ll long long int
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()
 
const int Inf = 1e3;
const ll mod = 1e9+7;
const ll INF = 1e18;
 
using ld = double;
const ld eps = 1e-6, inf = numeric_limits<ld>::max(), pi = acos(-1);
bool geq(ld a, ld b){return a-b >= -eps;}     //a >= b
bool leq(ld a, ld b){return b-a >= -eps;}     //a <= b
bool ge(ld a, ld b){return a-b > eps;}        //a > b
bool le(ld a, ld b){return b-a > eps;}        //a < b
bool eq(ld a, ld b){return abs(a-b) <= eps;}  //a == b
bool neq(ld a, ld b){return abs(a-b) > eps;}  //a != b
 
struct point3{
    ld x, y, z;
    point3(): x(0), y(0), z(0){}
    point3(ld x, ld y, ld z): x(x), y(y), z(z){}
 
    point3 operator+(const point3 & p) const{return point3(x + p.x, y + p.y, z + p.z);}
    point3 operator-(const point3 & p) const{return point3(x - p.x, y - p.y, z - p.z);}
    point3 operator*(const ld & k) const{return point3(x * k, y * k, z * k);}
    point3 operator/(const ld & k) const{return point3(x / k, y / k, z / k);}
 
    point3 operator+=(const point3 & p){*this = *this + p; return *this;}
    point3 operator-=(const point3 & p){*this = *this - p; return *this;}
    point3 operator*=(const ld & p){*this = *this * p; return *this;}
    point3 operator/=(const ld & p){*this = *this / p; return *this;}
 
    ld dot(const point3 & p) const{return x * p.x + y * p.y + z * p.z;}
    point3 cross(const point3 & p) const{return {y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x};}
    ld norm() const{return x * x + y * y + z * z;}
    ld length() const{return sqrtl(x * x + y * y + z * z);}
    point3 unit() const{return (*this) / length();}
    ld angle(const point3 & p) const{ld a = atan2l(((*this).cross(p)).norm(), (*this).dot(p)); if(le(a, 0)) a += pi * 2; return a;}
 
    bool operator==(const point3 & p) const{return eq(x, p.x) && eq(y, p.y) && eq(z, p.z);}
    bool operator!=(const point3 & p) const{return !(*this == p);}
    bool operator<(const point3 & p) const{return le(x, p.x) || (eq(x, p.x) && le(y, p.y)) || (eq(x, p.x) && eq(y, p.y) && le(z, p.z));}
    bool operator>(const point3 & p) const{return ge(x, p.x) || (eq(x, p.x) && ge(y, p.y)) || (eq(x, p.x) && eq(y, p.y) && ge(z, p.z));}
    bool operator<=(const point3 & p) const{return !(*this > p);}
    bool operator>=(const point3 & p) const{return !(*this < p);}
};
 
//istream &operator>>(istream &is, point3 & p){return is >> p.x >> p.y >> p.z;}
//ostream &operator<<(ostream &os, const point3 & p){return os << "(" << p.x << ", " << p.y << " , " << p.z << ")";}
 
struct plane{
    point3 n; ld d;
    plane(): n(0, 0, 0), d(0){}
    plane(point3 n, ld d): n(n), d(d){}
    plane(point3 n, point3 p): n(n), d(n.dot(p)){}
    plane(point3 p1, point3 p2, point3 p3): plane((p2 - p1).cross(p3 - p1), p1.dot((p2 - p1).cross(p3 - p1))){}
 
    ld side(const point3 & p) const{return ((*this).n).dot(p) - (*this).d;}
    ld dist(const point3 & p) const{return fabsl((*this).side(p)) / ((*this).n).length();}
    plane translate(const point3 & p) const{return {(*this).n, (*this).d + ((*this).n).dot(p)};}
    plane shift(const ld & t) const{return {(*this).n, (*this).d + t * ((*this).n).length()};}
    point3 proj(const point3 & p) const{return p - (*this).n * (*this).side(p) / ((*this).n).length();}
    point3 refl(const point3 & p) const{return p - (*this).n * (*this).side(p) / ((*this).n).length() * 2;}
    //Minium angle between two planes
    ld anglep(const plane & p) const{return ((*this).n).angle(p.n);}
    //Angle between a plane and line is pi / 2 - angle(between two vector directionals) ans parallel and perpendicular is exactly of plane to plane
 
    bool operator|(const plane & p) const{const point3 r = ((*this).n).cross(p.n); return (r.x == 0 && r.y == 0 && r.z == 0);} //Parallel
    bool operator/(const plane & p) const{return ((*this).n).dot(p.n) == 0;} //Perpendicular
 
};
 
struct line{
    point3 p, q;
    line(): p(0, 0, 0), q(0, 0, 0){}
    line(point3 p, point3 q): p(p), q(q){}
 
    ld dist(const point3 & p) const{return sqrtl((((*this).q).cross(((*this).p) - p)).norm() / ((*this).q).norm());}
    point3 proj(const point3 & p) const{return (*this).p + ((*this).q) * (((*this).q).dot(p - (*this).p)) / ((*this).q).norm();}
    point3 refl(const point3 & p) const{return (*this).proj(p) * 2 - p;}
    //assuming intersection exist and p.n and d not collinear (p.n).dot(d) != 0
    point3 inter(const plane & p) const{return (*this).p - ((*this).q) * (p.side((*this).p)) / (p.n).dot((*this).q);}
    //Minium angle between two lines
    ld angle(const line & l) const{return ((*this).q).angle(l.q);}
    //Angle between a line and plane is pi / 2 - angle(between two vector directionals) ans parallel and perpendicular is exactly of line to line
 
    bool operator|(const line & l) const{const point3 r = ((*this).q).cross(l.q); return (r.x == 0 && r.y == 0 && r.z == 0);} //Parallel
    bool operator/(const line & l) const{return ((*this).q).dot(l.q) == 0; } //Parallel
};
 
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
 
struct coords{
    point3 o, dx, dy, dz;
 
    coords(point3 p, point3 q, point3 r) : o(p) {
        dx = (q - p).unit();
        dz = ((dx).cross(r - p)).unit();
        dy = dz.cross(dx);
    }
 
    coords(point3 p, point3 q, point3 r, point3 s) : 
     o(p), dx(q - p), dy(r - p), dz(s - p) {}
    //Return a point in 2D
    point pos2d(point3 & p) {return point((p - o).dot(dx), (p - o).dot(dy));}
    point3 pos3d(point3 p) {return point3((p - o).dot(dx), (p - o).dot(dy), (p - o).dot(dz));}
};
 
int sgn(ld x){
    if(ge(x, 0)) return 1;
    if(le(x, 0)) return -1;
    return 0;
}
 
vector<pair<point, int>> convexHull(vector<pair<point, int>> & P){
	sort(P.begin(), P.end());
	vector<pair<point, int>> L, U;
	for(int i = 0; i < P.size(); i++){
		while(L.size() >= 2 && leq((L[L.size() - 2].fi - P[i].fi).cross(L[L.size() - 1].fi - P[i].fi), 0)){
			L.pop_back();
		}
		L.push_back(P[i]);
	}
	for(int i = P.size() - 1; i >= 0; i--){
		while(U.size() >= 2 && leq((U[U.size() - 2].fi - P[i].fi).cross(U[U.size() - 1].fi - P[i].fi), 0)){
			U.pop_back();
		}
		U.push_back(P[i]);
	}
	L.pop_back();
	U.pop_back();
	L.insert(L.end(), U.begin(), U.end());
	return L;
}
 
//Auxiliar points
point3 zero;
 
point3 closestPointLs(const line & l1, const line & l2){
    point3 n = (l2.q).cross((l1.q).cross(l2.q));
    return l1.p + l1.q * (l2.p - l1.p).dot(n) / (l1.q).dot(n);
}
 
int infointersectSegments(point3 & a, point3 & b, point3 & c, point3 & d){
    line l1(a, b - a), l2(c, d - c);
    if((l1.q).cross(l2.q) == zero) return 0;
    point3 it = closestPointLs(l1, l2);
    if(it == a || it == b || it == c || it == d) return -1;
    if(it > min(a, b) && it < max(a, b) && it > min(c, d) && it < max(c, d)) return 1;
    return 0;
}
 
//Assuming that they intersect and not parallel, (pl.n).cross(ql.n) != {0,0,0}
line intersectionPlanes(const plane & pl, const plane & ql){
    point3 q = (pl.n).cross(ql.n);
    point3 p = (ql.n * pl.d - pl.n * ql.d).cross(q) / q.norm();
    return line(p, q);
}
 
bool point3InLine(const point3 & a, const point3 & v, const point3 & p){
	return eq(line(a, v).dist(p), 0);
}
 
bool point3InSegment(const point3 & a, const point3 & b, const point3 & p){
    if(p == a || p == b) return 1;
	return point3InLine(a, b - a, p) && p >= min(a, b) && p <= max(a, b);
}
 
int intersectSegmentPlaneInfo(const point3 & a, const point3 & b, const plane & p){
    //if(leq(p.dist(a), 0) && leq(p.dist(b), 0)) return -1;
    point3 d = b - a;
    if(eq(d.dot(p.n), 0)){
        if(eq(p.dist(a), 0)) return -1;
        else return 0;
    }else{
        if(point3InSegment(a, b, line(a, b - a).inter(p))) return 1;
    }
    return 0;
}
 
ld volume(vector<vector<point3>> & P){
    int n = P.size();
    ld ans = 0;
    for(int i = 0; i < n; i++){
        point3 res;
        int m = P[i].size();
        for(int j = 0; j < m; j++){
            res += (P[i][j]).cross(P[i][(j + 1) % m]);
        }
        ans += (P[i][0]).dot(res);
    }
    return (ans) / 6.0;
}
 
void reorient(vector<vector<point3>> & P){
    int n = P.size();
    vector<vector<pair<int, bool>>> g(n);
    map<pair<point3, point3>, int> es;
    for(int i = 0; i < n; i++){
        int m = P[i].size();
        for(int j = 0; j < m; j++){
            int l = j + 1; if(l == m) l = 0;
            point3 a = P[i][j], b = P[i][l];
            if(es.find({a, b}) != es.end()){
                int v = es[{a, b}];
                g[i].pb({v, true});
                g[v].pb({i, true});
            }else if(es.find({b, a}) != es.end()){
                int v = es[{b, a}];
                g[i].pb({v, false});
                g[v].pb({i, false});
            }else es[{a, b}] = i;
        }
    }
    vector<bool> vis(n, false), flip(n, false);
    queue<int> q;
    q.push(0);
    while(!q.empty()){
        int u = q.front();
        q.pop();
        vis[u] = true;
        for(auto e : g[u]){
            if(!vis[e.fi]){
                vis[e.fi] = true;
                flip[e.fi] = (flip[u] ^ e.se);
                q.push(e.fi);
            }
        }
    }
    for(int i = 0; i < n; i++){
        if(flip[i]) reverse(all(P[i]));
    }
    return;
}
 
bool equalFace(vector<point3> A, vector<point3> B){
    int n = A.size(), m = B.size();
    sort(all(A)); sort(all(B));
    if(n != m) return false;
    for(int i = 0; i < n; i++) if(A[i] != B[i]) return false;
    return true;
}
 
int vls(bool a, bool b){
    if(a && b) return 0;
    if(!a && !b) return 2;
    if(a) return 1;
    return -1;
}
 
int sidePolyhedroPlane(vector<vector<point3>> P, plane p){
    int n = P.size(); bool sup = false, infr = false;
    for(int i = 0; i < n; i++){
        for(int l = 0; l < P[i].size(); l++){
            auto u = p.side(P[i][l]);
            if(ge(u, 0)) sup = true;
            if(le(u, 0)) infr = true;
        }
    }
    return vls(sup, infr);
}
 
int sideFacePlane(vector<point3> F, plane p){
    int n = F.size(); bool sup = false, infr = false;
    for(int i = 0; i < n; i++){
        auto u = p.side(F[i]);
        if(ge(u, 0)) sup = true;
        if(le(u, 0)) infr = true;
    }
    return vls(sup, infr);
}
 
vector<point3> clear(vector<point3> P){
    int n = P.size();
    if(n <= 2) return {};
    sort(all(P)); P.erase(unique(all(P)), P.end());
    n = P.size();
    if(n <= 2) return {};
    vector<pair<point, int>> aux;
    int j = 0;
    while((P[1] - P[0]).cross(P[j] - P[0]) == zero && j < n) j++;
    if(j == n) return {};
    coords base(P[0], P[1], P[j]);
    for(int i = 0; i < n; i++){
        aux.pb({base.pos2d(P[i]), i});
    }
    aux = convexHull(aux);
    if(aux.size() < 3) return {};
    vector<point3> res;
    for(int i = 0; i < aux.size(); i++){
        res.pb(P[aux[i].se]);
    }
    return res;
}
 
vector<vector<point3>> cutConvexPolyhedro(vector<vector<point3>> & P, plane p){
    int n = P.size();
    //if(sidePolyhedroPlane(P, p) > 0) return P;
    if(n < 4) return {};
    vector<vector<point3>> ans;
    vector<point3> at;
    for(int i = 0; i < n; i++){
        int s = sideFacePlane(P[i], p);
        if(s != 0){
            if(s > 0) ans.pb(P[i]); 
            if(s == -1) for(int j = 0; j < P[i].size(); j++) if(eq(p.side(P[i][j]), 0)) at.pb(P[i][j]);
            continue;
        }
        int m = P[i].size();
        vector<point3> res;
        for(int j = 0; j < m; j++){
            int l = j + 1; if(l == m) l = 0;
            if(geq(p.side(P[i][j]), 0)) res.pb(P[i][j]);
            if(eq(p.side(P[i][j]), 0)) at.pb(P[i][j]);
            if(intersectSegmentPlaneInfo(P[i][j], P[i][l], p) == 1){
                auto it = line(P[i][j], P[i][l] - P[i][j]).inter(p);
                if(it != P[i][j] && it != P[i][l]) res.pb(it);
                at.pb(it);
            }
        }
        if((int) res.size() > 2) ans.pb(res);
    }
    at = clear(at);
    if((int) at.size() > 2) ans.pb(at);
    if((int) ans.size() < 4) return {};
    return ans;
}
 
vector<vector<point3>> inicialize(){
    vector<vector<point3>> faces;
    faces.pb({point3(-Inf, -Inf, -Inf), point3(Inf, -Inf, -Inf), point3(Inf, Inf, -Inf), point3(-Inf, Inf, -Inf)});
    faces.pb({point3(-Inf, -Inf, -Inf), point3(-Inf, -Inf, Inf), point3(-Inf, Inf, Inf), point3(-Inf, Inf, -Inf)});
    faces.pb({point3(-Inf, -Inf, -Inf), point3(Inf, -Inf, -Inf), point3(Inf, -Inf, Inf), point3(-Inf, -Inf, Inf)});
    faces.pb({point3(-Inf, Inf, -Inf), point3(Inf, Inf, -Inf), point3(Inf, Inf, Inf), point3(-Inf, Inf, Inf)});
    faces.pb({point3(Inf, -Inf, -Inf), point3(Inf, Inf, -Inf), point3(Inf, Inf, Inf), point3(Inf, -Inf, Inf)});
    faces.pb({point3(-Inf, -Inf, Inf), point3(Inf, -Inf, Inf), point3(Inf, Inf, Inf), point3(-Inf, Inf, Inf)});
    return faces;
}
 
vector<vector<point3>> cubo(point3 i0, int d){
    point3 i1(i0.x + d, i0.y, i0.z), i2(i0.x + d, i0.y + d, i0.z), i3(i0.x, i0.y + d, i0.z), i4(i0.x, i0.y, i0.z + d),
     i5(i0.x + d, i0.y, i0.z + d), i6(i0.x + d, i0.y + d, i0.z + d), i7(i0.x, i0.y + d, i0.z + d);
    vector<vector<point3>> ans;
    ans.pb({i0, i1, i2, i3}); ans.pb({i4, i5, i6, i7}); ans.pb({i0, i1, i5, i4}); 
    ans.pb({i3, i2, i6, i7}); ans.pb({i0, i3, i7, i4}); ans.pb({i1, i2, i6, i5});
    reorient(ans); if(ge(volume(ans), 0)) for(auto & e : ans) reverse(all(e));
    return ans;
}
 
plane divSpace(point3 a, point3 b){
    point3 m = (a + b) / 2;
    return plane(a - b, m);
}
 
vector<vector<point3>> voronoi(vector<point3> P, int j, int f){
    auto ans = inicialize();
    for(int i = 0; i < P.size(); i++){
        if(i == j) continue;
        auto p = divSpace(P[j], P[i]);
        if(f){ p.n *= -1; p.d *= -1; }
        ans = cutConvexPolyhedro(ans, p);
    }
    return ans;
}
 
vector<vector<point3>> intersectionConvexPolyehdros(vector<vector<point3>> A, vector<vector<point3>> B){
    int n = B.size();
    for(int i = 0; i < n; i++){
        A = cutConvexPolyhedro(A, plane(B[i][0], B[i][1], B[i][2]));
        if(A.size() < 4) return {};
    }
    return A;
}
 
int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int n, q; cin>>n>>q;
    vector<point3> points(n);
    for(int i = 0; i < n; i++){
        int a, b, c; cin>>a>>b>>c;
        points[i] = point3(a, b, c);
    }
    vector<vector<vector<point3>>> v(n);
    vector<vector<vector<point3>>> iv(n);
    for(int i = 0; i < n; i++){
        v[i] = voronoi(points, i, 0);
    }
    for(int i = 0; i < n; i++){
        iv[i] = voronoi(points, i, 1);
    }
    vector<vector<ld>> cb;
    for(int i = 0; i < q; i++){
        int a, b, c, d; cin>>a>>b>>c>>d; cb.pb({a, b, c, d});
    }
    vector<ld> sum(n), res(n);
    for(int i = 0; i < q; i++){
        ld a = cb[i][0], b = cb[i][1], c = cb[i][2], d = cb[i][3];
        auto cube = cubo(point3(a, b, c), d);
        ld v0 = (ld) d * d * d;
        for(int j = 0; j < n; j++){
            auto t = intersectionConvexPolyehdros(v[j], cube);
            if(t.size() > 3){
                reorient(t); sum[j] += abs(volume(t)) / v0;
            }
        }
        for(int j = 0; j < n; j++){
            auto t = intersectionConvexPolyehdros(iv[j], cube);
            if(t.size() > 3){
                reorient(t); res[j] += abs(volume(t)) / v0;
            }
        }
    }
    vector<vector<vector<vector<point3>>>> v2(n);
    for(int i = 0; i < n; i++){
        vector<int> s;
        int k = v[i].size();
        for(int j = 0; j < n; j++){
            if(i == j) continue;
            for(int l = 0; l < k; l++){
                if(eq(plane(v[i][l][0], v[i][l][1], v[i][l][2]).dist((points[i] + points[j]) / 2), 0)){
                    s.pb(j); break;
                }
            }
        }
        int t = s.size();
        for(int j = 0; j < t; j++){
            auto r = v[i];
            for(int l = 0; l < t; l++){
                if(j == l) continue;
                r = cutConvexPolyhedro(r, divSpace(points[s[j]], points[s[l]]));
            }
            v2[s[j]].pb(r);
        }
    }
    
    vector<vector<vector<vector<point3>>>> iv2(n);
    for(int i = 0; i < n; i++){
        if(!iv[i].size()) continue;
        for(int j = 0; j < n; j++){
            if(i == j) continue;
            auto r = iv[i];
            for(int l = 0; l < n; l++){
                if(j == l || l == i) continue;
                r = cutConvexPolyhedro(r, divSpace(points[l], points[j]));
            }
            if(r.size()) iv2[j].pb(r);
        }
    }
    
    for(int i = 0; i < q; i++){
        ld a = cb[i][0], b = cb[i][1], c = cb[i][2], d = cb[i][3];
        auto cube = cubo(point3(a, b, c), d);
        ld v0 = (ld) d * d * d;
        for(int j = 0; j < n; j++){
            for(int l = 0; l < v2[j].size(); l++){
                auto t = intersectionConvexPolyehdros(v2[j][l], cube);
                if(t.size() < 4) continue;
                reorient(t);
                sum[j] += abs(volume(t)) / v0 * 2;
            }
        }
    }
    
    for(int i = 0; i < q; i++){
        ld a = cb[i][0], b = cb[i][1], c = cb[i][2], d = cb[i][3];
        auto cube = cubo(point3(a, b, c), d);
        ld v0 = (ld) d * d * d;
        for(int j = 0; j < n; j++){
            for(int l = 0; l < iv2[j].size(); l++){
                auto t = intersectionConvexPolyehdros(iv2[j][l], cube);
                if(t.size() < 4) continue;
                reorient(t);
                res[j] += abs(volume(t)) / v0 * 2;
            }
        }
    }
    
    /*ld e = 0, f = 0;
    for(auto g : sum) e += g; for(auto g : res) f += g;
    cout<<e<<" "<<f<<'\n';*/
    for(int i = 0; i < n; i++) cout<<setprecision(9)<<sum[i] - res[i]<<'\n';
}
