#include <bits/stdc++.h>
using namespace std;
 
// Holi c:
 
#define ll long long int
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()
 
const int Inf = 1e9;
const ll mod = 1e9+7;
const ll INF = 1e18;
 
using ld = long double;
const ld eps = 1e-9, inf = numeric_limits<ld>::max(), pi = acos(-1);
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
    ld angle(const point3 & p) const{return acos(max((ld)-1.0, min((ld)1.0, (*this).dot(p) / (*this).length() / p.length())));}
 
    bool operator==(const point3 & p) const{return eq(x, p.x) && eq(y, p.y) && eq(z, p.z);}
    bool operator!=(const point3 & p) const{return !(*this == p);}
    bool operator<(const point3 & p) const{return le(x, p.x) || (eq(x, p.x) && le(y, p.y)) || (eq(x, p.x) && eq(y, p.y) && le(z, p.z));}
    bool operator>(const point3 & p) const{return ge(x, p.x) || (eq(x, p.x) && ge(y, p.y)) || (eq(x, p.x) && eq(y, p.y) && ge(z, p.z));}
};
 
istream &operator>>(istream &is, point3 & p){return is >> p.x >> p.y >> p.z;}
ostream &operator<<(ostream &os, const point3 & p){return os << "(" << p.x << ", " << p.y << " , " << p.z << ")";}

struct plane{
    point3 n; ld d;
    plane(): n(0, 0, 0), d(0){}
    plane(point3 n, ld d): n(n), d(d){}
    plane(point3 p1, point3 p2, point3 p3): plane((p2 - p1).cross(p3 - p1), p1.dot((p2 - p1).cross(p3 - p1))){}
    
    ld side(const point3 & p) const{return ((*this).n).dot(p) - (*this).d;}
    ld dist(const point3 & p) const{return abs((*this).side(p)) / ((*this).n).length();}
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

point3 zero(0, 0, 0);

int sgn(ld x){
    if(ge(x, 0)) return 1;
    if(le(x, 0)) return -1;
    return 0;
}

line perpPlane(const plane & pl, const point3 & p){
    return line(p, pl.n);
}

plane perpLine(const line & l, const point3 & p){
    return plane(l.q, (l.q).dot(p));
}

//Assuming that they intersect and not parallel, (pl.n).cross(ql.n) != {0,0,0}
line intersectionPlanes(const plane & pl, const plane & ql){
    point3 q = (pl.n).cross(ql.n);
    point3 p = (ql.n * pl.d - pl.n * ql.d).cross(q) / q.norm();
    return line(p, q);
}

ld distanceLines(const line & l1, const line & l2){
    point3 n = (l1.q).cross(l2.q);
    if(n == zero)
        return l1.dist(l2.p);
    return abs((l1.p - l2.p).dot(n)) / n.length();
}

//Closest point to l2 at l1
point3 closestPointLs(const line & l1, const line & l2){
    point3 n = (l2.q).cross((l1.q).cross(l2.q));
    return l1.p + l1.q * (l2.p - l1.p).dot(n) / (l1.q).dot(n);
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
    point pos2d(point3 p) {return point((p - o).dot(dx), (p - o).dot(dy));}
    point3 pos3d(point3 p) {return point3((p - o).dot(dx), (p - o).dot(dy), (p - o).dot(dz));}
};

// lim = 5000, p = 0.1 & pass = 0.995 works okey
pair<point3, ld> miniumSphereEnclosing(vector<point3> P, const point3 & init, const int & lim, const & ld p, const ld & pass){
    int n = P.size();
    ld ans = INF;
    for(int i = 0; i < lim; i++){
        point3 res = P[0];
        for(int j = 1; j < n; j++) if((res - init).length() < (P[j] - init).length()) res = P[j];
        ans = min(ans, (res - init).length());
        init += (res - init) * p;
        p *= pass;
    }
    return {init, ans};
}

ld areaFace(vector<point3> P){
    int n = P.size();
    point3 ans;
    for(int i = 0; i < n; i++){
        ans += (P[i]).cross(P[(i + 1) % n]); 
    }
    return (ans.length()) / 2.0;
}

ld volume(vector<vector<point3>> P){
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
    return abs(ans) / 6.0;
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

plane divSpace(point3 a, point3 b){
    point3 m = (a + b) / 2;
    return plane(a - b, m);
}

//Runs in O(n^4), too slow, easy implementation
vector<vector<point3>> convexHull4n(vector<point3> P){
    int n = P.size();
    vector<vector<point3>> ans;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            for(int l = j + 1; l < n; l++){
                bool flag = true; int f = 0;
                plane pl(P[i], P[j], P[l]);
                if(pl.n == zero) continue;
                for(int l1 = 0; l1 < n; l1++){
                    if(l1 == i || l1 == j || l1 == l) continue;
                    if(!f) f = sgn(pl.side(P[l1]));
                    else if(f != sgn(pl.side(P[l1]))) flag = false;
                    if(!flag) break;
                }
                if(flag) ans.pb({P[i], P[j], P[l]});
            }
        }
    }
    return ans;
}

//Runs in O(n^2)
vector<vector<point3>> convexHull2n(vector<point3> P){
    int n = P.size();
    vector<vector<point3>> ans;
    vector<pair<int, pair<int, pair<int, point3>>>> res;
    vector<vector<bool>> flags(n, vector<bool>(n, true));
    auto add = [&](int i, int j, int l){
        res.pb({i, {j, {l, (P[j] - P[i]).cross(P[l] - P[i])}}});
        flags[i][j] = flags[j][l] = flags[l][i] = false;
    };
    add(0, 1, 2); add(0, 2, 1);
    for(int i = 3; i < n; i++){
        vector<pair<int, pair<int, pair<int, point3>>>> sup_res;
        for(auto & f : res){
            if((P[i] - P[f.fi]).dot(f.se.se.se) > 0)
                flags[f.fi][f.se.fi] = flags[f.se.fi][f.se.se.fi] = flags[f.se.se.fi][f.fi] = true;
        else
            sup_res.pb(f);
        }
        res.clear();
        for(auto & f : sup_res){
            if(flags[f.se.fi][f.fi]) add(f.se.fi, f.fi, i);
            if(flags[f.se.se.fi][f.se.fi]) add(f.se.se.fi, f.se.fi, i);
            if(flags[f.fi][f.se.se.fi]) add(f.fi, f.se.se.fi, i);
        }
        res.insert(res.end(), all(sup_res));
    }
    for(auto & i : res){
        ans.pb({P[i.fi], P[i.se.fi], P[i.se.se.fi]});
    }
    return ans;
}

int vls(bool a, bool b){
    if(a && b) return 0; // Intersection
    if(!a && !b) return 2; // Inside
    if(a) return 1; // Above
    return -1; // Below
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

vector<point3> createConvexFace(vector<point3> P){
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
    at = createConvexFace(at);
    if((int) at.size() > 2) ans.pb(at);
    if((int) ans.size() < 4) return {};
    return ans;
}

/*Coming soon...
    ConvexHull in nlogn
    Spherical Geometry
    ¿Half spaces? :O
*/
