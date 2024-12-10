/*
  Problem: https://codeforces.com/gym/101673
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
const ld eps = 1e-9, pi = acos(-1);
bool geq(ld a, ld b){return a - b >= -eps;}
bool leq(ld a, ld b){return b - a >= -eps;}
bool ge(ld a, ld b){return a - b > eps;}
bool le(ld a, ld b){return b - a > eps;}
bool eq(ld a, ld b){return abs(a - b) <= eps;}
bool neq(ld a, ld b){return abs(a - b) > -eps;}

struct point{
    ld x, y;
    point(): x(0), y(0){}
    point(ld x, ld y): x(x), y(y){}

    point operator+(const point & p) const{return point(x + p.x, y + p.y);}
    point operator-(const point & p) const{return point(x - p.x, y - p.y);}
    point operator*(const ld & p) const{return point(x * p, y * p);}
    point operator/(const ld & p) const{return point(x / p, y / p);}

    point operator+=(const point & p){*this = *this + p; return *this;}
    point operator-=(const point & p){*this = *this - p; return *this;}
    point operator*=(const ld & p){*this = *this * p; return *this;}
    point operator/=(const ld & p){*this = *this / p; return *this;}

    point perp() const{return point(-1 * y, x);}
    ld dot(const point & p) const{return x * p.x + y * p.y;}
    ld cross(const point & p) const{return x * p.y - y * p.x;}
    ld norm() const{return x * x + y * y;}
    ld length() const{return sqrtl(x * x + y * y);}
    point unit() const{return (*this) / length();}

    bool operator==(const point & p) const{return eq(x, p.x) && eq(y, p.y);}
    bool operator!=(const point & p) const{return !(*this == p);}
    bool operator<(const point & p) const{return le(x, p.x) || (eq(x, p.x) && le(y, p.y));}
    bool operator>(const point & p) const{return ge(x, p.x) || (eq(x, p.x) && ge(y, p.y));}

};

struct line{
    point a, v;
    line(): a(), v(){}
    line(const point & a, const point & v): a(a), v(v){}
};

int sgn(ld x){
    if(ge(x, 0)) return 1;
    if(le(x, 0)) return -1;
    return 0;
}

bool pointInLine(const point & a, const point & v, const point & p){
    return eq((p - a).cross(v), 0);
}

bool pointInSegment(const point & a, const point & b, const point & p){
    return pointInLine(a, b - a, p) && leq((a - p).dot(b - p), 0);
}

int intersectSegmentsInfo(const point & a, const point & b, const point & c, const point & d){
    point v1 = b - a, v2 = d - c;
    int t = sgn(v1.cross(c - a)), u = sgn(v1.cross(d - a));
    if(t == u){
        if(t == 0){
            if(pointInSegment(a, b, c) || pointInSegment(a, b, d) || pointInSegment(c, d, a) || pointInSegment(c, d, b)){
                return -1;
            }else{
                return 0;
            }
        }else{
            return 0;
        }
    }else{
        return sgn(v2.cross(a - c)) != sgn(v2.cross(b - c));
    }
}

point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
    ld det = v1.cross(v2);
    return a1 + v1 * ((a2 - a1).cross(v2) / det);
}

bool lexCompare(const point & a, const point & b){
    if(neq(a.x, b.x))
        return a.x < b.x;
    return a.y < b.y;
}

char segmentType(line seg, point o){
    if(eq(seg.a.x, seg.v.x))
        return 0;
    if(!lexCompare(seg.a, seg.v))
        swap(seg.a, seg.v);
    return (seg.v - seg.a).cross(o - seg.a) > 0 ? 1 : -1;
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
    for(int i = 0; i < n * 3; i++){
        if(!segmentsType[i]) continue;
        ld l = segments[i].a.x, r = segments[i].v.x;
        vector<pair<ld, int>> evs;
        for(int j = 0; j < n * 3; j++){
            if(!segmentsType[i] || i == j) continue;
            ld l1 = segments[j].a.x, r1 = segments[j].v.x;
            if(geq(l1, r) || geq(l, r1)) continue;
            ld coml = max(l, l1), comr = min(r, r1);
            auto f = intersectSegmentsInfo(segments[i].a, segments[i].v, segments[j].a, segments[j].v);
            if(f == 0){
                ld yl1 = k[j] * coml + t[j], yl = k[i] * coml + t[i];
                if(le(yl1, yl) == (segmentsType[i] == 1)){
                    int evTy = -segmentsType[i] * segmentsType[j];
                    evs.emplace_back(coml, evTy);
                    evs.emplace_back(comr, -evTy);
                }
            }else if(f == 1){
                auto u = intersectLines(segments[i].a, segments[i].v - segments[i].a, segments[j].a, segments[j].v - segments[j].a);
                ld yl = k[i] * coml + t[i], yl1 = k[j] * coml + t[j];
                int evTy = -segmentsType[i] * segmentsType[j];
                if(le(yl1, yl) == (segmentsType[i] == 1)){
                    evs.emplace_back(coml, evTy);
                    evs.emplace_back(u.x, -evTy);
                }
                yl = k[i] * comr + t[i], yl1 = k[j] * comr + t[j];
                if(le(yl1, yl) == (segmentsType[i] == 1)){
                    evs.emplace_back(u.x, evTy);
                    evs.emplace_back(comr, -evTy);
                }
            }else{
                if(segmentsType[i] != segmentsType[j] || j > i){
                    evs.emplace_back(coml, -2);
                    evs.emplace_back(comr, 2);
                }
            }
        }
        evs.emplace_back(l, 0);
        sort(all(evs));
        int j = 0, balance = 0;
        while(j < evs.size()){
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
    return ans / 2;
}

vector<pair<point, pair<point, point>>> triangulate(vector<point> & P){
    int n = P.size();
    vector<int> next(n);
    for(int i = 0; i < n - 1; i++) next[i] = i + 1;
    auto is_ear = [&](int i, int j, int k){
        if(sgn((P[j] - P[i]).cross(P[k] - P[i])) <= 0) return false;
        for(int l = next[k]; l != i; l = next[l])
            if(sgn((P[i] - P[l]).cross(P[j] - P[l])) >= 0 &&
             sgn((P[j] - P[l]).cross(P[k] - P[l])) >= 0 &&
             sgn((P[k] - P[l]).cross(P[i] - P[l])) >= 0) return false;
        return true;
    };
    vector<pair<point, pair<point, point>>> res;
    for(int i = 0; next[next[i]] != i;){
        if(is_ear(i, next[i], next[next[i]])){
            res.pb({P[i], {P[next[i]], P[next[next[i]]]}});
            next[i] = next[next[i]];
        }else i = next[i];
    }
    return res;
}

int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int n; cin>>n;
    ld res1 = 0, res2 = 0;
    vector<pair<point, pair<point, point>>> ans;
    for(int i = 0; i < n; i++){
        int m; cin>>m;
        vector<point> t(m);
        for(int j = 0; j < m; j++){
            int a, b; cin>>a>>b;
            t[j] = point(a, b);
        }
        ld e = area(t);
        res1 += abs(e);
        if(e < 0){
            reverse(all(t));
        }
        auto u = triangulate(t);
        for(int j = 0; j < u.size(); j++){
            ans.pb({u[j].fi, {u[j].se.fi, u[j].se.se}});
        }
    }
    res2 = areaUnionTriangles(ans);
    cout<<setprecision(20)<<res1<<" "<<res2;
}
