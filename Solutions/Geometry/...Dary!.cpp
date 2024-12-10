/*
  Problem: https://codeforces.com/contest/696/problem/F
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

const int Inf = 1e3;
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

int pointInConvexPolygon(const vector<point> & P, const point & p, int right){
	if(p < P[0] || P[right] < p) return 0;
	int orientation = sgn((P[right] - P[0]).cross(p - P[0]));
	if(orientation == 0){
		if(p == P[0] || p == P[right]) return -1;
		return (right == 1 || right + 1 == P.size()) ? -1 : 1;
	}else if(orientation < 0){
		auto r = lower_bound(P.begin() + 1, P.begin() + right, p);
		int det = sgn((p - r[-1]).cross(r[0] - r[-1])) - 1;
		if(det == -2) det = 1;
		return det;
	}else{
		auto l = upper_bound(P.rbegin(), P.rend() - right - 1, p);
		int det = sgn((p - l[0]).cross((l == P.rbegin() ? P[0] : l[-1]) - l[0])) - 1;
		if(det == -2) det = 1;
		return det;
	}
}

point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
	//lines a1+tv1, a2+tv2
	//assuming that they intersect
	ld det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
}

pair<bool, point> f(point a, point b, point c, point d, vector<point> & P, int right, ld m){
    point v1 = (b - a).unit() * m, v2 = (d - c).unit() * m;
    point a1 = a + v1.perp(), c1 = c + v2.perp();
    auto it = intersectLines(a1, b - a, c1, d - c);
    if(pointInConvexPolygon(P, it, right) != 0) return {true, it};
    return {false, it};
}

int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int n; cin>>n;
    vector<point> v(n);
    for(int i = 0; i < n; i++){
        int a, b; cin>>a>>b;
        v[i] = point(a, b);
    }
    if(n == 3){
        cout<<0<<'\n'<<v[0].x<<" "<<v[0].y<<'\n'<<v[1].x<<" "<<v[1].y;
        return 0;
    }else if(n == 4){
        cout<<0<<'\n'<<v[0].x<<" "<<v[0].y<<'\n'<<v[2].x<<" "<<v[2].y;
        return 0;
    }
    rotate(v.begin(), min_element(v.begin(), v.end()), v.end());
    int right = max_element(v.begin(), v.end()) - v.begin();
    ld ans = INF; 
    point res1, res2;
    for(int i = 0; i < n; i++){
        int j = i, k = (i - 1 + n) % n, k1 = (i - 2 + n) % n, s = (j + 1) % n;
        while(1){
            int j1 = (j + 1) % n, s1 = (s + 1) % n;
            if(eq((v[i] - v[k]).cross(v[j1] - v[s1]), 0)) break;
            auto it = intersectLines(v[i], v[i] - v[k], v[j1], v[s1] - v[j1]);
            if(ge((v[j1] - v[i]).cross(it - v[i]), 0)) break;
            j = j1; s = s1;
        }
        int s2 = (s + 1) % n;
        // lines v[i], v[k] with v[j], v[s] y v[k], v[k1] with v[s], v[s2]
        ld l = 0, r = (v[i] - v[j]).length() / 2;
        point sup;
        for(int w = 0; w < 100; w++){
            ld m = (l + r) / 2;
            auto z = f(v[k], v[i], v[j], v[s], v, right, m);
            if(z.fi){
                r = m;
                sup = z.se;
            } else l = m;
        }
        point sep;
        ld l1 = 0, r1 = (v[k1] - v[s2]).length() / 2;
        for(int w = 0; w < 100; w++){
            ld m = (l1 + r1) / 2;
            auto z = f(v[s], v[s2], v[k1], v[k], v, right, m);
            if(z.fi){
                r1 = m;
                sep = z.se;
            }else l1 = m;
        }
        if(max(l, l1) < ans){
            ans = max(l, l1);
            res1 = sup; res2 = sep;
        }
    }
    cout<<setprecision(20)<<ans<<'\n'<<res1.x<<" "<<res1.y<<'\n'<<res2.x<<" "<<res2.y;
}
