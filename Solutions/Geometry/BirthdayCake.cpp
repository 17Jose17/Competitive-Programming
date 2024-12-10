/*
  Problem: https://codeforces.com/gym/105216/problem/B
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
const ll INF = 1e12;
const int maxn = 1e5;

using ld = int;
const ld eps = 0, inf = 1e9, pi = acos(-1);
// For use with integers, just set eps=0 and everything remains the same
bool geq(ld a, ld b){return a-b >= -eps;}     //a >= b
bool leq(ld a, ld b){return b-a >= -eps;}     //a <= b
bool ge(ld a, ld b){return a-b > eps;}        //a > b
bool le(ld a, ld b){return b-a > eps;}        //a < b
bool eq(ld a, ld b){return fabs(a-b) <= eps;}  //a == b
bool neq(ld a, ld b){return fabs(a-b) > eps;}  //a != b

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

pair<vector<point>, vector<point>> convexHull(vector<point> P){
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
	sort(all(U));
	return {L, U};
}

int main(){
    ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int n, k; cin>>n>>k; k--;
    vector<point> v(n);
    for(int i = 0; i < n; i++){
        int a, b; cin>>a>>b; v[i] = point(a, b);
    }
    auto ss = convexHull(v);
    auto tap = ss.se, fad = ss.fi;
    int mat[5005][5005] = {};
    int l1 = 0, l2 = 0;
    if(fad[l1].x == fad[l1 + 1].x) l1++; if(tap[l2].x == tap[l2 + 1].x) l2++;
    for(int i = fad[0].x; i <= fad[fad.size() - 1].x; i++){
        if(fad[l1 + 1].x < i) l1++; if(tap[l2 + 1].x < i) l2++;
        for(int l = 0; l <= 5000; l++){
            if((fad[l1 + 1] - fad[l1]).cross(point(i, l) - fad[l1]) >= 0 &&
             (tap[l2 + 1] - tap[l2]).cross(point(i, l) - tap[l2]) <= 0) mat[i % k][l % k]++;
        }
    }
    int ans = 0;
    for(int i = 0; i <= 5000; i++){
        for(int l = 0; l <= 5000; l++){
            if(mat[i][l] > 1) ans += mat[i][l];
        }
    }
    cout<<ans;
}
