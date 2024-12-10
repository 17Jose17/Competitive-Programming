/*
  Problem: https://codeforces.com/gym/104614/problem/H
  Explication:
*/

#include <bits/stdc++.h>
using namespace std;

// Holi c:

#define ll long long int
//#define ld long double
#define ii __int128
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()

const int Inf = 1e9;
const ll mod = 998244353;
const ll INF = 1e18;
const int maxn = 3e5 + 5;

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

bool pointInLine(const point & a, const point & v, const point & p){
	//line a+tv, point p
	return eq((p - a).cross(v), 0);
}

bool pointInSegment(const point & a, const point & b, const point & p){
	//segment ab, point p
	return pointInLine(a, b - a, p) && leq((a - p).dot(b - p), 0);
}


point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
	//lines a1+tv1, a2+tv2
	//assuming that they intersect
	ld det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
}

int intersectLinesInfo(const point & a1, const point & v1, const point & a2, const point & v2){
	//lines a1+tv1 and a2+tv2
	ld det = v1.cross(v2);
	if(eq(det, 0)){
		if(eq((a2 - a1).cross(v1), 0)){
			return -1; //infinity points
		}else{
			return 0; //no points
		}
	}else{
		return 1; //single point
	}
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

ld f1(point c, ld r, point s0, point s1){
    point bi = s0 - s1; bi = bi.perp();
    auto at = intersectLines(c, bi, s0, s1 - s0);
    ld res = INF;
    if(at.x + eps >= min(s0.x, s1.x) && at.x - eps <= max(s0.x, s1.x)) res = min(res, (at - c).length());
    res = min(res, min((s0 - c).length(), (s1 - c).length()));
    if(res <= r + eps) return 0;
    return res - r;
}

ld f(point c, ld r0, point vt, point s0, point s1, ld in, ld fi){
    ld l = 0, r = Inf, ans = INF;
    for(int i = 0; i < 300; i++){
        ld m1 = l + (r - l) / 3, m2 = r - (r - l) / 3;
        ld g1 = f1(c + vt * m1, r0, s0, s1), g2 = f1(c + vt * m2, r0, s0, s1);
        if(!g1){
            ans = 0;
            r = m1; continue;
        }
        if(!g2){
            ans = 0;
            r = m2; continue;
        }
        if(g1 >= g2) l = m1;
        else r = m2;
    }
    point se = c + vt * l;
    point bi = s0 - s1; bi = bi.perp();
    auto at = intersectLines(se, bi, s0, s1 - s0);
    if(!(at.x + eps >= min(s0.x, s1.x) && at.x - eps <= max(s0.x, s1.x))){
        if((s0 - se).length() <= (s1 - se).length()) at = s0; else at = s1;
    }
    if(at.x - eps > fi || at.x + eps < in) return INF;
    if(!ans) return l;
    return INF;
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	int n; cin>>n; vector<point> v(n + 1);
	for(int i = 0; i <= n; i++){
	    int a, b; cin>>a>>b; v[i] = point(a, b);
	}
	int cx, sx, sy, dx, dy, v1; ld r; cin>>cx>>sx>>sy>>r>>dx>>dy>>v1;
	point ini(sx, sy); point vt(dx, dy); vt = vt.unit() * v1; point cam;
	for(int i = 0; i < n; i++){
	    if(cx >= v[i].x && cx <= v[i + 1].x){
	        cam = intersectLines(v[i], v[i + 1] - v[i], point(cx, 0), point(cx, Inf) - point(cx, 0));
	        break;
	    }
	}
	int in = v[0].x, fi = v[n].x;
	cam.x = cx;
	vector<point> side1, side2, izq, der;
	int j1 = 0; //izq.pb(v[0]);
	while(v[j1].x < cam.x) izq.pb(v[j1++]);
	if(v[j1] == cam) j1++;
	while(j1 <= n) der.pb(v[j1++]);
	int m = izq.size(), k = der.size();
	//Reconstuir el poligono en dos partes segun quienes p0.x < cam.x || p0.x > cam.x
	//Introducir todos aquellos puntos (puntos limites de la poli-linea e intersecciones de la linea de vision con los segmentos) que se puedan ver desde cam
	//side1
	for(int i = 0; i < k; i++){
	    side1.pb(der[i]);
	    for(int j = 0; j < k - 1; j++){
	        if(intersectLineSegmentInfo(cam, der[i] - cam, der[j], der[j + 1]) == 1){
	            auto it = intersectLines(cam, der[i] - cam, der[j], der[j + 1] - der[j]);
	            if(it != der[i] && it != der[i + 1]) side1.pb(it);
	        }
	    }
	}
	sort(all(side1));
	//side2
	for(int i = 0; i < m; i++){
	    side2.pb(izq[i]);
	    for(int j = 0; j < m - 1; j++){
	        if(intersectLineSegmentInfo(cam, izq[i] - cam, izq[j], izq[j + 1]) == 1){
	            auto it = intersectLines(cam, izq[i] - cam, izq[j], izq[j + 1] - izq[j]);
	            if(it != izq[j] && it != izq[j + 1]) side2.pb(it);
	        }
	    }
	} 
	sort(all(side2)); reverse(all(side2)); reverse(all(izq));
	//Limpiar los puntos, dado la anterior linea de vision mas proxima, verificar si se puede ver el punto
	vector<point> polyd, polyi;
	int j = 0, a = 1;
	for(int i = 0; i < side1.size(); i++){
	    if(j + a < k){
	        if(side1[i] == der[j + a]){
	            if((der[j] - cam).cross(side1[i] - cam) > eps) j += a, a = 1;
	            else a++;
	        } 
	    } 
	    if((der[j] - cam).cross(side1[i] - cam) >= -eps) polyd.pb(side1[i]);
	}
	j = 0; a = 1;
	for(int i = 0; i < side2.size(); i++){
	    if(j + a < m){
	        if(side2[i] == izq[j + a]){
	            if((izq[j] - cam).cross(side2[i] - cam) < -eps) j += a, a = 1;
	            else a++;
	        } 
	    } 
	    if((izq[j] - cam).cross(side2[i] - cam) <= eps) polyi.pb(side2[i]);
	}
	//Añadimos un punto final para representar el infinito
	if(polyd.size()) polyd.pb(cam + (polyd[polyd.size() - 1] - cam) * Inf);
	if(polyi.size()) polyi.pb(cam + (polyi[polyi.size() - 1] - cam) * Inf);
	//Verificar el tiempo que tarda la nube en cbocar con algun segmento del poliogno
	ld ans = INF;
	for(int i = 0; i < polyd.size(); i++){
	    point p0, p1 = polyd[i];
	    if(i) p0 = polyd[i - 1]; else p0 = cam;
	    if(p0 == p1) continue;
	    if(abs((p0 - p1).cross(vt)) <= eps){
	        if(abs((p0 - ini).cross(p1 - ini)) <= eps) ans = min(ans, min((abs(p0.y - ini.y) - r) / vt.y, (abs(p1.y - ini.y) - r) / vt.y));
	    }else{
	        ans = min(ans, f(ini, r, vt, p0, p1, in, fi));
	    }
	}
	for(int i = 0; i < polyi.size(); i++){
	    point p0, p1 = polyi[i];
	    if(i) p0 = polyi[i - 1]; else p0 = cam;
	    if(p0 == p1) continue;
	    if(abs((p0 - p1).cross(vt)) <= eps){
	        if(abs((p0 - ini).cross(p1 - ini)) <= eps) ans = min(ans, min((abs(p0.y - ini.y) - r) / vt.y, (abs(p1.y - ini.y) - r) / vt.y));
	    }else{
	        ans = min(ans, f(ini, r, vt, p0, p1, in, fi));
	    }
	}
	if(ans == INF) cout<<-1;
	else cout<<setprecision(25)<<ans;
}
