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

	bool left(const point & p, const point & q){ 
        return (q - p).cross(*this - p) > eps;
    }
    
};

struct line{
    int a, b, c;
    line(): a(0), b(0), c(0){}
	line(int a, int b, int c): a(a), b(b), c(c){}
	line(point p0, point p1): line(p0.y - p1.y, p1.x - p0.x, (p0.y - p1.y) * p0.x * -1 - (p1.x - p0.x) * p0.y){}
	
	void norm(){
	    int g = __gcd(abs(a), __gcd(abs(b), abs(c)));
	    if(g != 0){
	        a /= g; b /= g; c /= g;
	    }
	    if(a < 0 || (a == 0 && b < 0)){
	        a *= -1; b *= -1; c *= -1;
	    }
	}
	
	bool operator==(const line & l) const{return a == l.a && b == l.b && c == l.c;}
	bool operator<(const line & l) const{
	    if(a != l.a) return a < l.a;
	    if(b != l.b) return b < l.b;
	    return c < l.c;
	}
};

ld distancePointSegment(point a, point b, point p){
    ld res = min((p - a).length(), (p - b).length());
    if((b - a).cross(p - a) && (b - a).dot(p - a) >= 0 && (a - b).dot(p - b) >= 0) res = min(res, distancePointLine(a, b - a, p));
    return res;
}

ld distancePointPolygon(vector<point> P, point p){
        ld res = INF; int n = P.size();
        for(int i = 0; i < n; i++){
                int j = i + 1; if(j == n) j = 0;
                res = min(res, distancePointSegment(P[i], P[j], p));
        }
        return res;
}

// O(n + m)
ld distanceConvexPolygons(vector<point> P1, vector<point> P2){
        int nb = P2.size();
        for(int i = 0; i < nb; i++) P2[i] *= -1;
        sort(all(P1)); P1 = convexHull(P1); sort(all(P2)); P2 = convexHull(P2);
        auto P = minkowskiSum(P1, P2);
        ld res = INF; int n = P.size();
        for(int i = 0; i < n; i++){
                int j = i + 1; if(j == n) j = 0;
                res = min(res, distancePointSegment(P[i], P[j], point(0, 0)));
        }
        return res;
}
