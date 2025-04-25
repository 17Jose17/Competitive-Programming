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

	ll cross(const point & p, const point & q) const{return (p - *this).cross(q - *this);}
	ll dot(const point & p, const point & q) const{return (p - *this).dot(q - *this);}
    
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

int sgn2(ld x) {
    return x < 0 ? -1 : 1;
}

struct ln {
    point p, pq;
    ln(const point &p, const point &q): p(p), pq(q - p) {}
    ln() {}

    bool operator/(const ln &l) {
        return abs(pq.unit().cross(l.pq.unit())) <= eps;
    }
    point operator ^ (const ln &l) {
        if (*this / l)
            return {Inf, Inf};

        return l.p + l.pq * ((p - l.p).cross(pq) / (l.pq.cross(pq)));
    }
};

struct halfplane : public ln {
    ld angle;

    halfplane() {}
    halfplane(const point &a, const point &b) {
        p = a;
        pq = b - a;
        angle = atan2(pq.y, pq.x);
    }

    bool operator < (const halfplane &b) const {
        return angle < b.angle;
    }
    bool out(const point &q) const {
        return pq.cross(q - p) < -eps;
    }
};

vector<point> intersect(vector<halfplane> b) {
    vector<point> bx = {
        {Inf, Inf},
        {-Inf, Inf},
        {-Inf, -Inf},
        {Inf, -Inf}
    };

    forn(i, 4) b.pb(halfplane(bx[i], bx[(i + 1) % 4]));

    sort(all(b));

    int n = b.size(), q = 1, h = 0;
    vector<halfplane> c(b.size() + 10);

    forn(i, n) {
        while (q < h && b[i].out(c[h] ^ c[h - 1]))
            h--;

        while (q < h && b[i].out(c[q] ^ c[q + 1]))
            q++;

        c[++h] = b[i];

        if (q < h && abs(c[h].pq.cross(c[h - 1].pq)) < eps) {
            if (c[h].pq.dot(c[h - 1].pq) <= 0)
                return {};

            h--;

            if (b[i].out(c[h].p))
                c[h] = b[i];
        }
    }

    while (q < h - 1 && c[q].out(c[h] ^ c[h - 1]))
        h--;

    while (q < h - 1 && c[h].out(c[q] ^ c[q + 1]))
        q++;

    if (h - q <= 1)
        return {};

    c[h + 1] = c[q];

    vector<point> s;

    forr(i, q, h) s.pb(c[i] ^ c[i + 1]);

    return s;
}

ld distancePointSegment(point a, point b, point p) {
    if ((b - a).dot(p - a) < 0)
        return (a - p).length();
    if ((a - b).dot(p - b) < 0)
        return (b - p).length();
    return distancePointLine(a, b - a, p);
}
 
ld distSegments(point p0, point p1, point p2, point p3){
    ld ans = INF;
    ans = min(ans, distancePointSegment(p0, p1, p2));
    ans = min(ans, distancePointSegment(p0, p1, p3));
    ans = min(ans, distancePointSegment(p2, p3, p0));
    ans = min(ans, distancePointSegment(p2, p3, p1));
    return ans;
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

vector<point> tangentsPointPolygon(const vector<point> & P, const vector<vector<point>> & Ps, const point & p){
	int n = P.size(), m = Ps[0].size(), k = Ps[1].size();
	
	int lk = m; if(Ps[0][m - 1] == Ps[1][0]) lk--; 

	auto tang = [&](int l, int r, ld w, int kl) -> int {
		int res = min(l, r);
	        while(l <= r){
			int m = (l + r) / 2;
			ld a = (P[(m + kl) % n] - p).cross(P[(m + 1 + kl) % n] - p) * w, b = (P[(m + kl) % n] - p).cross(P[(m - 1 + n + kl) % n] - p) * w;
			if(geq(a, 0) && geq(b, 0)) return m;
			if(geq(a, 0)) r = m - 1, res = m;
			else l = m + 1;
	        }
	        return res;
    	};

	auto bs = [&](int l, int r, const vector<point> & A, ld w) -> int {
	        int res = l;
	        ld w1 = p.x * w;
	        while(l <= r){
			int m = (l + r) / 2;
			if(ge(A[m].x * w, w1)) r = m - 1;
			else res = m, l = m + 1;
		}
	        return res;
	};
    
	point left = p, rigth = p;

	int t1 = bs(0, m - 1, Ps[0], 1), t2 = bs(0, k - 1, Ps[1], -1);
	
	auto u1 = tang(0, t1, -1, 0), u2 = tang(0, t2, -1, lk);
	auto v1 = tang(t1, m - 1, 1, 0), v2 = tang(t2, k - 1, 1, lk);
	
	if(leq((P[u1] - p).cross(P[(u1 - 1 + n) % n] - p), 0) && leq((P[u1] - p).cross(P[(u1 + 1) % n] - p), 0)) left = P[u1];
	if(leq((P[(lk + u2) % n] - p).cross(P[(lk + u2 - 1 + n) % n] - p), 0) && leq((P[(lk + u2) % n] - p).cross(P[(lk + u2 + 1) % n] - p), 0)) left = P[(lk + u2) % n];
	
	if(geq((P[v1] - p).cross(P[(v1 - 1 + n) % n] - p), 0) && geq((P[v1] - p).cross(P[(v1 + 1) % n] - p), 0)) rigth = P[v1];
	if(geq((P[(lk + v2) % n] - p).cross(P[(lk - 1 + n + v2) % n] - p), 0) && geq((P[(lk + v2) % n] - p).cross(P[(lk + 1 + v2) % n] - p), 0)) rigth = P[(lk + v2) % n];
    
	return {left, rigth};
}

vector<vector<point>> trian(vector<point> & P){
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
    vector<vector<point>> res;
    for(int i = 0; next[next[i]] != i;){
        if(is_ear(i, next[i], next[next[i]])){
            res.pb({P[i], P[next[i]], P[next[next[i]]]});
            next[i] = next[next[i]];
        }else i = next[i];
    }
    return res;
}

ld segmentCircular(point c, point a, point b){
    ld r = (c - a).length();
    point v1 = a - c, v2 = b - c;
    ld tet = acosl(v1.dot(v2) / v1.length() / v2.length());
    point rt = v1.rotate(tet);
    if(!(abs(rt.x - v2.x) <= .5 && abs(rt.y - v2.y) <= .5)) tet = pi * 2 - tet;
    return (a - c).length() * tet;
}

struct line{
    point a, v;
    line(): a(), v(){}
    line(const point & a, const point & v): a(a), v(v){}
};

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

ld areaUnionTriangles(vector<vector<point>> & P){
    int n = P.size();
    vector<line> segments(n * 3); vector<char> segmentsType(n * 3);
    for(int i = 0; i < n; i++){
        point a = P[i][0], b = P[i][1], c = P[i][2];
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

int sgn(ll x) { return leq(x, 0) ? eq(x, 0) ? 0 : -1 : 1; }

using pt = point;

struct Edge {
    pt l, r;
};

bool edge_cmp(Edge* edge1, Edge* edge2)
{
    const pt a = edge1->l, b = edge1->r;
    const pt c = edge2->l, d = edge2->r;
    int val = sgn(a.cross(b, c)) + sgn(a.cross(b, d));
    if (val != 0)
        return val > 0;
    val = sgn(c.cross(d, a)) + sgn(c.cross(d, b));
    return val < 0;
}

enum EventType { DEL = 2, ADD = 3, GET = 1, VERT = 0 };

struct Event {
    EventType type;
    int pos;
    bool operator<(const Event& event) const { return type < event.type; }
};

vector<Edge*> sweepline(vector<Edge*> planar, vector<pt> queries)
{
    using pt_type = decltype(pt::x);

    auto s =
        set<pt_type, std::function<bool(const pt_type&, const pt_type&)>>(le);
    for (pt p : queries)
        s.insert(p.x);
    for (Edge* e : planar) {
        s.insert(e->l.x);
        s.insert(e->r.x);
    }

    int cid = 0;
    auto id =
        map<pt_type, int, std::function<bool(const pt_type&, const pt_type&)>>(
            le);
    for (auto x : s)
        id[x] = cid++;

    // create events
    auto t = set<Edge*, decltype(*edge_cmp)>(edge_cmp);
    auto vert_cmp = [](const pair<pt_type, int>& l,
                       const pair<pt_type, int>& r) {
        if (!eq(l.first, r.first))
            return le(l.first, r.first);
        return l.second < r.second;
    };
    auto vert = set<pair<pt_type, int>, decltype(vert_cmp)>(vert_cmp);
    vector<vector<Event>> events(cid);
    for (int i = 0; i < (int)queries.size(); i++) {
        int x = id[queries[i].x];
        events[x].push_back(Event{GET, i});
    }
    for (int i = 0; i < (int)planar.size(); i++) {
        int lx = id[planar[i]->l.x], rx = id[planar[i]->r.x];
        if (lx > rx) {
            swap(lx, rx);
            swap(planar[i]->l, planar[i]->r);
        }
        if (lx == rx) {
            events[lx].push_back(Event{VERT, i});
        } else {
            events[lx].push_back(Event{ADD, i});
            events[rx].push_back(Event{DEL, i});
        }
    }

    vector<Edge*> ans(queries.size(), nullptr);
    for (int x = 0; x < cid; x++) {
        sort(events[x].begin(), events[x].end());
        vert.clear();
        for (Event event : events[x]) {
            if (event.type == DEL) {
                t.erase(planar[event.pos]);
            }
            if (event.type == VERT) {
                vert.insert(make_pair(
                    min(planar[event.pos]->l.y, planar[event.pos]->r.y),
                    event.pos));
            }
            if (event.type == ADD) {
                t.insert(planar[event.pos]);
            }
            if (event.type == GET) {
                auto jt = vert.upper_bound(
                    make_pair(queries[event.pos].y, planar.size()));
                if (jt != vert.begin()) {
                    --jt;
                    int i = jt->second;
                    if (geq(max(planar[i]->l.y, planar[i]->r.y),
                           queries[event.pos].y)) {
                        ans[event.pos] = planar[i];
                        continue;
                    }
                }
                Edge* e = new Edge;
                e->l = e->r = queries[event.pos];
                auto it = t.upper_bound(e);
                if (it != t.begin())
                    ans[event.pos] = *(--it);
                delete e;
            }
        }

        for (Event event : events[x]) {
            if (event.type != GET)
                continue;
            if (ans[event.pos] != nullptr &&
                eq(ans[event.pos]->l.x, ans[event.pos]->r.x))
                continue;

            Edge* e = new Edge;
            e->l = e->r = queries[event.pos];
            auto it = t.upper_bound(e);
            delete e;
            if (it == t.begin())
                e = nullptr;
            else
                e = *(--it);
            if (ans[event.pos] == nullptr) {
                ans[event.pos] = e;
                continue;
            }
            if (e == nullptr)
                continue;
            if (e == ans[event.pos])
                continue;
            if (id[ans[event.pos]->r.x] == x) {
                if (id[e->l.x] == x) {
                    if (ge(e->l.y, ans[event.pos]->r.y))
                        ans[event.pos] = e;
                }
            } else {
                ans[event.pos] = e;
            }
        }
    }
    return ans;
}

struct DCEL {
    struct Edge {
        pt origin;
        Edge* nxt = nullptr;
        Edge* twin = nullptr;
        int face;
    };
    vector<Edge*> body;
};

vector<pair<int, int>> point_location(DCEL planar, vector<pt> queries)
{
    vector<pair<int, int>> ans(queries.size());
    vector<Edge*> planar2;
    map<intptr_t, int> pos;
    map<intptr_t, int> added_on;
    int n = planar.body.size();
    for (int i = 0; i < n; i++) {
        if (planar.body[i]->face > planar.body[i]->twin->face)
            continue;
        Edge* e = new Edge;
        e->l = planar.body[i]->origin;
        e->r = planar.body[i]->twin->origin;
        added_on[(intptr_t)e] = i;
        pos[(intptr_t)e] =
            le(planar.body[i]->origin.x, planar.body[i]->twin->origin.x)
                ? planar.body[i]->face
                : planar.body[i]->twin->face;
        planar2.push_back(e);
    }
    auto res = sweepline(planar2, queries);
    for (int i = 0; i < (int)queries.size(); i++) {
        if (res[i] == nullptr) {
            ans[i] = make_pair(1, -1);
            continue;
        }
        pt p = queries[i];
        pt l = res[i]->l, r = res[i]->r;
        if (eq(p.cross(l, r), 0) && leq(p.dot(l, r), 0)) {
            ans[i] = make_pair(0, added_on[(intptr_t)res[i]]);
            continue;
        }
        ans[i] = make_pair(1, pos[(intptr_t)res[i]]);
    }
    for (auto e : planar2)
        delete e;
    return ans;
}

vector<DCEL::Edge*> edges(k), twins(k);
for(int j = 0; j < k; j++){
    edges[j] = new DCEL::Edge{v[j]};
    twins[j] = new DCEL::Edge{v[(j + 1) % k]};
    edges[j]->twin = twins[j];
    twins[j]->twin = edges[j];
    edges[j]->face = i;
    twins[j]->face = -1;
}

for(int j = 0; j < k; j++){
    edges[j]->nxt = edges[(j + 1) % k];
    twins[j]->nxt = twins[(j - 1 + k) % k];
}

for(int j = 0; j < k; j++){
    dcel.body.push_back(edges[j]);
    dcel.body.push_back(twins[j]);
}

auto res = point_location(dcel, ps);
