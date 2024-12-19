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
