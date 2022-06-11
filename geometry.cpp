/* never forget collinearity corner cases */

#include <bits/stdc++.h>
#define MOD 1000000007
#define pb push_back
#define F first
#define S second

using namespace std;
typedef long long ll;
typedef long double ld;
typedef pair<double, double> pdd;
typedef pair<int, int> pii;
const double PI = acos(-1.0);
const double eps = 1e-9;

const ll MX = 2e5 + 9;
const ll inf = 1e18;

template <typename T> int sgn(T x) {
    return (T(0) < x) - (x < T(0));
}

struct P {
    int x, y;
    
    void read() {
        cin >> x >> y;
    }
    
    void Out() {
        cout << "(x, y) " << x << " " << y << "\n";
    }
    
    P (): x(0), y(0) {}
    
    P (int a, int b) : x(a), y(b) {}
    
    P operator - (const P& b) const {
        return P {x - b.x, y - b.y};
    }
    
    void operator -= (const P& b) {
        x -= b.x, y -= b.y;
    }
    
    P operator + (const P& b) const {
        return P {x + b.x, y + b.y};
    }
    
    void operator += (const P& b) {
        x += b.x, y += b.y;
    }
    
    ll operator * (const P& b) {
        return ll (x * b.y) - ll (y * b.x);
    }
    
    bool operator != (const P &b) {
        return x != b.x || y != b.y;
    }
    
    bool operator == (const P &b) {
        return x == b.x && y == b.y;
    }
    
    P operator * (const double& b) {
        return P (x * b, y * b);
    }
    
    P operator / (const double& b) {
        return P (x / b, y / b);
    }
    
    double triangle (const P &b, const P& c) {
        return (b - *this) * (c - *this);
    }
    
    friend bool operator < (const P& a, const P& b) {
        return make_pair(a.x, a.y) < make_pair(b.x, b.y);
    }
    
    int sq (P p) {
        return p.x * p.x + p.y * p.y;
    }
    
    double absP (P p) {
        return sqrt(sq(p));
    }
};

typedef pair<P, P> ppp;

int sq (P p) {
    return p.x * p.x + p.y * p.y;
}

double abs (P p) {
    return sqrt (sq(p));
}

P perp (P p) {
    return P (-p.y, p.x);
}

P translate (P v, P p) {
    return p + v;
}

P scale (P c, double factor, P p) {
    return c + (p - c) * factor;
}

P rot (P p, double a) {
    return P (p.x * cos (a) - p.y * sin (a),
        p.x * sin(a) + p.y * cos (a));
}

int dot (P v, P w) {
    return v.x * w.x + v.y * w.y;
}

struct Line {
    P v; int c; // c is int but you might need to change it! 
    
    Line (P _v, int _c) : v(_v), c(_c){} //from direction of v and offset c
    Line (int a, int b, int _c) : v(P(b, -a)), c(_c) {} //from equation ax + by = c
    Line (P p, P q) : v(P(q - p)), c(v * p) {} // from points P and Q
    
    int side (P p) {
        return v * p - c;
    }
    
    double dist (P p) {
        return abs (side(p)) / abs(v);
    }
    
    double sqDist (P p) {
        return side (p) * side (p) / (double) sq(v);
    }
    
    Line perpThrough (P p) {
        return Line (p, p + perp(v));
    }
    
    bool cmpProj (P p, P q) { 
        //to sort points in order they appear on the line
        //also compares two points by their orthogonal projection
        
        return dot (v, p) < dot (v, q);
    }
    
    //translate a line by vector 
    Line translate (P t) {
        return Line (v, c + v * t);
    }
    
    // translate with fixed distance to the left
    Line shiftLeft (double dist) {
        return Line (v, c + dist * abs(v));
    }
    
    P proj (P p) {//projection
        return p - perp(v) * (side (p) / sq (v));
    }
    
    P refl (P p) {//reflection
        return p - perp(v) * 2 * (side(p) / sq(v));
    }
};

typedef pair <P, P> Segment;

ll Sign (P& p, Segment seg) {
    ll x = seg.F.triangle (p, seg.S);
    return x == 0ll ? 0ll : (x / abs (x));
}

double angle (P v, P w) {
    double cosTheta = dot (v, w) / abs (v) / abs (w);
    return acos (max (-1.0, min(1.0, cosTheta)));
}

bool SegSegIntersection (Segment seg1, Segment seg2) {
    P A = seg1.F, B = seg1.S, C = seg2.F, D = seg2.S;
    if (A.triangle (C, D) == 0 && B.triangle (C, D) == 0) {
        for (int rep = 0; rep < 2; rep++) {
            if (max (A.x, B.x) < min (C.x, D.x) || max (A.y, B.y) < min (C.y, D.y))
                return false;
                
            swap (A, C), swap (B, D);
        }
        
        return true;
    }
    
    for (int rep = 0; rep < 2; rep++) {
        ll s1 = Sign (A, Segment {C, D}), s2 = Sign (B, Segment {C, D});
        if (s1 == s2 && s1 != 0)
            return false;
            
        swap (A, C), swap (B, D);
    }
    
    return true;
}

Line bisector (Line l1, Line l2, bool interior) {
    assert (l1.v * l2.v != 0);
    double sign = interior ? 1 : -1;
    
    return Line (l2.v / abs(l2.v) + l1.v/ abs(l1.v) * sign, 
        l2.c/abs(l2.v) + l1.c/(abs(l1.v) * sign));
}

bool PointInSegment (P& p, Segment seg) {
    if (p.triangle (seg.F, seg.S) != 0)
        return false;
        
    // p is somewhere inside the bounding box of the segment
    return min (seg.F.x, seg.S.x) <= p.x && p.x <= max (seg.F.x, seg.S.x)
        && min (seg.F.y, seg.S.y) <= p.y && p.y <= max (seg.F.y, seg.S.y);
}

int quadRoots (double a, double b, double c, pdd &out) {
    assert (a != 0);
    double disc = b * b - 4 * a * c;
    
    if (disc < 0)
        return 0;
        
    double sum = (b >= 0.0) ? (-b-sqrt(disc)) : (-b+sqrt(disc));
    out = {sum / (2 * a), sum == 0 ? 0 : (2 * c) / sum};
    return 1 + (disc > 0);
    
}

double orient (P a, P b, P c) {
    return a.triangle (b, c);
}

bool inDisk (P a, P b, P p) {
    return dot (a - p, b - p) <= 0;
}

bool onSegment (P a, P b, P p) {
    //this the SECOND WAY to check if point on segment
    return a.triangle (b, p) == 0 && inDisk (a, b, p);
}

bool properInter (P a, P b, P c, P d, P& out) {
    //[A, B], [C, D]
    double oa = orient (c, d, a), 
        ob = orient (c, d, b), 
        oc = orient (a, b, c),
        od = orient (a, b, d);
        
    if (oa * ob < 0 && oc * od < 0) {
        out = (a * ob - b * oa) / (ob - oa);
        return true;
    }
    
    return false;
}

set <P> inters (P a, P b, P c, P d) {
    P out;
    
    if (properInter (a, b, c, d, out)) {
        set <P> s;
        s.insert (out);
        return s;
    }
    
    set <P> s;
    if (onSegment(c, d, a)) s.insert (a);
    if (onSegment(c, d, b)) s.insert (b);
    if (onSegment(a, b, c)) s.insert (c);
    if (onSegment(a, b, d)) s.insert (d);
    
    return s;
}

double segPoint (P a, P b, P p) {
    /* either the projection is inside the segment or 
     * the distance is to one of the endpoints */
    if (a != b) {
        Line l(a, b);
        if (l.cmpProj (a, p) && l.cmpProj (p, b))
            return l.dist (p);
    }
    
    //distance to either of points
    return min (abs(p - a), abs(p - b));
}

double SegSegDistance (P a, P b, P c, P d) {
    P dummy;
    if (properInter (a, b, c, d, dummy))
        return 0; // answer is 0 if segments intersect properly
        
    return min (min(segPoint(a, b, c), segPoint(a, b, d)),
        min(segPoint(c, d, a), segPoint(c, d, b)));
    /*otherwise the shortest distance is to one of the endpoints*/
}

class PointInPolygon {
    //Technique #1: Imaginary Vertical Ray
    bool chk1 (P &p, vector <P>& Polygon) {
        //Try 3e9, be careful with hacking testcases
        //take care of overflow
        //another way is to imagine the upper point (q) tilted by epsilon to the right
        //and care about p.triangle(A,B) < 0 && (A.x <= p.x) && (B.x <= p.x)
        P q (p.x + 1, 1e9 + 1);
        int cnt = 0;
    
        for (int i = 0; i < (int) Polygon.size(); i++) {
            P A = Polygon[i], B = Polygon[((i + 1) < (int) Polygon.size()) ? (i + 1) : 0];
            if (PointInSegment (p, Segment {A, B}))
                return 0;
            
            cnt += SegSegIntersection (Segment {A, B}, Segment {p, q});
        }
        
        return (cnt & 1);
    }
    
    bool above (P a, P p) {//true if P at least as high as A (blue part)
        return p.y >= a.y;
    }
    
    bool crossesRay (P a, P p, P q) {//check if [PQ] crosses ray from A
        return (above (a, q) - above (a, p)) * orient (a, p, q) > 0;
    }
    
    //Technique #2: Imaginary Horizontical Ray
    bool chk2 (vector <P> p, P a, bool strict = true) {
        //if strict, return false when A is on the boundary
        int cnt = 0;
        for (int i = 0, n = p.size(); i < n; i++) {
            if (onSegment (p[i], p[(i + 1) % n], a))
                return !strict;
            
            cnt += crossesRay (a, p[i], p[(i + 1) % n]);
        }
        
        return (cnt & 1);
    }
};

/* winding number with floats */
double angleTravelled (P a, P p, P q) {
    double ampli = angle (p - a, q - a);
    if (orient (a, p, q) > 0)
        return ampli;
    
    return -ampli;
}

int windingNumber (vector <P> p, P a) {
    double ampli = 0;
    for (int i = 0, n = p.size(); i < n; i++)
        ampli += angleTravelled (a, p[i], p[(i + 1) % n]);
        
    return round (ampli / (2 * PI));
}
/* winding number with floats */

bool up (P p) {
    assert (p.x != 0 || p.y != 0); //argument of 0 0 is undefined
    
    return p.y > 0 || (p.y == 0 && p.x < 0);
}

/* winding number INTEGER COORDINATES */
struct Angle {
    P d; int t = 0;
    Angle () {}
    Angle (P pp, int tt) {
        d = pp;
        t = tt;
    }
    
    Angle t180(){
        return Angle(d * (-1), t + up (d));
    }
    
    Angle t360() {
        return Angle(d, t + 1);
    }
    
    friend bool operator < (Angle a, Angle b) {
        return make_tuple (a.t, up(a.d), 0) < make_tuple (b.t, up(b.d), a.d * b.d);
    }
};

Angle moveTo (Angle a, P newD) {
    assert (!onSegment(a.d, newD, P(0, 0)));
    
    Angle b(newD, a.t);
    if (a.t180() < b)
        b.t--;
    if (b.t180() < a)
        b.t++;
    
    return b;
}

int windingNumberInteger (vector <P> p, P a) {
    Angle A (p.back() - a, 0);
    for (P d: p)
        A = moveTo (A, d - a);
        
    return A.t;
}

/* winding number INTEGER COORDINATES */

#include <bits/stdc++.h>
using namespace std;

P circumCenter (P a, P  b, P c) {
    // it is a circle that passes throught triangle (a, b, c);
    b -= a, c -= a;
    assert (b * c != 0);
    
    return a + perp(b * sq(c) - c * sq(b)) / (b * c) / 2;
    // to find r find the distanct from center to a or b or c
}

int circleLine (P o, double r, Line l, ppp &out) {
    double h2 = r * r - l.sqDist(o);
    if (h2 >= 0) {//the line touches the circle
        P p = l.proj(o);//point P
        P h = l.v*sqrt(h2)/abs(l.v);//vector parallel to l, of length h
        
        out = {p - h, p + h};
    }
    
    return 1 + sgn (h2);
}

int circleCircle (P o1, double r1, P o2, double r2, ppp &out) {
    P d = o2 - o1;
    double d2 = sq(d);
    
    if (d2 == 0) {
        assert (r1 != r2);
        return 0;
    }
    
    double pd = (d2 + r1*r1 - r2*r2)/2;
    double h2 = r1*r1 - pd*pd/d2;//h^2
    if (h2 >= 0) {
        P p = o1 + d*pd/d2, h = perp(d)*sqrt(h2/d2);
        out = {p-h, p+h};
    }
    
    return 1 + sgn(h2);
}

int tangents (P o1, double r1, P o2, double r2, bool inner, vector <ppp> &out) {
    if (inner) r2 = -r1;
    P d = o2 - o1;
    double dr = r1 - r2, d2 = sq(d), h2 = d2 - dr * dr;
    if (d2 == 0 || h2 < 0) {
        assert (h2 != 0);
        return 0;
    }
    
    for (double sign: {-1, 1}) {
        P v = (d * dr + perp(d) * sqrt(h2) * sign) / d2;
        out.pb ({o1 + v*r1, o2 + v*r2});
    }
    
    return 1 + (h2 > 0); 
    /* -conveniently, the same code can be used to find the tangent to a 
    circle passing through a point by setting t2 to 0 (in which case the value 
    of inner does not matter)
    - if there are 2 tangents, it fills out with two pairs of points: the pairs
    of tangency points on each circle (P1, P2), for each of the tangents
    - if there is 1 tangent, the circles are tangent to each other at some
    point P, out just contains P 4 times, and the tangent line can be found as 
    line(o1,p).perpThrough(p)
    - if there are 0 tangents, it does nothing.
    - if the circles are identical, it aborts.
    */
}

vector <P> convexHull (vector <P> points) {
    sort (points.begin(), points.end());
    vector <P> hull;
    
    for (int rep = 0; rep < 2; rep++) {
        const int S = hull.size();
        
        for (P C: points) {
            while ((int) hull.size() - S >= 2) {
                P A = hull.end() [-2];
                P B = hull.end() [-1];
                
                if (A.triangle (B, C) <= 0ll) {
                    // take care of <= or < for consecutive points on the same line
                    //B on the left of C
                    //good
                    break;
                }
                
                hull.pop_back();
            }
            
            hull.push_back(C);
        }
        
        hull.pop_back();
        reverse (points.begin(), points.end());
    }
    
    return hull;
}

void polarSort (vector <P> &v, P b) {
    sort (v.begin(), v.end(), [b] (P v1, P v2) {
        //modify equal angles
        v1 -= b;
        v2 -= b;
        return make_pair (up (v1), 0ll) < make_pair (up(v2), (v1) * (v2));
    });
    
    return;
}

void solve () {
    
    return;
}

int main () {
    int t = 1;
    cin >> t;
    
    while (t--)
        solve ();
        
    return 0;
}
