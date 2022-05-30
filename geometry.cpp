/*
 *  never forget collinearity
 */


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
const double pi = acos(-1.0);
const double eps = 1e-9;

const ll MX = 2e5 + 9;
const ll inf = 1e18;

struct P {
    ll x, y;
    
    void read() {
        cin >> x >> y;
    }
    
    void Write() {
        cout << "(x, y) " << x << " " << y << "\n";
    }
    
    P (): x(0), y(0) {}
    
    P (ll a, ll b) : x(a), y(b) {}
    
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
    
    ll triangle (const P &b, const P& c) {
        return (b - *this) * (c - *this);
    }
    
    bool operator < (const P& b) {
        return make_pair(x, y) < make_pair(b.x, b.y);
    }
};

typedef pair <P, P> Segment;

ll Sign (P& p, Segment seg) {
    ll x = seg.F.triangle (p, seg.S);
    return x == 0ll ? 0ll : (x / abs (x));
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

bool PointInSegment (P& p, Segment seg) {
    if (p.triangle (seg.F, seg.S) != 0)
        return false;
        
    // p is somewhere inside the bounding box of the segment
    return min (seg.F.x, seg.S.x) <= p.x && p.x <= max (seg.F.x, seg.S.x)
        && min (seg.F.y, seg.S.y) <= p.y && p.y <= max (seg.F.y, seg.S.y);
}

int PointInPolygon (P& p, vector <P>& Polygon) {
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
    
    return (cnt % 2) ? 1 : -1;
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

bool up (P p) {
    assert (p.x != 0 || p.y != 0); //argument of 0 0 is undefined
    
    return p.y > 0 || (p.y == 0 && p.x < 0);
}

void polarSort (vector <P> &v) {
    sort (v.begin(), v.end(), [] (P v1, P v2) {
        //modify equal angles
        return make_pair (up (v1), 0ll) < make_pair (up(v2), v1 * v2);
    });
    
    return;
}

void solve () {
    
    
    return;
}

int main () {
    int t = 1;
    //cin >> t;
    
    while (t--)
        solve ();
        
    return 0;
}
