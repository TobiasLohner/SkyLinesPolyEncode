#include "SkyLinesPolyEncoder.h"

#include <stack>
#include <sstream>

#include <stdio.h>

using namespace std;

SkyLinesPolyEncoder::SkyLinesPolyEncoder(int numLevels, int zoomFactor, double threshold, bool forceEndpoints)
    : numLevels(numLevels), zoomFactor(zoomFactor), threshold(threshold), forceEndpoints(forceEndpoints) 
    {
        zoomLevelBreaks = new double[numLevels];
        for (int i=0; i<numLevels; i++) {
            zoomLevelBreaks[i] = threshold * pow((double)zoomFactor, numLevels - i - 1);
        }
    }

SkyLinesPolyEncoder::~SkyLinesPolyEncoder() {
    delete[] zoomLevelBreaks;
}

vector<int> SkyLinesPolyEncoder::dpEncode(vector<vector<double> >& points, char *type) {
    unsigned i, maxLoc = 0;
    stack<pair<int, int> > stack;
    double *dists = new double[points.size()];
    fill(&dists[0], &dists[points.size()], 0.0);
    double temp, maxDist;

    double absMaxDist_squared = 0.0, absMaxDist;
    double threshold_squared = pow(threshold, 2);



    // use normal douglas peucker distance (perpendicular to segment)
    // or use simple distance calculation from adjacent points
    list<size_t> points_dp, points_simple;

    for (i = 0; i < sizeof(type)/sizeof(type[0]); i++) {
      if (type[i] == 'd')
        points_simple.push_back(i);
      else if (type[i] == 'p')
        points_dp.push_back(i);
      else
        break;
    }

    for (i; i < points[0].size(); i++)
      points_dp.push_back(i);

    // simplify using Douglas-Peucker
    if (points.size() > 2) {
        stack.push(pair<unsigned, unsigned>(0, (points.size() - 1)));

        while (stack.size() > 0) {
            pair<int, int> current = stack.top();
            stack.pop();
            maxDist = 0;

            for (int i = current.first + 1; i < current.second; i++) {
                temp = max(distance_dp(points[i], points[current.first], points[current.second], points_dp),
                  distance_simple(points[i], points[current.first], points[current.second], points_simple));

                if (temp > maxDist) {
                    maxDist = temp;
                    maxLoc = i;
                }

            }

            if (maxDist > absMaxDist_squared) {
                absMaxDist_squared = maxDist;
            }

            if (maxDist > threshold_squared) {
                dists[maxLoc] = sqrt(maxDist);
                stack.push(pair<int, int>(current.first, maxLoc));
                stack.push(pair<int, int>(maxLoc, current.second));
            }
        }
    }

    absMaxDist = sqrt(absMaxDist_squared);

    vector<int> r = classify(points.size(), dists, absMaxDist);
    
    delete[] dists;
    return r;
}

/**
 * distance(p0, p1, p2) computes the distance between the point p0 and the
 * segment [p1,p2]. This could probably be replaced with something that is a
 * bit more numerically stable.
 */
double SkyLinesPolyEncoder::distance_dp(vector<double>& p0, vector<double>& p1, vector<double>& p2, list<size_t>& points) {
    double u, out = 0.0;
    double u_nom = 0.0, u_denom = 0.0;

    if (p1 == p2) {
      for (list<size_t>::iterator i = points.begin(); i != points.end(); i++) {
        out += pow(p2[*i] - p0[*i], 2);
      }
    } else {
        for (list<size_t>::iterator i = points.begin(); i != points.end(); i++) {
          u_nom += (p0[*i] - p1[*i]) * (p2[*i] - p1[*i]);
        }
        for (list<size_t>::iterator i = points.begin(); i != points.end(); i++) {
          u_denom += pow(p2[*i] - p1[*i], 2);
        }

        u = u_nom / u_denom;
        
        if (u <= 0) {
          for (list<size_t>::iterator i = points.begin(); i != points.end(); i++) {
            out += pow(p0[*i] - p1[*i], 2);
          }
        } else if (u >= 1) {
          for (list<size_t>::iterator i = points.begin(); i != points.end(); i++) {
            out += pow(p0[*i] - p2[*i], 2);
          }
        } else if (0 < u && u < 1) {
          for (list<size_t>::iterator i = points.begin(); i != points.end(); i++) {
            out += pow(p0[*i] - p1[*i] - u * (p2[*i] - p1[*i]), 2);
          }
        }
    }

    return out;
}

double SkyLinesPolyEncoder::distance_simple(vector<double>& p0, vector<double>& p1, vector<double>& p2, list<size_t>& points) {
    double out = 0.0;

    for (list<size_t>::iterator i = points.begin(); i != points.end(); i++) {
      out += sqrt(abs(p1[*i] - p0[*i])) + sqrt(abs(p2[*i] - p0[*i]));
    }

    out = pow(out, 2)/4;

    return pow(out, 2);
}

string SkyLinesPolyEncoder::encodeSignedNumber(int num) {
    int sgn_num = num << 1;
    if (num < 0) {
        sgn_num = ~(sgn_num);
    }
    return (encodeNumber(sgn_num));
}

string SkyLinesPolyEncoder::encodeNumber(int num) {
    ostringstream encodeString;

    while (num >= 0x20) {
        int nextValue = (0x20 | (num & 0x1f)) + 63;
        encodeString << ((char) (nextValue));
        num >>= 5;
    }

    num += 63;
    encodeString << ((char) (num));

    return encodeString.str();
}


auto_ptr<pair<string, string> > SkyLinesPolyEncoder::encode(vector<pair<double,double> >& points, vector<int>& levels) {
    ostringstream encodedLevels;
    ostringstream encodedPoints;

    int plat = 0;
    int plng = 0;

    size_t n_points = points.size();
    for (size_t i=0; i<n_points; i++) {
        if (levels[i] != -1) {
            encodedLevels << encodeNumber(levels[i]);
            
            pair<double, double> point = points[i];

            int late5 = floor1e5(point.second);
            int lnge5 = floor1e5(point.first);

            int dlat = late5 - plat;
            int dlng = lnge5 - plng;

            plat = late5;
            plng = lnge5;

            encodedPoints << encodeSignedNumber(dlat);
            encodedPoints << encodeSignedNumber(dlng);
        }
    }
    
    auto_ptr<pair<string, string> > r(new pair<string,string>);
    r->first = encodedPoints.str();
    r->second = encodedLevels.str();
    return r;
}

string SkyLinesPolyEncoder::encodeList(list<int>& points) {
    ostringstream encodedList;

    int val = 0;

    for (list<int>::iterator i = points.begin(); i != points.end(); i++) {
        int dval = *i - val;
        val = *i;

        encodedList << encodeSignedNumber(dval);
     }

     return encodedList.str();
}

vector<int> SkyLinesPolyEncoder::classify(size_t n_points, const double dists[], double absMaxDist) {
    vector<int> r;

    if (forceEndpoints) {
        r.push_back(numLevels - 1);
    } else {
        r.push_back(numLevels - computeLevel(absMaxDist) - 1);
    }

    if (n_points > 1) {
        for (size_t i=1; i<n_points-1; i++) {

            if (dists[i] != 0.0)
              r.push_back(numLevels - computeLevel(dists[i]) - 1);
            else
              r.push_back(-1);
        }

        if (forceEndpoints) {
            r.push_back(numLevels - 1);
        } else {
            r.push_back(numLevels - computeLevel(absMaxDist) - 1);
        }
    }

    return r;
}


/**
 * This computes the appropriate zoom level of a point in terms of it's
 * distance from the relevant segment in the DP algorithm. Could be done in
 * terms of a logarithm, but this approach makes it a bit easier to ensure
 * that the level is not too large.
 */
int SkyLinesPolyEncoder::computeLevel(double absMaxDist) {
    int lev = 0;
    if (absMaxDist > threshold) {
        while (absMaxDist < zoomLevelBreaks[lev]) {
            lev++;
        }
    }
    return lev;
}
