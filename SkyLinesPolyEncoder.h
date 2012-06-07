#ifndef SKYLINESPOLYENCODER_H_
#define SKYLINESPOLYENCODER_H_

#include <vector>
#include <list>
#include <string>
#include <cmath>
#include <memory>

class SkyLinesPolyEncoder {

private:
    int numLevels;
    int zoomFactor;
    double threshold;
    bool forceEndpoints;
    double *zoomLevelBreaks;
    
    void _buildZoomLevelBreaks();
    double distance_dp(std::vector<double>& p0, std::vector<double>& p1, std::vector<double>& p2, std::list<size_t>& points);
    double distance_simple(std::vector<double>& p0, std::vector<double>& p1, std::vector<double>& p2, std::list<size_t>& points);
    inline int floor1e5(double coordinate) { return (int)floor(coordinate * 1e5); }
    std::string encodeSignedNumber(int num);
    std::string encodeNumber(int num);
    int computeLevel(double absMaxDist);
    std::vector<int> classify(size_t n_points, const double dists[], double absMaxDist);
    
public:
    SkyLinesPolyEncoder(int numLevels=18, int zoomFactor=2, double threshold=0.00001, bool forceEndpoints=true);
    ~SkyLinesPolyEncoder();

    std::vector<int> dpEncode(std::vector<std::vector<double> >& points, char *type);
    std::auto_ptr<std::pair<std::string, std::string> > encode(std::vector<std::pair<double,double> >& points, std::vector<int>& levels);
    std::string encodeList(std::list<int>& points);
    int getNumLevels() { return numLevels; }
    int getZoomFactor() { return zoomFactor; }
};

#endif /*SKYLINESPOLYENCODER_H_*/
