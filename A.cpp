#include "A.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <boost/math/distributions.hpp>
using namespace std;
const double PI = acos(-1);
const int MAX_NUM = 1000000000;

Coordinate generate_periphery(double r0, vector<double>p0){
    Coordinate temp;
    int num_points = r0*500+1;
    double stride_theta = 2 * PI / (num_points-1);
    for(int i = 0; i < num_points-1; i++){
        temp.x.emplace_back(p0[0] + r0 * cos(i * stride_theta));
        temp.y.emplace_back(p0[1] + r0 * sin(i * stride_theta));
    }
//    prevent precision problem
    temp.x.emplace_back(p0[0] + r0 * cos(0));
    temp.y.emplace_back(p0[1] + r0 * sin(0));
    return temp;
}

Coordinate generate_all(double r0, vector<double>p0){
    Coordinate temp{{p0[0]}, {p0[1]}};
    while(r0 > 0.01){
        int num_points = r0*100+1;
        double stride_theta = 2 * PI / (num_points-1);
        for(int i = 0; i < num_points; i++){
            temp.x.emplace_back(p0[0] + r0 * cos(i * stride_theta)) ;
            temp.y.emplace_back(p0[1] + r0 * sin(i * stride_theta)) ;
        }
        r0 -=0.01;
    }
    return temp;
}

Coordinate generate_all_inside(double r0, vector<double> p0, double r, vector<double> origin){
    Coordinate temp{{p0[0]}, {p0[1]}};
    while(r0 > 0.01){
        int num_points = r0*100+1;
        double stride_theta = 2 * PI / (num_points-1);
        for(int i = 0; i < num_points; i++){
            double x = p0[0] + r0 * cos(i * stride_theta);
            double y = p0[1] + r0 * sin(i * stride_theta);
            if(len(x - origin[0], y - origin[1]) > r) continue;
            temp.push(x, y);
        }
        r0 -=0.01;
    }
    return temp;
}

Coordinate generate_ellipse(double a, double b){
    Coordinate temp;
    double x = a, y = 0;
    double num = a / 0.01;
    for(int i = 0; i < num; i++){
        temp.x.emplace_back(x); temp.y.emplace_back(y);
        x -= 0.01;
        y = sqrt((1 - x*x/(a*a)) * (b*b));
    }
    x = 0; y = b;
    for(int i = 0; i < num; i++){
        temp.x.emplace_back(x); temp.y.emplace_back(y);
        x -= 0.01;
        y = sqrt((1 - x*x/(a*a)) * (b*b));
    }
    x = -a; y = 0;
    for(int i = 0; i < num; i++){
        temp.x.emplace_back(x); temp.y.emplace_back(y);
        x += 0.01;
        y = -sqrt((1 - x*x/(a*a)) * (b*b));
    }
    x = 0; y = -b;
    for(int i = 0; i < num; i++){
        temp.x.emplace_back(x); temp.y.emplace_back(y);
        x += 0.01;
        y = -sqrt((1 - x*x/(a*a)) * (b*b));
    }
    temp.x.emplace_back(a); temp.y.emplace_back(0);
    return temp;
}

double len(double x, double y){
    return pow(x*x+y*y, 0.5);
}

Coordinate FPS(const int &required_num, const Coordinate &original_coo, const Coordinate &all_coo){
    int original_length = original_coo.x.size();
    int all_length = all_coo.x.size();
    int FPS_length = required_num;
//    cout << "number of all point: " << all_length
//         << ", number of original point: " << original_length << endl;
    if(FPS_length < 0) {
        cout << "Invalid input FPS_required_num" << endl;
        exit(1);
    }
    Coordinate ans;
    vector<double> min_dis;
    for(int j = 0; j < all_length; j++){
        min_dis.emplace_back(len(all_coo.x[j] - original_coo.x[0],
                                 all_coo.y[j] - original_coo.y[0]));
    }

    //max distance in the min distance
    vector<double>::iterator max = max_element(min_dis.begin(), min_dis.end());
    auto subindex = distance(begin(min_dis), max);
//    std::cout << "Max len is " << *max<< " at index " << subindex << endl;
    ans.push(all_coo.x[subindex], all_coo.y[subindex]);
//    cout << all_coo.x_set[subindex] << " " << all_coo.y_set[subindex] << endl;

    while(--FPS_length>0){
        for(int j = 0; j < all_length; j++){
            double new_dis = len(all_coo.x[j] - *(ans.x.end() - 1),
                                 all_coo.y[j] - *(ans.y.end() - 1));
            if(new_dis < min_dis[j]) min_dis[j] = new_dis;
        }
        //max distance in the min distance
        max = max_element(min_dis.begin(), min_dis.end());
        subindex = distance(begin(min_dis), max);
//        std::cout << "Max len is " << *max<< " at index " << subindex << endl;
        ans.push(all_coo.x[subindex], all_coo.y[subindex]);
    }
    return ans;
}

//square of 3 points
double square(vector<double>a, vector<double>b, vector<double>c){
    double squ = a[0]*b[1] + c[0]*a[1] + b[0]*c[1] - c[0]*b[1] - b[0]*a[1] - a[0]*c[1];
    return squ;
}

void sort_vec(vector<vector<double>> &points){
    int len = points.size();
    point poi[len];
    for(int i=0; i<len; i++){
        poi[i].x = points[i][0];
        poi[i].y = points[i][1];
    }
    sort(poi,poi+len);
    points.clear(); vector<vector<double>>().swap(points);
    for(int i=0; i<len; i++){
        points.emplace_back(vector<double> {poi[i].x, poi[i].y});
    }
}

vector<int> visit{};
void calConvexHull(int start, int end, int flag, const vector<vector<double>> &points){
    //flag = 0 for upper hull, 1 for lower hull
    int sub = -1;
    double mul, max_squ = 0;
    switch (flag) {
        case 0 : mul = -0.5; break;
        case 1 : mul = 0.5; break;
        default: {
            cout << "Invalid input calConvexHull_flag" << endl;
            exit(1);
        }
    }
    for (int i = start+1; i < end; i++){
        double squ = mul * square(points[start], points[i], points[end]);
        if (squ > max_squ){
            max_squ = squ;
            sub = i;
        }
    }
    if (sub != -1){
        visit.emplace_back(sub);
        calConvexHull(start, sub, flag, points);
        calConvexHull(sub, end, flag, points);
    }
}

double cal_angle(vector<double>a, vector<double>b, vector<double>c){
    double numerator = (b[0]-a[0])*(c[0]-a[0]) + (b[1]-a[1])*(c[1]-a[1]);
    double denominator = pow((b[0]-a[0])*(b[0]-a[0])+(b[1]-a[1])*(b[1]-a[1]), 0.5) * pow((c[0]-a[0])*(c[0]-a[0])+(c[1]-a[1])*(c[1]-a[1]), 0.5);
    return ((c[1] < a[1]) ?  -1*acos(numerator/denominator) + 2*PI : acos(numerator/denominator));
}

void set_angle(const point &origin, const point &endpoint, point &a){
    a.angle = cal_angle(vector<double>{origin.x, origin.y},
                        vector<double>{endpoint.x, endpoint.y}, vector<double>{a.x, a.y});
}

bool com_angle(const point &a, const point &b){
    return a.angle < b.angle;
}

void sort_angle(const point &origin, const point &endpoint, vector<vector<double>> &points){
    int len = points.size();
    point poi[len];
    for(int i=0; i<len; i++){
        poi[i].x = points[i][0];
        poi[i].y = points[i][1];
        set_angle(origin, endpoint, poi[i]);
    }
    sort(poi,poi+len, com_angle);
    points.clear(); vector<vector<double>>().swap(points);
    for(int i=0; i<len; i++){
        points.emplace_back(vector<double>{poi[i].x, poi[i].y});
    }
    //notice that the first point is at end for closed plot, and z should be added ones
    points.emplace_back(vector<double>{poi[0].x, poi[0].y});
}

void plotConvexHull(vector<vector<double>> &hulls){
    sort_vec(hulls);
    point origin, endpoint;
    int sub_end = hulls.size()-1;
    origin.x = (hulls[sub_end][0]+hulls[0][0])/2; origin.y = (hulls[sub_end][1]+hulls[0][1])/2;;
    endpoint.x = origin.x + 1; endpoint.y = origin.y;
    calConvexHull(0, sub_end, 0, hulls);
    calConvexHull(0, sub_end, 1, hulls);
    vector<vector<double>> temp_hulls;
    temp_hulls.emplace_back(hulls[0]);
    for(auto i : visit){
        temp_hulls.emplace_back(hulls[i]);
    }
//    remember to clear visit, giving way to next plot
    visit.clear(); vector<int>().swap(visit);
    temp_hulls.emplace_back(*(hulls.end()-1));
//    cout << x_hull.size() << endl;
    sort_angle(origin, endpoint, temp_hulls);
//    x_hull.clear(); vector<double>().swap(x_hull);
//    y_hull.clear(); vector<double>().swap(y_hull);
//    z_hull.clear(); vector<double>().swap(z_hull);
    hulls = temp_hulls;
}

vector<double>k_line, b_line, x_l, x_r, y_l, y_r;
double  max_x = -MAX_NUM, min_x = MAX_NUM, max_y = -MAX_NUM, min_y = MAX_NUM;
Coordinate generate_interpolation(vector<vector<double>> points){
//    vector<vector<double>> contour;
    Coordinate contour;
    for(int i=1; i<points.size(); i++){
        contour.push(points[i-1][0], points[i-1][1]);
        double x1 = points[i-1][0]; double y1 = points[i-1][1];
        double x2 = points[i][0]; double y2 = points[i][1];
        if(x1 > max_x) max_x = x1; if(x1 < min_x) min_x = x1;
        if(y1 > max_y) max_y = y1; if(y1 < min_y) min_y = y1;
        if(x1 != x2){
            int num = len(x1 - x2, y1 - y2) / 0.01;
            //            if adjacent points are same as each other, continue
            if (num==0) continue;
            double delta = (x2-x1)/num;
            double k = (y1 - y2) / (x1 - x2);
            double b = (x1 * y2 - x2 * y1) / (x1 - x2);
            k_line.emplace_back(k); b_line.emplace_back(b);
            x_l.emplace_back(x1); x_r.emplace_back(x2);
//            y_l,y_r are not used
            y_l.emplace_back(y1); y_r.emplace_back(y2);
            while(num-->0){
                x1 += delta;
                contour.push(x1, k*x1+b);
            }
        }
        else{
            int num = len(0, y1 - y2) / 0.01;
            //            if adjacent points are same as each other, continue
            if (num==0) continue;
            double delta = (y2-y1)/num;
//            if x1 == x2, set k=b=0
            k_line.emplace_back(0); b_line.emplace_back(0);
            x_l.emplace_back(x1); x_r.emplace_back(x1);
            y_l.emplace_back(y1); y_r.emplace_back(y2);
            while(num-->0){
                y1 += delta;
                contour.push(x1, y1);
            }
        }
    }
    double x1 = (*(points.end()-1))[0];
    double y1 = (*(points.end()-1))[1];
    if(x1 > max_x) max_x = x1; if(x1 < min_x) min_x = x1;
    if(y1 > max_y) max_y = y1; if(y1 < min_y) min_y = y1;
    contour.x.emplace_back(x1);
    contour.y.emplace_back(y1);
    return  contour;
}

Coordinate generate_irregular(){
    Coordinate points;
    //  add pow for some precision problem on the first x and y
    double x1 = min_x+pow(10,-10); double x2 = max_x-pow(10,-10);
    double y1 = min_y+pow(10,-10); double y2 = max_y-pow(10,-10);
//    double width = *max_x - *min_x; double height = *max_y - *min_y;

    double x_temp = x1; double y_temp = y2;
    while(y_temp > y1){
        while(x_temp < x2){
            points.push(x_temp, y_temp);
            x_temp += 0.05;
        }
        x_temp = x1;
        y_temp -= 0.05;
    }
//    last row
    while(x_temp < x2){
        points.push(x_temp, y1);
        x_temp += 0.05;
    }
//    last column
    y_temp = y2;
    while(y_temp > y1){
        points.push(x2, y_temp);
        y_temp -= 0.05;
    }
//    last point
    points.push(x2, y1);

    Coordinate irregular;

//    calculate the number of points of intersection and push_back the points in  points of contour
    for(int j =0; j < points.x.size(); j++){
        int num_inter = 0;
        int vertex_count = 0;
        for(int i = 0; i < k_line.size(); i++){
            if(x_l[i] == x_r[i]){
                if(points.x[j]<=x_l[i] && points.y[j]<=max(y_l[i],y_r[i]) && points.y[j]>=min(y_l[i],y_r[i]))
                    num_inter++;
            }
//            if y1==y2
            else if(k_line[i] == 0) continue;
            else{
                double x_inter = (points.y[j] - b_line[i]) / k_line[i];
                if(x_inter <= max(x_l[i],x_r[i]) && x_inter >= min(x_l[i],x_r[i]) && x_inter >= points.x[j])
                    num_inter++;
            }
            if(points.y[j]==y_l[i] || points.y[j]==y_r[i]){
                vertex_count++;
            }
        }
        num_inter -= (vertex_count/2);
        if((num_inter%2)==1){
//            if(points.y[j]<5&&points.y[j]>3.5&&points.x[j]<-2) cout << points.x[j] << ", " << points.y[j] << endl;
            irregular.push(points.x[j], points.y[j]);
        }
    }

//    give way to next generate_interpolation
//    k_line.clear(); vector<double>().swap(k_line);
//    b_line.clear(); vector<double>().swap(b_line);
//    x_l.clear(); vector<double>().swap(x_l);
//    x_r.clear(); vector<double>().swap(x_r);
//    y_l.clear(); vector<double>().swap(y_l);
//    y_r.clear(); vector<double>().swap(y_r);
//    max_x = -MAX_NUM, min_x = MAX_NUM, max_y = -MAX_NUM, min_y = MAX_NUM;
    return irregular;
}

double sn_cal(const double &r0, const vector<double> &p0){
    Coordinate points = generate_all(r0, p0);
    int actual_square = 0;

//    calculate the number of points of intersection and push_back the points in  points of contour
    for(int j =0; j < points.x.size(); j++){
        int num_inter = 0;
        int vertex_count = 0;
        for(int i = 0; i < k_line.size(); i++){
            if(x_l[i] == x_r[i]){
                if(points.x[j]<=x_l[i] && points.y[j]<=max(y_l[i],y_r[i]) && points.y[j]>=min(y_l[i],y_r[i]))
                    num_inter++;
            }
//            if y1==y2
            else if(k_line[i] == 0) continue;
            else{
                double x_inter = (points.y[j] - b_line[i]) / k_line[i];
                if(x_inter <= max(x_l[i],x_r[i]) && x_inter >= min(x_l[i],x_r[i]) && x_inter >= points.x[j])
                    num_inter++;
            }
            if(points.y[j]==y_l[i] || points.y[j]==y_r[i]){
                vertex_count++;
            }
        }
        num_inter -= (vertex_count/2);
        if((num_inter%2)==1){
//            if(points.y[j]<5&&points.y[j]>3.5&&points.x[j]<-2) cout << points.x[j] << ", " << points.y[j] << endl;
            actual_square++;
        }
    }
    double sn = (1.0 * actual_square) / points.x.size();
    return sn;
}

vector<double> sn_all;
void sn_cal_all(vector<double> &all_square, const vector<double> &all_r, const vector<vector<double>> &centers){
    for(int i = 0; i < all_r.size(); i++){
        double sn = sn_cal(all_r[i], centers[i]);
        sn_all.emplace_back(sn);
        all_square[i] *= sn;
    }
}

vector<double> c_all;
void c_cal_all(const vector<double> &confidence){
    boost::math::chi_squared dist(2);
    for(int i = 0; i < sn_all.size(); i++){
//        double pro = boost::math::cdf(dist,9.21);
        if(sn_all[i] > 0.999) sn_all[i] = 0.999;
//        if(sn_all[i] <=0 || sn_all[i]>=1) cout << "sn!!! " << sn_all[i] << " " << i << endl;
        double c = boost::math::quantile(dist, confidence[i]);
        c_all.push_back(c);
    }
}

vector<double> probabilitySize(const vector<double> &all_num){
    vector<double> temp = all_num;
    double sum = accumulate(temp.begin(), temp.end(),0);
    for(int i=0; i<temp.size(); i++)
        temp[i] /= sum;
    return temp;
}

double probabilityCoordinate(const double &index, const double &squ0, const vector<double> &p0, const vector<double> &p){
//    double integral of Gaussian distrinution
//    double sn = 2 * PI * sigma * (1 - exp(-1*len(p0[0]-p[0],p0[1]-p[1])/2/sigma));
    double sigma = squ0 /(PI * c_all[index]);
    double temp = -1 * (pow(p0[0]-p[0], 2) + pow(p0[1]-p[1], 2)) / (2 * pow(sigma, 2));
    double pro = 1 / (sn_all[index] * 2 * PI * sigma) * exp(temp);
//    attention, pro may larger than MAX_DOUBLE
    if(pro > 1000000) pro = 1000000;
    return pro;
}

bool judgeNeighbour(const int &index, const vector<int> &mf_nearests_index, const vector<int> &mfs_assigned_index){
    for(auto i = mf_nearests_index.begin(); i != mf_nearests_index.end(); i++){
        if((*i) != -1) {
            if(mfs_assigned_index[*i] == index) return true;
        }
    }
    return false;
}

void clear(){
    //    give way to next generate_interpolation
    k_line.clear(); vector<double>().swap(k_line);
    b_line.clear(); vector<double>().swap(b_line);
    x_l.clear(); vector<double>().swap(x_l);
    x_r.clear(); vector<double>().swap(x_r);
    y_l.clear(); vector<double>().swap(y_l);
    y_r.clear(); vector<double>().swap(y_r);
    max_x = -MAX_NUM, min_x = MAX_NUM, max_y = -MAX_NUM, min_y = MAX_NUM;
    //   give way to next sn_cal_all
    sn_all.clear(); vector<double>().swap(sn_all);
    //   give way to next c_cal_all
    c_all.clear(); vector<double>().swap(c_all);
}
