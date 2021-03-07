// test!!!
#ifndef A_H_
#define A_H_
#include <vector>
#include <cmath>
namespace KDTreeSpace{};
extern std::vector<int> visit;
extern const double PI;
//for plot
struct Coordinate{
    std::vector<double>x, y;
    Coordinate() = default;
    Coordinate(std::vector<double>x_, std::vector<double>y_){
        x = x_; y = y_;
    }
    Coordinate &operator=(const Coordinate &coo){
        x = coo.x; y = coo.y;
    }
    void push(double x_, double y_){
        x.emplace_back(x_);
        y.emplace_back(y_);
    }
};
//for sort
struct point{
    double x, y, angle;
    bool operator < (const point &poi) const{
        if(x != poi.x) return x < poi.x;
        else return y < poi.y;
    }
};

Coordinate generate_all(double r0, std::vector<double>p0);
Coordinate generate_all_inside(double r0, std::vector<double>p0, double r, std::vector<double>origin);
Coordinate generate_periphery(double r0, std::vector<double>p0);
Coordinate generate_ellipse(double a, double b);
double len(double x, double y);
Coordinate FPS(const int &required_num, const Coordinate &original_coo, const Coordinate &all_coo);

//square of 3 points
double square(std::vector<double>a, std::vector<double>b, std::vector<double>c);
void sort_vec(std::vector<std::vector<double>>&points);
double cal_angle(std::vector<double>a, std::vector<double>b, std::vector<double>c);
void set_angle(const point &origin, const point &endpoint, point &a);
bool com_angle(const point &a, const point &b);
void sort_angle(const point &origin, const point &endpoint, std::vector<std::vector<double>> &points);
//flag = 0 for upper hull, 1 for lower hull
void calConvexHull(int start, int end, int flag, const std::vector<std::vector<double>> &points);
void plotConvexHull(std::vector<std::vector<double>> &hulls);

Coordinate generate_interpolation(std::vector<std::vector<double>> points);
Coordinate generate_irregular();

double sn_cal(const double &r0, const std::vector<double> &p0);
void sn_cal_all(std::vector<double> &all_square, const std::vector<double> &all_r, const std::vector<std::vector<double>> &centers);
void c_cal_all(const std::vector<double> &confidence);
std::vector<double> probabilitySize(const std::vector<double> &all_num);
double probabilityCoordinate(const double &index, const double &squ0, const std::vector<double> &p0, const std::vector<double> &p);
bool judgeNeighbour(const int &index, const std::vector<int> &mf_nearests_index, const std::vector<int> &mfs_assigned_index);

void clear();


namespace KDTreeSpace {
    template <typename T>
    class Point{
    public:
        Point(){}
        Point(Point p,int n,Point *parent =NULL){
            this->X = p.X;
            this->id = p.id;
            this->n=n;
            this->top = parent;
        }

        std::vector<T> X;
        int n = -1;//dimension

        bool    visited = false;
        double  distance = -1;

        Point *top   = NULL;
        Point *left  = NULL;
        Point *right = NULL;

        int id = -1;
    };

    template <typename T>
    class KDTree
    {
    public:
        KDTree(){}

        void feed(const std::vector<std::vector<T>> &points_t)
        {
            clear();
            if(points_t.size() <= 0)
                return;

            std::vector<Point<T>> points_tmp;
            for (size_t i =0;i<points_t.size();i++) {
                Point<T> p;
                p.X = points_t.at(i);
                p.id = i;
                points_tmp.emplace_back(p);
            }
            init(points_tmp,topPoint,0);
        }

        std::vector<int> findAdjacentK(const std::vector<T> &X,size_t k)
        {
            std::vector<Point<T>> J;
            if(topPoint == NULL){
                std::vector<int> empty_;
                return empty_;
            }
            Point<T> point_target;
            point_target.X = X;
            tryfindAdjacentAnchor(J,k,topPoint,&point_target);
            clearTempData(topPoint);
//            std::vector<std::vector<T>> points;
//            for(int i = 0; i < J.size(); i++){
//                points.push_back(J[i].X);
//            }
//            Coordinate points;
//            for(int i = 0; i < J.size(); i++){
//                points.push(J[i].X[0], J[i].X[1]);
//            }
            std::vector<int> points_id;
            for(int i = 0; i < J.size(); i++){
                points_id.emplace_back(J[i].id);
            }
            return points_id;
        }

        void clear(){
            Point<T> *anchor = topPoint;
            while (topPoint != NULL) {
                while(anchor->left != NULL){
                    anchor = anchor->left;
                }
                while(anchor->right != NULL){
                    anchor = anchor->right;
                }

                if(anchor->top == NULL){
                    delete anchor;
                    break;
                }

                anchor = anchor->top;
                if(anchor->right != NULL){
                    delete anchor->right;
                    anchor->right =NULL;
                }else{
                    delete anchor->left;
                    anchor->left =NULL;
                }
            }
            topPoint = NULL;
        }

//        //visualization
//        void toDot(const char *file_path)
//        {
//            ofstream file(file_path);
//            if(file.is_open()){
//                file<< "digraph\n{\nnode [shape = Mrecord];\n";
//                writeDot(file,topPoint);
//                file<<"}";
//                file.close();
//            }
//        }

    private:
//        void writeDot(ofstream &stream,Point<T> *p)
//        {
//            if(p == NULL){
//                return;
//            }
//            stream <<"p"<< p<<" [label = \"<f0> ("<<p->id<<")"<<p->X[0];
//            for(size_t i = 1; i < p->X.size();i++)
//            {
//                stream <<"| <f"<<i<<"> "<<p->X[i];
//            }
//            stream <<"\"];\n";
//
//            if(p->top){
//                stream<<"p"<<p->top<<":f"<<p->top->n<<" -> p"<<p<<";\n";
//            }
//            writeDot(stream,p->left);
//            writeDot(stream,p->right);
//        }

        void clearTempData(Point<T> *p_anchor)
        {
            if(p_anchor != NULL) {
                p_anchor->distance = -1;
                p_anchor->visited = false;
                clearTempData(p_anchor->left);
                clearTempData(p_anchor->right);
            }
        }
        //step1
        //down search according to p and split
        void tryfindAdjacentAnchor(std::vector<Point<T>> &J,size_t k,Point<T> *p_anchor,Point<T> *p_target)
        {
            while (true) {
                if(p_anchor->X[p_anchor->n] > p_target->X[p_anchor->n]){
                    if(p_anchor->left != NULL){
                        p_anchor = p_anchor->left;
                    }else{
                        break;
                    }
                }else{
                    if(p_anchor->right != NULL){
                        p_anchor = p_anchor->right;
                    }else {
                        break;
                    }
                }
            }
            whenIntoDown(J,k,p_anchor,p_target);
        }

        //step2
        void whenIntoDown(std::vector<Point<T>> &J,size_t k,Point<T> *p_anchor,Point<T> *p_target)
        {
            p_anchor->visited = true;
            checkWhetherReplace(J,k,p_anchor,p_target);
            climb(J,k,p_anchor,p_target);
        }

        double getMaxDistance(const std::vector<Point<T>> &J,size_t &index)
        {
            double distance = 0;
            for (size_t i = 0;i<J.size();i++) {
                if(J[i].distance > distance){
                    distance = J[i].distance;
                    index = i;
                }
            }
            return distance;
        }

        //J.size()<k,push_back the point
        //J.size()=k, if(max_distance J > p_anchor->distance), renewal the point that has max distance
        void checkWhetherReplace(std::vector<Point<T>> &J,size_t k,Point<T> *p_anchor,Point<T> *p_target)
        {
            double distance = 0;
            for(size_t i = 0;i<p_anchor->X.size();i++){
                distance += pow(p_anchor->X[i] - p_target->X[i],2);
            }

            p_anchor->distance = sqrt(distance);

            if(J.size() < k){
                J.emplace_back(*p_anchor);
            }else{
                size_t i = 0;
                double max_distance = getMaxDistance(J,i);
                if(max_distance > p_anchor->distance){
                    J[i] = *p_anchor;
                }
            }
        }

        //step3
        //climb up
        void climb(std::vector<Point<T>> &J,size_t k,Point<T> *p_anchor,Point<T> *p_target)
        {
            while (p_anchor->visited) {
                if(p_anchor == topPoint){
                    return;
                }
                p_anchor = p_anchor->top;
            }

            p_anchor->visited = true;

            //(1)
            checkWhetherReplace(J,k,p_anchor,p_target);

            //(2)
            //distance beween target and split
            double l_t_ah = p_target->X[p_anchor->n] - p_anchor->X[p_anchor->n];
            size_t index = 0;
            double max_distance = getMaxDistance(J,index);

            Point<T> *child_anchor = p_anchor->left->visited?p_anchor->right:p_anchor->left;
            if( (abs(l_t_ah) >= max_distance && J.size() >= k)||child_anchor == NULL ){
                //keep climbing
                climb(J,k,p_anchor,p_target);
            }else{
                tryfindAdjacentAnchor(J,k,
                                      child_anchor,
                                      p_target);
            }
        }

        void init(const std::vector<Point<T>> &points_t,Point<T> *parent, int n,bool is_left = true)
        {
            Point<T> *p =NULL;

            //sort in ascending order in n demensions
            std::vector<Point<T>> points_tmp(points_t);

            for(size_t i=0;i<points_tmp.size()-1;i++){
                for(size_t j=i+1;j<points_tmp.size();j++){
                    if( points_tmp[j].X.at(n) < points_tmp[i].X.at(n)){
                        Point<T> p_tmp = points_tmp[j];
                        points_tmp[j] = points_tmp[i];
                        points_tmp[i] = p_tmp;
                    }
                }
            }

            size_t anchor_index = points_tmp.size()/2;
            p = new Point<T>(points_tmp[anchor_index],n,parent);

            if(parent == NULL){
                topPoint = p;
            }else{
                if(is_left){
                    parent->left = p;
                }else{
                    parent->right = p;
                }
            }

            if(anchor_index > 0){
                std::vector<Point<T>> left;
                copy(points_tmp.begin(),points_tmp.begin()+anchor_index,back_inserter(left));
                init(left,p,(n+1)%p->X.size(),true);
            }

            if(anchor_index < (points_tmp.size() - 1)){
                std::vector<Point<T>> right;
                copy(points_tmp.begin()+anchor_index+1,points_tmp.end(),back_inserter(right));
                init(right,p,(n+1)%p->X.size(),false);
            }
        }

        Point<T> *topPoint = NULL;
    };
}

#endif //A_H_
