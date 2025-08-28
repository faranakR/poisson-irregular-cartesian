//
// Created by Faranak Rajabi on 8/27/25.
//

#ifndef GRID2D_H
#define GRID2D_H
template<typename T>
class Grid2D {
public:
    Grid2D(T xmin, T xmax, T ymin, T ymax, int nX, int nY) : xmin(xmin), ymin(ymin), xmax(xmax), ymax(ymax), nX(nX), nY(nY) {};
    T dx = static_cast<T> ((xmax - xmin) / (nX-1));
    T dy = static_cast<T> ((ymax - ymin) / (nY-1));


private:
    T xmin, xmax, ymin, ymax;
    // T dx, dy;
    int nX, nY;
};
#endif //GRID2D_H
